# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import json
import os
import re
import shlex

from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import Genome, SetOfGenomes
from pwchem.utils.utilsRNA import GATKPreparationMixin


class ProtImportGenomes(EMProtocol):
    """Import local reference FASTA files and optional GTF annotations.

    Each configured entry creates a ``Genome`` in the output ``SetOfGenomes``.
    Files are copied into the protocol extra directory, so subsequent protocols
    do not depend on the original input paths. Source files and their indexes
    are not modified. A samtools FASTA index is generated for each imported
    reference; sequence dictionaries and aligner indexes are prepared by the
    downstream protocols that require them.

    Input parameters
    ----------------
    genomesData : str
        JSON list configured through the genome-import wizard. At least one
        entry is required. Each entry is an object with the following fields:

        - ``fastaFile``: required path to an existing, non-empty, uncompressed
          reference FASTA. Its first line must contain a non-empty FASTA header.
        - ``gtfFile``: optional path to an existing, non-empty, uncompressed
          GTF annotation. Omit it or use an empty value when no annotation is
          available. GTF content and assembly compatibility are not validated.
        - ``scientificName``: optional scientific name, such as ``Mus musculus``.
        - ``assembly``: optional reference assembly, such as ``GRCm39``.
        - ``release``: optional annotation release. Text or an integer release
          is accepted and stored as text.

        Relative input paths are resolved against the execution working
        directory. Metadata values are trimmed; absent optional fields remain
        unset. Additional fields are ignored. FASTA and GTF files compressed
        with gzip or BGZF must be decompressed before import, because the
        downstream RNA-seq protocols expect uncompressed reference files.

    File names and indexing
    -----------------------
    Scientific name and assembly are used to form readable output filenames.
    Spaces and path separators are replaced by underscores; empty names use
    ``genome``. If several entries produce the same name, numeric suffixes
    keep their files separate. Labels retain the supplied metadata.

    An imported reference normally has a name such as
    ``Mus_musculus_GRCm39_genome.fa``. FASTA extensions ``.fa``, ``.fasta`` and
    ``.fna`` are retained; other input suffixes use ``.fa``. An annotation is
    stored as ``Mus_musculus_GRCm39_annotation.gtf``. The FASTA index is stored
    beside its reference with the additional ``.fai`` suffix.

    The local FASTA index is regenerated after copying the reference, including
    on a repeated execution. Existing source indexes are not reused. The
    generated index must be non-empty and have valid sequence names, lengths,
    line widths and offsets within the copied FASTA. These checks do not
    validate every reference base or establish biological reference identity.

    Output
    ------
    referenceGenomes : SetOfGenomes
        Imported references in configuration order. Each ``Genome`` stores
        its local FASTA path, optional GTF path, scientific name, assembly,
        release and the source ``Imported``. Entries with identical metadata
        remain separate objects with distinct file paths.

    Validation and requirements
    ---------------------------
    The protocol checks the JSON structure, supported field types, file
    existence, non-empty contents and the FASTA header before copying files.
    Invalid configurations are reported before a result is published. FASTA
    sequence parsing and indexing are performed by samtools; annotation format
    and compatibility with the FASTA remain the user's responsibility.

    samtools must be available in the RNA-seq environment configured by the
    plugin. The protocol directory must be writable and have enough space for
    copies of all references and annotations. Input directories require only
    read access; original files can be removed after a successful import.
    """


    _label = 'import genomes'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'genomesData',
            params.StringParam,
            default='[]',
            label='Genomes',
            help='Configure the genomes to import.'
        )

    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    def _getGenomesData(self):
        try:
            entries = json.loads(self.genomesData.get())
        except (TypeError, ValueError) as error:
            raise ValueError('Invalid genomes configuration: expected a JSON list.') from error
        if not isinstance(entries, list):
            raise ValueError('Genomes configuration must be a JSON list.')
        if not entries:
            raise ValueError('At least one genome must be configured.')

        return [self._normalizeGenomeEntry(entry, index) for index, entry in enumerate(entries, start=1)]

    @classmethod
    def _normalizeGenomeEntry(cls, entry, index):
        if not isinstance(entry, dict):
            raise ValueError('Genome {} must be a JSON object.'.format(index))
        fields = ('scientificName', 'assembly', 'release', 'fastaFile', 'gtfFile')
        return {field: cls._normalizeGenomeValue(entry.get(field), field, index) for field in fields}

    @staticmethod
    def _normalizeGenomeValue(value, field, index):
        if value is None:
            value = ''
        if field == 'release' and type(value) is int:
            value = str(value)
        if not isinstance(value, str) or '\x00' in value:
            raise ValueError('Genome {}: {} must be text{}.'.format(
                index, field, ' or an integer' if field == 'release' else ''))
        if field.endswith('File'):
            return os.path.abspath(value) if value else ''
        return value.strip()

    @staticmethod
    def _getGenomeStem(genomeData):
        parts = [genomeData['scientificName'], genomeData['assembly']]
        name = '_'.join(part for part in parts if part)
        return re.sub(r'[^\w.-]+', '_', name).strip('._-') or 'genome'

    def _createGenome(self, genomeData, safeName):
        genome = Genome()
        for field, setter in (('scientificName', genome.setScientificName),
                              ('assembly', genome.setAssembly), ('release', genome.setRelease)):
            if genomeData[field]:
                setter(genomeData[field])
        genome.setSource('Imported')

        fastaFile = genomeData['fastaFile']
        extension = os.path.splitext(fastaFile)[1]
        if extension.lower() not in ('.fa', '.fasta', '.fna'):
            extension = '.fa'
        importedFasta = self._getExtraPath(safeName + '_genome' + extension)
        copyFile(fastaFile, importedFasta)
        self._createFastaIndex(importedFasta)
        genome.setFastaFile(importedFasta)

        if genomeData['gtfFile']:
            importedGtf = self._getExtraPath(safeName + '_annotation.gtf')
            copyFile(genomeData['gtfFile'], importedGtf)
            genome.setGtfFile(importedGtf)

        label = genomeData['scientificName'] or 'Genome'
        if genomeData['assembly']:
            label += ' ({})'.format(genomeData['assembly'])
        genome.setObjLabel(label)
        return genome

    def createOutputStep(self):
        try:
            genomesData = self._getGenomesData()
        except ValueError as error:
            raise RuntimeError(str(error)) from error
        errors = self._validateFiles(genomesData)
        if errors:
            raise RuntimeError('Invalid genome inputs: ' + '; '.join(errors))

        os.makedirs(self._getExtraPath(), exist_ok=True)
        referenceGenomes = SetOfGenomes().create(outputPath=self._getPath())
        referenceGenomes.setObjLabel('Reference genomes')
        usedNames = set()
        try:
            for genomeData in genomesData:
                stem = self._getGenomeStem(genomeData)
                safeName, suffix = stem, 2
                while safeName.casefold() in usedNames:
                    safeName = '{}_{}'.format(stem, suffix)
                    suffix += 1
                usedNames.add(safeName.casefold())
                referenceGenomes.append(self._createGenome(genomeData, safeName))
            self._defineOutputs(referenceGenomes=referenceGenomes)
        except Exception:
            referenceGenomes.close()
            raise

    def _createFastaIndex(self, fastaFile):
        """Build an index for the copied FASTA instead of reusing an older one."""
        fastaIndex = fastaFile + '.fai'
        if os.path.lexists(fastaIndex):
            os.remove(fastaIndex)
        Plugin.runCondaCommand(self, 'faidx ' + shlex.quote(fastaFile), RNASEQ_DIC, 'samtools')
        try:
            GATKPreparationMixin.readFastaIndex(fastaIndex)
        except (OSError, ValueError) as error:
            raise RuntimeError('samtools did not create a valid FASTA index: {}'.format(fastaIndex)) from error

    @classmethod
    def _validateFiles(cls, genomesData):
        errors = []
        for index, genomeData in enumerate(genomesData, start=1):
            for field, label in (('fastaFile', 'FASTA'), ('gtfFile', 'GTF')):
                prefix = 'Genome {}: {}'.format(index, label)
                errors.extend(cls._validateGenomeFile(genomeData[field], prefix, field == 'fastaFile'))
        return errors

    @classmethod
    def _validateGenomeFile(cls, path, prefix, isFasta):
        if not path:
            return [prefix + ' file is required.'] if isFasta else []
        try:
            if not os.path.isfile(path):
                return ['{} file does not exist: {}'.format(prefix, path)]
            if os.path.getsize(path) == 0:
                return ['{} file is empty: {}'.format(prefix, path)]
            return cls._validateGenomeHeader(path, prefix, isFasta)
        except OSError as error:
            return ['{} file could not be read: {}'.format(prefix, error)]

    @staticmethod
    def _validateGenomeHeader(path, prefix, isFasta):
        with open(path, 'rb') as handle:
            signature = handle.read(2)
            if signature == b'\x1f\x8b' or path.lower().endswith('.gz'):
                return ['{} must be uncompressed: {}'.format(prefix, path)]
            if isFasta:
                handle.seek(0)
                header = handle.readline().rstrip(b'\r\n')
                if not header.startswith(b'>') or not header[1:].strip():
                    return ['{} has an invalid FASTA header: {}'.format(prefix, path)]
        return []

    def _validate(self):
        try:
            return self._validateFiles(self._getGenomesData())
        except ValueError as error:
            return [str(error)]
