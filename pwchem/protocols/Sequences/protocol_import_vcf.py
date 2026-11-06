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
import zlib

from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import VCFFile, SetOfVCFFiles
from pwchem.utils.utilsRNA import GATKPreparationMixin


class ProtImportVCF(EMProtocol):
    """
    Import one or more known-variant VCF files.

    This protocol imports local Variant Call Format (VCF) files into
    Scipion and creates a ``SetOfVCFFiles`` that can be used as input for
    downstream genomic and variant-processing protocols.

    Multiple VCF files can be imported in a single protocol execution.
    Each imported file is represented by an independent ``VCFFile`` object
    containing the VCF path and the available variant metadata.

    Existing indexes (.tbi, .csi, .idx) can be imported. If no index is
    specified, the protocol looks for a neighboring index and, if needed,
    creates one for supported VCF formats.


    ---------------------------------------------------------------------
    Input parameters
    ---------------------------------------------------------------------

    vcfsData : str
        JSON list describing the VCF files to import. At least one entry
        is required, and every entry must be an object. Metadata and paths
        must be strings; optional fields can be omitted or set to null.
        Relative paths are resolved against the execution working directory.

        Each entry can contain the following fields:

        scientificName
            Scientific name of the organism associated with the variants.

            For example::

                Homo sapiens
                Mus musculus

        assembly
            Genome assembly associated with the VCF.

            For example::

                GRCh38
                GRCm39

        vcfFile
            Path to the local VCF file to import.

            This field is required. Supported extensions are ``.vcf`` and
            ``.vcf.gz``. Compressed inputs must use BGZF, as required for
            genomic indexing; ordinary gzip is rejected.

        indexFile
            Optional path to an existing VCF index.

            BGZF VCFs accept Tabix (``.tbi``) or CSI (``.csi``) indexes.
            Plain VCFs accept Tribble (``.idx``) indexes. Other extensions
            and incompatible VCF/index combinations are rejected.

        The ``vcfsData`` parameter is normally configured through the
        corresponding Scipion wizard rather than edited manually.


    ---------------------------------------------------------------------
    Imported VCF
    ---------------------------------------------------------------------

    Each configured VCF file is copied into the protocol output directory.

    Imported files are renamed using the scientific name and, when
    available, the genome assembly. Metadata are trimmed, and characters
    unsuitable for filenames, including path separators, are replaced by
    underscores. Missing names use ``variants``. Numeric suffixes keep
    entries with identical or equivalent names from overwriting each other.

    For example, a VCF configured as::

        scientificName = Mus musculus
        assembly = GRCm39

    is stored using a name based on::

        Mus_musculus_GRCm39_variants

    The original VCF extension is preserved. Compressed ``.vcf.gz`` files
    therefore remain compressed after import. The original VCF and index
    files are read without modification; downstream objects reference the
    local copies rather than depending on the source paths.

    Each imported VCF is represented by a ``VCFFile`` object with the
    following metadata:

    scientificName
        Scientific name provided during import, when available.

    assembly
        Genome assembly provided during import, when available.

    source
        Set to ``Imported``.

    database
        Set to ``Local``.

    variantType
        Set to ``known``.

    isCompressed
        Determined from whether the imported VCF filename ends in ``.gz``.


    ---------------------------------------------------------------------
    VCF index
    ---------------------------------------------------------------------

    An existing VCF index can optionally be imported together with the
    corresponding VCF.

    For Tabix indexes (``.tbi``), the imported index is stored as::

        <imported_vcf>.tbi

    For CSI indexes (``.csi``), the imported index is stored as::

        <imported_vcf>.csi

    Tribble indexes (``.idx``) are stored as ``<imported_vcf>.idx``
    beside a plain VCF. Index types are not converted during copying.

    The resulting index path is stored in the corresponding ``VCFFile``
    object.

    A supplied index must exist and be non-empty. A compatible index is
    reused only when it is at least as recent as the source VCF and its
    binary signature matches the expected index format. An unusable supplied
    index generates a warning, then neighboring compatible indexes are tried.
    For BGZF VCFs, Tabix is preferred to CSI when both are usable.

    If no usable index is available, the protocol generates a .tbi index for
    BGZF-compressed .vcf.gz using bcftools, or an .idx index for plain .vcf
    using GATK. Previously imported local indexes are removed before copying
    or generating the new index, so they cannot mask a failed execution or
    remain associated with a replaced VCF. Generated signatures are checked
    before publication. The input VCF is never modified.

    Signature and timestamp checks do not inspect every index entry or prove
    that a supplied index belongs to this exact VCF. The user must provide
    correctly paired source files. In particular, a CSI index remains CSI;
    GATK preprocessing can generate a local Tribble/Tabix index when it
    requires a different index type.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputVCFs : SetOfVCFFiles
        Set containing all imported VCF files.

        A ``SetOfVCFFiles`` is created independently of the number of VCF
        files configured for import.

        Each ``VCFFile`` contains:

        - Path to the copied VCF file.
        - Scientific name, when provided.
        - Genome assembly, when provided.
        - Source set to ``Imported``.
        - Database set to ``Local``.
        - Variant type set to ``known``.
        - Compression status.
        - VCF index path (imported or generated).

        Each object is labelled using the scientific name and genome
        assembly when available.


    ---------------------------------------------------------------------
    Workflow
    ---------------------------------------------------------------------

    The protocol performs the following steps for each configured VCF:

    1. Read the VCF configuration from ``vcfsData``.

    2. Generate an output filename using the scientific name and genome
       assembly when available.

    3. Copy the VCF file into the protocol output directory.

    4. Create a ``VCFFile`` object and assign the available metadata.

    5. Copy an existing index (explicit or adjacent) or generate one,
       then associate it with the ``VCFFile``.

    6. Add the resulting ``VCFFile`` to the output ``SetOfVCFFiles``.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - The VCF configuration is a JSON list of objects with valid field types.
    - At least one VCF has been configured.
    - Every configured entry contains a VCF file.
    - Every VCF file exists, is non-empty and has a supported extension.
    - Compressed VCFs have a valid initial BGZF block.
    - A VCF version declaration and valid column header are present.
    - Declared contig names and sample names are not duplicated.
    - Every specified index exists, is non-empty and has a compatible suffix.

    Scientific name, genome assembly and VCF index are optional.

    Header checks reuse the same routines as GATK input preparation. They
    accept sites-only VCFs and valid headers without variant records; they
    do not scan all variants, verify REF alleles against a FASTA or establish
    assembly compatibility. Index generation retains the external tools'
    validation, including ordering requirements. No reference FASTA is
    required for this import protocol.


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Multiple VCF files can be imported in a single protocol execution.

    - Input VCF files are copied into the protocol output directory.

    - The original VCF files are not modified.

    - Both BGZF-compressed and uncompressed VCF files can be imported.

    - ``.vcf.gz`` extensions are preserved when generating the imported
      filename.

    - Imported VCF files are classified as known variants. This metadata
      label does not certify their suitability as BQSR known sites, a
      Mutect2 germline population resource or a panel of normals. Those
      uses require appropriately prepared biological resources.

    - The source of imported VCF objects is stored as ``Imported`` and
      their database as ``Local``.

    - Scientific name and genome assembly metadata are optional.

    - VCF index files are optional inputs; missing indexes are generated.

    - Existing ``.tbi`` and ``.csi`` indexes are associated directly with
      the copied VCF using the standard VCF index naming convention.

    - Index generation requires bcftools (.vcf.gz) or GATK (.vcf) in
      the RNASEQ_DIC environment. BGZF VCFs must be sorted.

    - The protocol always returns a ``SetOfVCFFiles`` to provide a
      consistent output type for downstream Scipion protocols.
    """

    VCF_GZ_EXTENSION = '.vcf.gz'

    _label = 'import vcf'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'vcfsData',
            params.StringParam,
            default='[]',
            label='VCFs',
            help='Configure the VCF files to import.'
        )

    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    def _getVCFsData(self):
        try:
            entries = json.loads(self.vcfsData.get())
        except (TypeError, ValueError) as error:
            raise ValueError('Invalid VCF configuration: expected a JSON list.') from error
        if not isinstance(entries, list):
            raise ValueError('VCF configuration must be a JSON list.')
        if not entries:
            raise ValueError('At least one VCF must be configured.')

        return [self._normalizeVCFEntry(entry, index) for index, entry in enumerate(entries, start=1)]

    @classmethod
    def _normalizeVCFEntry(cls, entry, index):
        if not isinstance(entry, dict):
            raise ValueError('VCF {} must be a JSON object.'.format(index))
        fields = ('scientificName', 'assembly', 'vcfFile', 'indexFile')
        return {field: cls._normalizeVCFValue(entry.get(field), field, index) for field in fields}

    @staticmethod
    def _normalizeVCFValue(value, field, index):
        if value is None:
            value = ''
        if not isinstance(value, str) or '\x00' in value:
            raise ValueError('VCF {}: {} must be text.'.format(index, field))
        if field.endswith('File'):
            return os.path.abspath(value) if value else ''
        return value.strip()

    @staticmethod
    def _getVCFStem(vcfData):
        parts = [vcfData['scientificName'] or 'variants', vcfData['assembly']]
        name = '_'.join(part for part in parts if part)
        return re.sub(r'[^\w.-]+', '_', name).strip('._-') or 'variants'

    def _createVCF(self, vcfData, safeName):
        source = vcfData['vcfFile']
        extension = self.VCF_GZ_EXTENSION if source.endswith(self.VCF_GZ_EXTENSION) else '.vcf'
        importedVCF = self._getExtraPath(safeName + '_variants' + extension)
        copyFile(source, importedVCF)
        importedIndex = self._ensureVCFIndex(source, importedVCF, vcfData['indexFile'])

        vcf = VCFFile(filename=importedVCF)
        if vcfData['scientificName']:
            vcf.setScientificName(vcfData['scientificName'])
        if vcfData['assembly']:
            vcf.setAssembly(vcfData['assembly'])
        vcf.setSource('Imported')
        vcf.setDatabase('Local')
        vcf.setVariantType('known')
        vcf.setIsCompressed(extension == self.VCF_GZ_EXTENSION)
        vcf.setIndexFile(importedIndex)
        label = vcfData['scientificName'] or 'VCF'
        if vcfData['assembly']:
            label += ' ({})'.format(vcfData['assembly'])
        vcf.setObjLabel(label)
        return vcf

    @staticmethod
    def _indexSuffixes(vcfFile):
        return ('.tbi', '.csi') if vcfFile.endswith('.gz') else ('.idx',)

    def _isUsableIndex(self, indexFile, vcfFile):
        suffix = os.path.splitext(indexFile)[1].lower()
        signature = {'.tbi': b'TBI\x01', '.csi': b'CSI\x01', '.idx': b'TIDX'}.get(suffix)
        return (suffix in self._indexSuffixes(vcfFile)
                and GATKPreparationMixin.isCurrentIndex(indexFile, vcfFile)
                and GATKPreparationMixin.hasIndexSignature(indexFile, signature))

    def _findAdjacentVCFIndex(self, vcfFile):
        for suffix in self._indexSuffixes(vcfFile):
            candidate = vcfFile + suffix
            if self._isUsableIndex(candidate, vcfFile):
                return candidate
        return None

    def _ensureVCFIndex(self, sourceVCF, importedVCF, indexFile):
        """Copy a usable source index or generate one for the current local VCF."""
        for suffix in ('.tbi', '.csi', '.idx'):
            path = importedVCF + suffix
            if os.path.lexists(path):
                os.remove(path)
        if indexFile and not self._isUsableIndex(indexFile, sourceVCF):
            self.warning('Ignoring a stale or invalid VCF index: {}'.format(indexFile))
            indexFile = ''
        if not indexFile:
            indexFile = self._findAdjacentVCFIndex(sourceVCF)
        if indexFile:
            importedIndex = importedVCF + os.path.splitext(indexFile)[1].lower()
            copyFile(indexFile, importedIndex)
            return importedIndex

        if importedVCF.endswith(self.VCF_GZ_EXTENSION):
            args = 'index -f -t ' + shlex.quote(importedVCF)
            expectedIndex = importedVCF + '.tbi'
            program = 'bcftools'
        else:
            args = 'IndexFeatureFile -I ' + shlex.quote(importedVCF)
            expectedIndex = importedVCF + '.idx'
            program = 'gatk'
        Plugin.runCondaCommand(self, args, RNASEQ_DIC, program)
        if not self._isUsableIndex(expectedIndex, importedVCF):
            raise RuntimeError('A valid index was not generated for {} (expected {}).'.format(
                importedVCF, expectedIndex))
        return expectedIndex

    def createOutputStep(self):
        try:
            vcfsData = self._getVCFsData()
        except ValueError as error:
            raise RuntimeError(str(error)) from error
        errors = self._validateFiles(vcfsData)
        if errors:
            raise RuntimeError('Invalid VCF inputs: ' + '; '.join(errors))

        os.makedirs(self._getExtraPath(), exist_ok=True)
        outputVCFs = SetOfVCFFiles().create(outputPath=self._getPath())
        outputVCFs.setObjLabel('Known variants VCFs')
        usedNames = set()
        try:
            for vcfData in vcfsData:
                stem = self._getVCFStem(vcfData)
                safeName, suffix = stem, 2
                while safeName.casefold() in usedNames:
                    safeName = '{}_{}'.format(stem, suffix)
                    suffix += 1
                usedNames.add(safeName.casefold())
                outputVCFs.append(self._createVCF(vcfData, safeName))
            self._defineOutputs(outputVCFs=outputVCFs)
        except Exception:
            outputVCFs.close()
            raise

    def _validateFiles(self, vcfsData):
        errors = []
        for index, vcfData in enumerate(vcfsData, start=1):
            prefix = 'VCF {}: '.format(index)
            path = vcfData['vcfFile']
            if not path:
                errors.append(prefix + 'VCF file is required.')
                continue
            supported = path.endswith(('.vcf', self.VCF_GZ_EXTENSION))
            if not supported:
                errors.append(prefix + 'VCF must use .vcf or .vcf.gz extension.')
            for field, label in (('vcfFile', 'VCF'), ('indexFile', 'index')):
                errors.extend(self._validateVcfFile(vcfData[field], field, label, path, supported, prefix))
        return errors

    def _validateVcfFile(self, filename, field, label, path, supported, prefix):
        if not filename:
            return []
        errors = []
        try:
            if not os.path.isfile(filename):
                errors.append('{}{} file does not exist: {}'.format(prefix, label, filename))
            elif os.path.getsize(filename) == 0:
                errors.append('{}{} file is empty: {}'.format(prefix, label, filename))
            elif field == 'vcfFile' and supported:
                GATKPreparationMixin.readVCFContigs(filename)
            elif field == 'indexFile':
                if os.path.splitext(filename)[1].lower() not in self._indexSuffixes(path):
                    errors.append(prefix + 'Index extension is incompatible with the VCF format.')
        except (OSError, EOFError, UnicodeError, RuntimeError, zlib.error) as error:
            errors.append(prefix + str(error))
        return errors

    def _validate(self):
        try:
            return self._validateFiles(self._getVCFsData())
        except ValueError as error:
            return [str(error)]
