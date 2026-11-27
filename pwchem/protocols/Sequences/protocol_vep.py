# **************************************************************************
# *
# * Authors:    Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import os
import re
import shutil
import tempfile

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import BooleanParam, IntParam, PathParam, PointerParam, StringParam

from pwchem.objects import VariantAnnotation
from pwchem.utils.utilsRNA import GATKPreparationMixin, iterVepRows


class ProtVEP(GATKPreparationMixin, EMProtocol):
    """Annotate a VCF offline with Ensembl VEP and select missense effects.

    Implements the reference VEP command with cache, UniProt, protein,
    symbol, transcript-version and TSL annotations. The output remains the
    native tab-delimited VEP format, including its comment/header lines.
    By default filter_vep selects ``Consequence is missense_variant``;
    outputAllAnnotations always retains the complete annotation. Disable
    selectMissense to pass all consequences to outputAnnotations.

    Supply an installed Ensembl cache root, species, assembly and cache
    version. cacheDirectory must contain species/version_assembly/info.txt.
    Human and mouse use homo_sapiens/GRCh38 and mus_musculus/GRCm39,
    respectively. Species and assembly must agree with VCF metadata when
    available. Cache resources are not downloaded by this protocol. Use a
    cache release compatible with the installed VEP (116 by default).
    A cache made for another genome cannot be substituted by renaming it.

    Existing VCF filters are retained as annotations: this step does not
    discard non-PASS calls, known variants or germline variants. Filtering
    somatic calls belongs upstream. Empty result tables with a valid VEP
    header are legitimate. Variants may have several transcript effects;
    one annotation row is not necessarily one unique genomic mutation.

    outputAnnotations and outputAllAnnotations are VariantAnnotation
    objects containing VEP tables and species/assembly/cache metadata.
    They connect to the mutation/peptide protocol without altering headers
    with sed. Local output files are replaced only after successful checks.
    Requires ensembl-vep and filter_vep in the RNA-seq environment.
    """

    _label = 'VEP annotation / missense selection'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVCF', PointerParam, pointerClass='VCFFile', label='Variants:')
        form.addParam('cacheDirectory', PathParam, label='VEP cache root:')
        form.addParam('species', StringParam, default='homo_sapiens', label='Ensembl species:')
        form.addParam('assembly', StringParam, default='GRCh38', label='Genome assembly:')
        form.addParam('cacheVersion', IntParam, default=116, label='Cache release:')
        form.addParam('selectMissense', BooleanParam, default=True, label='Select missense effects:')
        form.addParallelSection(threads=4, mpi=1)

    def _getCommandFile(self):
        return self._getExtraPath('vep_commands.txt')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.annotateStep)
        self._insertFunctionStep(self.createOutputStep)

    def _validate(self):
        errors = []
        obj = self.inputVCF.get()
        self._validateVcf(obj, errors)
        species, assembly = self.species.get() or '', self.assembly.get() or ''
        if not re.fullmatch(r'[a-z][a-z0-9_]+', species) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', assembly):
            errors.append('Use an Ensembl species name and a valid assembly name.')
        cache = os.path.join(self.cacheDirectory.get() or '', species,
                             '{}_{}'.format(self.cacheVersion.get(), assembly), 'info.txt')
        if not os.path.isfile(cache):
            errors.append('VEP cache info.txt not found: {}'.format(cache))
        if self.cacheVersion.get() < 1 or self.numberOfThreads.get() < 1:
            errors.append('Cache release and thread count must be positive.')
        return errors

    def _validateVcf(self, obj, errors):
        if obj is None or not obj.getFileName() or not os.path.isfile(obj.getFileName()):
            errors.append('An existing input VCF is required.')
        else:
            try:
                self.readVCFContigs(obj.getFileName())
            except (OSError, RuntimeError, ValueError, EOFError) as error:
                errors.append(str(error))
            if obj.getAssembly() and obj.getAssembly() != self.assembly.get():
                errors.append('VCF and VEP assemblies differ.')
            if obj.getScientificName() and obj.getScientificName().lower().replace(' ', '_') != self.species.get():
                errors.append('VCF and VEP species differ.')

    def annotateStep(self):
        os.makedirs(self._getExtraPath(), exist_ok=True)
        errors = self._validate()
        if errors:
            raise RuntimeError('; '.join(errors))
        with open(self._getCommandFile(), 'w') as commandFile:
            commandFile.write('')
        with tempfile.TemporaryDirectory(prefix='vep_', dir=self._getExtraPath()) as work:
            raw = os.path.join(work, 'all.tsv')
            self.runTool('vep', ['--input_file', self.inputVCF.get().getFileName(),
                '--output_file', raw, '--format', 'vcf', '--tab', '--cache', '--offline',
                '--dir_cache', os.path.abspath(self.cacheDirectory.get()),
                '--cache_version', self.cacheVersion.get(), '--species', self.species.get(),
                '--assembly', self.assembly.get(), '--uniprot', '--protein', '--symbol',
                '--transcript_version', '--tsl', '--fork', self.numberOfThreads.get(),
                '--no_stats', '--force_overwrite'])
            annotationCount = sum(1 for _ in iterVepRows(raw))
            self.info('Validated {} VEP annotations.'.format(annotationCount))
            selected = os.path.join(work, 'missense.tsv')
            if self.selectMissense.get():
                self.runTool('filter_vep', ['-i', raw, '-o', selected,
                    '--filter', 'Consequence is missense_variant', '--force_overwrite'])
                for row in iterVepRows(selected):
                    if 'missense_variant' not in row['Consequence'].replace('&', ',').split(','):
                        raise RuntimeError('filter_vep produced a non-missense effect.')
            destination = self._getExtraPath('annotation')
            if os.path.isdir(destination):
                shutil.rmtree(destination)
            shutil.move(work, destination)

    def createOutputStep(self):
        obj = self.inputVCF.get()
        raw = VariantAnnotation(filename=self._getExtraPath('annotation', 'all.tsv'),
            scientificName=obj.getScientificName() or self.species.get().replace('_', ' ').capitalize(),
            assembly=self.assembly.get(), cacheVersion=self.cacheVersion.get(),
            command=self.readCommandHistory())
        selected = raw.clone()
        if self.selectMissense.get():
            selected.setFileName(self._getExtraPath('annotation', 'missense.tsv'))
        self._defineOutputs(outputAllAnnotations=raw, outputAnnotations=selected)
        self._defineSourceRelation(obj, raw)
        self._defineSourceRelation(obj, selected)

    def _summary(self):
        return ['Offline VEP; missense selection: {}.'.format(self.selectMissense.get())]
