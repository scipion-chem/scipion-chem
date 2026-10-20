
# **************************************************************************
# *
# * Authors:    Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import shlex

from pyworkflow.protocol.params import (
    BooleanParam, IntParam, PointerParam, StringParam
)
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import VCFFile


class ProtGATKMutect2(EMProtocol):
    """Call candidate somatic variants with GATK Mutect2.

    Supports tumor-only and tumor-normal calling. Optionally applies
    FilterMutectCalls. Input alignments must refer to the same reference
    assembly and contain valid sample names in their read groups.

    The output VCF contains candidate somatic variants, not necessarily
    confirmed somatic mutations. Additional filtering and validation may
    be required, especially for RNA-seq and tumor-only analyses.
    """

    _label = 'GATK Mutect2'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam(
            'inputTumor', PointerParam, pointerClass='AlignmentFile',
            label='Tumor alignment:',
            help='Tumor BAM (AlignmentFile) previously prepared for variant calling.'
        )
        form.addParam(
            'useNormal', BooleanParam, default=False,
            label='Use matched normal:',
            help='Enable tumor-normal variant calling.'
        )
        form.addParam(
            'inputNormal', PointerParam, pointerClass='AlignmentFile',
            condition='useNormal', label='Normal alignment:',
            help='Matched normal BAM with the same reference assembly.'
        )
        form.addParam(
            'normalSample', StringParam, default='',
            condition='useNormal', label='Normal sample name:',
            help='Sample name (SM) from the normal BAM read groups. '
                 'Must match the BAM header exactly.'
        )

        form.addSection(label='Optional resources')
        form.addParam(
            'useGermlineResource', BooleanParam, default=False,
            label='Use germline resource:'
        )
        form.addParam(
            'germlineVCFs', PointerParam, pointerClass='SetOfVCFFiles',
            condition='useGermlineResource', label='Germline VCFs:'
        )
        form.addParam(
            'germlineIndex', IntParam, default=0,
            condition='useGermlineResource', label='Germline VCF index:',
            help='Zero-based index of the resource in the VCF set.'
        )
        form.addParam(
            'usePanelOfNormals', BooleanParam, default=False,
            label='Use panel of normals:'
        )
        form.addParam(
            'ponVCFs', PointerParam, pointerClass='SetOfVCFFiles',
            condition='usePanelOfNormals', label='Panel of normals VCFs:'
        )
        form.addParam(
            'ponIndex', IntParam, default=0,
            condition='usePanelOfNormals', label='Panel of normals VCF index:',
            help='Zero-based index of the PoN in the VCF set.'
        )

        form.addSection(label='Variant filtering')
        form.addParam(
            'filterCalls', BooleanParam, default=True,
            label='Run FilterMutectCalls:',
            help='Apply GATK FilterMutectCalls to Mutect2 calls. '
                 'Variants failing filters remain in the VCF with FILTER tags.'
        )
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareInputsStep)
        self._insertFunctionStep(self.mutect2Step)
        if self.filterCalls.get():
            self._insertFunctionStep(self.filterMutectCallsStep)
        self._insertFunctionStep(self.createOutputStep)

    def prepareInputsStep(self):
        tumor = self.inputTumor.get()
        reference = self._getReferenceFasta(tumor)
        self._ensureFastaIndex(reference)
        self._ensureSequenceDictionary(reference)

        if self.useGermlineResource.get():
            self._ensureFeatureIndex(
                self._getSelectedVcfPath(self.germlineVCFs, self.germlineIndex)
            )
        if self.usePanelOfNormals.get():
            self._ensureFeatureIndex(
                self._getSelectedVcfPath(self.ponVCFs, self.ponIndex)
            )

    def mutect2Step(self):
        tumor = self.inputTumor.get()
        reference = self._getReferenceFasta(tumor)
        args = [
            'Mutect2', '-R', reference,
            '-I', tumor.getFileName(),
            '-O', self._getUnfilteredVcf(),
            '--native-pair-hmm-threads', str(self.numberOfThreads.get())
        ]
        if self.useNormal.get():
            args += ['-I', self.inputNormal.get().getFileName(),
                     '-normal', self.normalSample.get().strip()]
        if self.useGermlineResource.get():
            args += ['--germline-resource', self._getSelectedVcfPath(
                self.germlineVCFs, self.germlineIndex)]
        if self.usePanelOfNormals.get():
            args += ['--panel-of-normals', self._getSelectedVcfPath(
                self.ponVCFs, self.ponIndex)]

        self._runGatk(args)
        self._checkNonempty(self._getUnfilteredVcf(), 'Mutect2 VCF')
        if self.filterCalls.get():
            self._checkNonempty(self._getStatsFile(), 'Mutect2 statistics')

    def filterMutectCallsStep(self):
        args = [
            'FilterMutectCalls', '-R',
            self._getReferenceFasta(self.inputTumor.get()),
            '-V', self._getUnfilteredVcf(),
            '-O', self._getFilteredVcf()
        ]
        self._runGatk(args)
        self._checkNonempty(self._getFilteredVcf(), 'Filtered Mutect2 VCF')

    def createOutputStep(self):
        tumor = self.inputTumor.get()
        reference = self._getReferenceFasta(tumor)
        unfiltered = self._makeVcf(self._getUnfilteredVcf(), reference)
        outputs = {'outputUnfilteredVCF': unfiltered}
        if self.filterCalls.get():
            filtered = self._makeVcf(self._getFilteredVcf(), reference)
            outputs['outputVCF'] = filtered
        else:
            outputs['outputVCF'] = unfiltered.clone()

        self._defineOutputs(**outputs)
        for output in outputs.values():
            self._defineSourceRelation(self.inputTumor, output)
            if self.useNormal.get():
                self._defineSourceRelation(self.inputNormal, output)
            if self.useGermlineResource.get():
                self._defineSourceRelation(self.germlineVCFs, output)
            if self.usePanelOfNormals.get():
                self._defineSourceRelation(self.ponVCFs, output)

    @staticmethod
    def _makeVcf(path, reference):
        vcf = VCFFile()
        vcf.setFileName(path)
        index = next((path + ext for ext in ('.idx', '.tbi', '.csi')
                      if os.path.isfile(path + ext)), None)
        if index and hasattr(vcf, 'setIndexFile'):
            vcf.setIndexFile(index)
        if hasattr(vcf, 'setReferenceFasta'):
            vcf.setReferenceFasta(reference)
        return vcf

    def _getUnfilteredVcf(self):
        return self._getExtraPath('mutect2_unfiltered.vcf.gz')

    def _getFilteredVcf(self):
        return self._getExtraPath('mutect2_filtered.vcf.gz')

    def _getStatsFile(self):
        return self._getUnfilteredVcf() + '.stats'

    @staticmethod
    def _getReferenceFasta(alignment):
        if alignment is not None and hasattr(alignment, 'getReferenceFasta'):
            return alignment.getReferenceFasta()
        return None

    @staticmethod
    def _getSelectedVcfPath(pointer, indexParam):
        vcfSet = pointer.get()
        if vcfSet is None or vcfSet.getSize() == 0:
            raise RuntimeError('A non-empty SetOfVCFFiles is required.')
        index = indexParam.get()
        if index < 0 or index >= vcfSet.getSize():
            raise RuntimeError('VCF index {} is out of range (size {}).'.format(
                index, vcfSet.getSize()))
        for i, vcf in enumerate(vcfSet):
            if i == index:
                path = vcf.getFileName()
                if not path or not os.path.isfile(path):
                    raise RuntimeError('Selected VCF does not exist: {}'.format(path))
                return path
        raise RuntimeError('Unable to retrieve the selected VCF.')

    def _ensureFastaIndex(self, reference):
        if os.path.isfile(reference + '.fai'):
            return
        Plugin.runCondaCommand(
            self, 'faidx {}'.format(shlex.quote(reference)),
            RNASEQ_DIC, 'samtools')
        self._checkNonempty(reference + '.fai', 'FASTA index')

    @staticmethod
    def _getDictionaryPath(reference):
        root, extension = os.path.splitext(reference)
        if extension.lower() in ('.fa', '.fasta', '.fna'):
            return root + '.dict'
        return reference + '.dict'

    def _ensureSequenceDictionary(self, reference):
        dictionary = self._getDictionaryPath(reference)
        if os.path.isfile(dictionary):
            return
        self._runGatk(['CreateSequenceDictionary', '-R', reference,
                       '-O', dictionary])
        self._checkNonempty(dictionary, 'Sequence dictionary')

    def _ensureFeatureIndex(self, vcf):
        # GATK accepts Tribble indexes for plain VCF and tabix indexes for
        # block-gzip-compressed VCF. A CSI index is not used here as the
        # only available index, matching the existing ProtGATK convention.
        candidates = [vcf + '.idx', vcf + '.tbi']
        if any(os.path.isfile(path) for path in candidates):
            return
        self._runGatk(['IndexFeatureFile', '-I', vcf])
        if not any(os.path.isfile(path) for path in candidates):
            raise RuntimeError('No GATK-compatible index for VCF: {}'.format(vcf))

    def _runGatk(self, arguments):
        command = ' '.join(shlex.quote(str(arg)) for arg in arguments)
        Plugin.runCondaCommand(self, command, RNASEQ_DIC, 'gatk')
        with open(self._getExtraPath('gatk_mutect2_commands.txt'), 'a') as handle:
            handle.write('gatk ' + command + '\n')

    @staticmethod
    def _checkNonempty(path, label):
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise RuntimeError('{} is missing or empty: {}'.format(label, path))

    def _validate(self):
        errors = []
        tumor = self.inputTumor.get()
        if tumor is None:
            return ['A tumor AlignmentFile is required.']
        try:
            self._checkNonempty(tumor.getFileName(), 'Tumor BAM')
        except (RuntimeError, TypeError):
            errors.append('Tumor BAM is missing or empty.')
        reference = self._getReferenceFasta(tumor)
        if not reference or not os.path.isfile(reference):
            errors.append('Tumor alignment must contain an existing reference FASTA.')

        if self.useNormal.get():
            normal = self.inputNormal.get()
            if normal is None:
                errors.append('A normal AlignmentFile is required.')
            else:
                try:
                    self._checkNonempty(normal.getFileName(), 'Normal BAM')
                except (RuntimeError, TypeError):
                    errors.append('Normal BAM is missing or empty.')
                normalReference = self._getReferenceFasta(normal)
                if not normalReference or not os.path.isfile(normalReference):
                    errors.append('Normal alignment must contain an existing reference FASTA.')
                elif reference and os.path.realpath(normalReference) != os.path.realpath(reference):
                    errors.append('Tumor and normal must use the same reference FASTA.')
                if normal.getFileName() == tumor.getFileName():
                    errors.append('Tumor and normal BAM paths must differ.')
            if not self.normalSample.get().strip():
                errors.append('Enter the normal sample name from the BAM SM tag.')

        for enabled, pointer, index, label in (
            (self.useGermlineResource, self.germlineVCFs,
             self.germlineIndex, 'Germline resource'),
            (self.usePanelOfNormals, self.ponVCFs,
             self.ponIndex, 'Panel of normals')
        ):
            if enabled.get():
                try:
                    self._getSelectedVcfPath(pointer, index)
                except (RuntimeError, AttributeError) as error:
                    errors.append('{}: {}'.format(label, error))
        return errors

    def _summary(self):
        summary = [
            'Calling mode: {}'.format(
                'Tumor-normal' if self.useNormal.get() else 'Tumor-only'),
            'FilterMutectCalls: {}'.format(
                'enabled' if self.filterCalls.get() else 'disabled')
        ]
        if hasattr(self, 'outputVCF'):
            summary.append('Output VCF: {}'.format(self.outputVCF.getFileName()))
        return summary

    def _methods(self):
        if not hasattr(self, 'outputVCF'):
            return []
        return [
            'Candidate somatic variants were called with GATK Mutect2 '
            'using {} analysis{}. The result was stored in a VCF file.'
            .format('tumor-normal' if self.useNormal.get() else 'tumor-only',
                    ' and FilterMutectCalls' if self.filterCalls.get() else '')
        ]