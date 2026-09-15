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
import shlex

from pyworkflow.protocol.params import BooleanParam, PointerParam
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC



class ProtPicard(EMProtocol):
    """
    Process BAM alignments using one or more Picard operations.

    The protocol accepts an AlignmentFile object and allows several Picard
    operations to be applied sequentially to the same BAM file.

    Available operations
    --------------------

    AddOrReplaceReadGroups
        Adds or replaces read-group metadata in the BAM header. Read groups
        identify the biological sample, sequencing library and sequencing
        platform. This information is required by several downstream tools,
        particularly GATK.

        The RNA-seq workflow uses fixed values for the technical read-group
        fields:

            RGID = id
            RGLB = library
            RGPL = ILLUMINA
            RGPU = machine

        RGSM is obtained automatically from the sample name stored in the
        input AlignmentFile.

    MarkDuplicates
        Identifies reads that are likely to be duplicates generated during
        library preparation or sequencing. Duplicates are marked in the BAM
        file but are retained.

        A Picard metrics file describing duplicate statistics is also created.

    ReorderSam
        Reorders the BAM sequence dictionary and records according to the
        reference genome stored in the input AlignmentFile.

        This operation requires the reference FASTA associated with the
        alignment.

    Execution order
    ---------------

    When more than one operation is selected, operations are always executed
    in the following order:

        AddOrReplaceReadGroups
            -> MarkDuplicates
            -> ReorderSam

    The user selects which operations are enabled, but not their order.

    The final BAM is indexed with samtools and returned as an AlignmentFile.
    Intermediate BAM files remain inside the protocol working directory.
    """

    _label = 'Picard processing'

    # -------------------------------------------------------
    # Parameters
    # -------------------------------------------------------

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'inputAlignment',
            PointerParam,
            pointerClass='AlignmentFile',
            label='Input alignment: ',
            help=(
                'Input BAM represented as an AlignmentFile object. '
                'The output of the RNA-seq alignment protocol can be used '
                'directly.'
            )
        )

        form.addSection(label='Picard operations')

        form.addParam(
            'addReadGroups',
            BooleanParam,
            default=True,
            label='Add or replace read groups: ',
            help=(
                'Run Picard AddOrReplaceReadGroups. This adds read-group '
                'metadata to the BAM header, including sample, library and '
                'sequencing-platform information. Read groups are required '
                'by several downstream GATK tools. The sample name is taken '
                'automatically from the input AlignmentFile.'
            )
        )

        form.addParam(
            'markDuplicates',
            BooleanParam,
            default=True,
            label='Mark duplicates: ',
            help=(
                'Run Picard MarkDuplicates. Reads that are likely to be PCR '
                'or optical duplicates are marked in the BAM file. They are '
                'not removed. A duplication metrics file is also generated.'
            )
        )

        form.addParam(
            'reorderBam',
            BooleanParam,
            default=False,
            label='Reorder BAM: ',
            help=(
                'Run Picard ReorderSam. This makes the BAM sequence order '
                'match the reference genome associated with the AlignmentFile. '
                'In the original RNA-seq workflow this step is normally used '
                'after the GATK BAM-processing steps.'
            )
        )

    # -------------------------------------------------------
    # Steps
    # -------------------------------------------------------

    def _insertAllSteps(self):
        if self.addReadGroups.get():
            self._insertFunctionStep(self.addReadGroupsStep)

        if self.markDuplicates.get():
            self._insertFunctionStep(self.markDuplicatesStep)

        if self.reorderBam.get():
            self._insertFunctionStep(self.reorderBamStep)

        self._insertFunctionStep(self.indexFinalBamStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------------------------------------
    # Picard operations
    # -------------------------------------------------------

    def addReadGroupsStep(self):
        inputAlignment = self.inputAlignment.get()
        inputBam = inputAlignment.getFileName()
        outputBam = self._getAddReadGroupsBam()

        self._validateBamFile(inputBam)

        sampleName = self._getInputSampleName(inputAlignment) or 'sample'

        args = (
            'AddOrReplaceReadGroups '
            'I={inputBam} '
            'O={outputBam} '
            'SO=coordinate '
            'RGID=id '
            'RGLB=library '
            'RGPL=ILLUMINA '
            'RGPU=machine '
            'RGSM={sampleName}'
        ).format(
            inputBam=self._quote(inputBam),
            outputBam=self._quote(outputBam),
            sampleName=self._quote(sampleName)
        )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'picard'
        )

        self._validateBamFile(outputBam)
        self._appendCommand('picard {}'.format(args))

    def markDuplicatesStep(self):
        inputBam = self._getInputForMarkDuplicates()
        outputBam = self._getMarkDuplicatesBam()
        metricsFile = self._getMarkDuplicatesMetricsFile()

        self._validateBamFile(inputBam)

        args = (
            'MarkDuplicates '
            'I={inputBam} '
            'O={outputBam} '
            'M={metricsFile} '
            'REMOVE_DUPLICATES=false '
            'CREATE_INDEX=false '
            'VALIDATION_STRINGENCY=SILENT'
        ).format(
            inputBam=self._quote(inputBam),
            outputBam=self._quote(outputBam),
            metricsFile=self._quote(metricsFile)
        )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'picard'
        )

        self._validateBamFile(outputBam)

        if not os.path.isfile(metricsFile):
            raise RuntimeError(
                'Picard duplication metrics file was not created: {}'
                .format(metricsFile)
            )

        self._appendCommand('picard {}'.format(args))

    def reorderBamStep(self):
        inputAlignment = self.inputAlignment.get()
        inputBam = self._getInputForReorderSam()
        outputBam = self._getReorderedBam()
        referenceFasta = self._getReferenceFasta(inputAlignment)

        self._validateBamFile(inputBam)

        if not referenceFasta:
            raise RuntimeError(
                'ReorderSam requires a reference FASTA stored in the '
                'input AlignmentFile.'
            )

        if not os.path.isfile(referenceFasta):
            raise RuntimeError(
                'Reference FASTA does not exist: {}'.format(referenceFasta)
            )

        args = (
            'ReorderSam '
            'I={inputBam} '
            'O={outputBam} '
            'R={referenceFasta} '
            'CREATE_INDEX=false'
        ).format(
            inputBam=self._quote(inputBam),
            outputBam=self._quote(outputBam),
            referenceFasta=self._quote(referenceFasta)
        )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'picard'
        )

        self._validateBamFile(outputBam)
        self._appendCommand('picard {}'.format(args))

    # -------------------------------------------------------
    # Final BAM
    # -------------------------------------------------------

    def indexFinalBamStep(self):
        finalBam = self._getFinalBam()

        self._validateBamFile(finalBam)

        Plugin.runCondaCommand(
            self,
            'index {}'.format(self._quote(finalBam)),
            RNASEQ_DIC,
            'samtools'
        )

        baiFile = finalBam + '.bai'

        if not os.path.isfile(baiFile):
            raise RuntimeError(
                'BAM index was not created: {}'.format(baiFile)
            )

    def createOutputStep(self):
        inputAlignment = self.inputAlignment.get()

        finalBam = self._getFinalBam()
        finalBai = finalBam + '.bai'

        outputAlignment = inputAlignment.clone()
        outputAlignment.setFileName(finalBam)
        outputAlignment.setIndexFile(finalBai)
        outputAlignment.setIsIndexed(True)

        if hasattr(outputAlignment, 'setIsSorted'):
            outputAlignment.setIsSorted(True)

        if hasattr(outputAlignment, 'setCommand'):
            outputAlignment.setCommand(self._readCommandHistory())

        self._defineOutputs(
            outputAlignment=outputAlignment
        )

        self._defineSourceRelation(
            self.inputAlignment,
            outputAlignment
        )

    # -------------------------------------------------------
    # BAM flow
    # -------------------------------------------------------

    def _getInputForMarkDuplicates(self):
        if self.addReadGroups.get():
            return self._getAddReadGroupsBam()

        return self.inputAlignment.get().getFileName()

    def _getInputForReorderSam(self):
        if self.markDuplicates.get():
            return self._getMarkDuplicatesBam()

        if self.addReadGroups.get():
            return self._getAddReadGroupsBam()

        return self.inputAlignment.get().getFileName()

    def _getFinalBam(self):
        if self.reorderBam.get():
            return self._getReorderedBam()

        if self.markDuplicates.get():
            return self._getMarkDuplicatesBam()

        if self.addReadGroups.get():
            return self._getAddReadGroupsBam()

        return self.inputAlignment.get().getFileName()

    # -------------------------------------------------------
    # Paths
    # -------------------------------------------------------

    def _getAddReadGroupsBam(self):
        return self._getExtraPath('add_read_groups.bam')

    def _getMarkDuplicatesBam(self):
        return self._getExtraPath('mark_duplicates.bam')

    def _getReorderedBam(self):
        return self._getExtraPath('reordered.bam')

    def _getMarkDuplicatesMetricsFile(self):
        return self._getExtraPath('mark_duplicates_metrics.txt')

    def _getCommandFile(self):
        return self._getExtraPath('picard_commands.txt')

    # -------------------------------------------------------
    # Command history
    # -------------------------------------------------------

    def _appendCommand(self, command):
        with open(self._getCommandFile(), 'a') as outputFile:
            outputFile.write(command + '\n')

    def _readCommandHistory(self):
        commandFile = self._getCommandFile()

        if not os.path.isfile(commandFile):
            return ''

        with open(commandFile) as inputFile:
            commands = [
                line.strip()
                for line in inputFile
                if line.strip()
            ]

        return ' && '.join(commands)

    # -------------------------------------------------------
    # Alignment metadata
    # -------------------------------------------------------

    @staticmethod
    def _getInputSampleName(alignment):
        if hasattr(alignment, 'getSampleName'):
            return alignment.getSampleName()

        return ''

    @staticmethod
    def _getReferenceFasta(alignment):
        if hasattr(alignment, 'getReferenceFasta'):
            return alignment.getReferenceFasta()

        return None

    # -------------------------------------------------------
    # Validation helpers
    # -------------------------------------------------------

    @staticmethod
    def _validateBamFile(bamFile):
        if not bamFile:
            raise RuntimeError('The BAM path is empty.')

        if not os.path.isfile(bamFile):
            raise RuntimeError(
                'BAM file does not exist: {}'.format(bamFile)
            )

        if os.path.getsize(bamFile) == 0:
            raise RuntimeError(
                'BAM file is empty: {}'.format(bamFile)
            )

    @staticmethod
    def _quote(value):
        return shlex.quote(str(value))

    def _validate(self):
        errors = []

        alignment = self.inputAlignment.get()

        if alignment is None:
            errors.append(
                'An input AlignmentFile object is required.'
            )
            return errors

        inputBam = alignment.getFileName()

        if not inputBam:
            errors.append(
                'The input AlignmentFile does not contain a BAM path.'
            )
        elif not os.path.isfile(inputBam):
            errors.append(
                'The input BAM file does not exist: {}'
                .format(inputBam)
            )
        elif os.path.getsize(inputBam) == 0:
            errors.append(
                'The input BAM file is empty: {}'
                .format(inputBam)
            )

        if not (
            self.addReadGroups.get()
            or self.markDuplicates.get()
            or self.reorderBam.get()
        ):
            errors.append(
                'Select at least one Picard operation.'
            )

        if self.reorderBam.get():
            referenceFasta = self._getReferenceFasta(alignment)

            if not referenceFasta:
                errors.append(
                    'ReorderSam requires a reference FASTA stored in '
                    'the input AlignmentFile.'
                )
            elif not os.path.isfile(referenceFasta):
                errors.append(
                    'The reference FASTA does not exist: {}'
                    .format(referenceFasta)
                )

        return errors

    # -------------------------------------------------------
    # Summary and methods
    # -------------------------------------------------------

    def _getSelectedOperationNames(self):
        operations = []

        if self.addReadGroups.get():
            operations.append('AddOrReplaceReadGroups')

        if self.markDuplicates.get():
            operations.append('MarkDuplicates')

        if self.reorderBam.get():
            operations.append('ReorderSam')

        return operations

    def _summary(self):
        if not hasattr(self, 'outputAlignment'):
            return [
                'Selected Picard operations: {}'.format(
                    ' -> '.join(self._getSelectedOperationNames())
                )
            ]

        summary = [
            'Picard operations: {}'.format(
                ' -> '.join(self._getSelectedOperationNames())
            ),
            'Output BAM: {}'.format(
                self.outputAlignment.getFileName()
            )
        ]

        if self.markDuplicates.get():
            summary.append(
                'Duplication metrics: {}'.format(
                    self._getMarkDuplicatesMetricsFile()
                )
            )

        return summary

    def _methods(self):
        if not hasattr(self, 'outputAlignment'):
            return []

        operations = ' followed by '.join(
            self._getSelectedOperationNames()
        )

        return [
            'The input BAM alignment was processed with Picard using {}. '
            'The resulting BAM file was indexed with samtools.'
            .format(operations)
        ]
