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

    This protocol accepts an ``AlignmentFile`` and applies selected Picard
    operations sequentially to the associated BAM file. It is intended for
    BAM preparation and processing within the RNA-seq workflow.

    The following Picard operations are supported:

    AddOrReplaceReadGroups
        Adds or replaces read-group metadata in the BAM header. Read groups
        identify the biological sample, sequencing library and sequencing
        platform and are required by several downstream tools, particularly
        GATK.

        The protocol uses the following read-group values::

            RGID = id
            RGLB = library
            RGPL = ILLUMINA
            RGPU = machine

        ``RGSM`` is obtained automatically from the sample name stored in
        the input ``AlignmentFile``. If no sample name is available,
        ``sample`` is used.

    MarkDuplicates
        Identifies reads that are likely to represent PCR or optical
        duplicates.

        Duplicate reads are marked but are not removed from the BAM file.
        A Picard duplication metrics file is also generated.

    ReorderSam
        Reorders the BAM according to the sequence dictionary of the
        reference genome associated with the input ``AlignmentFile``.

        This operation requires a valid reference FASTA stored in the
        alignment metadata.

        In the RNA-seq workflow, this operation can be used after the GATK
        BAM-processing steps.


    ---------------------------------------------------------------------
    Execution order
    ---------------------------------------------------------------------

    When several operations are selected, they are always executed in the
    following order::

        AddOrReplaceReadGroups
            -> MarkDuplicates
            -> ReorderSam

    Operations that are not selected are skipped while preserving this
    execution order.

    At least one Picard operation must be selected.


    ---------------------------------------------------------------------
    Input
    ---------------------------------------------------------------------

    inputAlignment : AlignmentFile
        Input BAM alignment.

        The BAM path is obtained from the ``AlignmentFile`` object.

        The sample name stored in the alignment is used by
        ``AddOrReplaceReadGroups``.

        The reference FASTA stored in the alignment is used by
        ``ReorderSam``.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputAlignment : AlignmentFile
        Final processed BAM alignment.

        The output object is created by cloning the input
        ``AlignmentFile`` and updating the BAM and BAM index paths. This
        preserves the metadata associated with the original alignment.

        The final BAM depends on the last selected Picard operation:

        - ``ReorderSam`` output, when enabled.
        - Otherwise, ``MarkDuplicates`` output, when enabled.
        - Otherwise, ``AddOrReplaceReadGroups`` output.

        The final BAM is indexed with samtools and the resulting ``.bai``
        file is associated with the output ``AlignmentFile``.

        When supported by the ``AlignmentFile`` object, the output is
        marked as indexed and sorted and the executed Picard commands are
        stored as command metadata.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - An input ``AlignmentFile`` has been provided.
    - The input BAM exists and is not empty.
    - At least one Picard operation has been selected.
    - A valid reference FASTA is available when ``ReorderSam`` is enabled.


    ---------------------------------------------------------------------
    Intermediate files
    ---------------------------------------------------------------------

    Depending on the selected operations, the protocol can generate:

    add_read_groups.bam
        BAM generated by ``AddOrReplaceReadGroups``.

    mark_duplicates.bam
        BAM generated by ``MarkDuplicates``.

    mark_duplicates_metrics.txt
        Duplication metrics generated by ``MarkDuplicates``.

    reordered.bam
        BAM generated by ``ReorderSam``.

    picard_commands.txt
        Internal record of the Picard commands executed by the protocol.


    ---------------------------------------------------------------------
    Requirements
    ---------------------------------------------------------------------

    The following programs must be available through the RNA-seq
    environment configured by the Scipion-Chem plugin:

    - Picard
    - samtools


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Picard operations are executed in a fixed order.

    - ``MarkDuplicates`` marks duplicate reads but does not remove them.

    - ``ReorderSam`` requires the reference FASTA associated with the
      input ``AlignmentFile``.

    - Intermediate BAM files remain in the protocol working directory.

    - The final BAM is always indexed with samtools.

    - The output ``AlignmentFile`` preserves the metadata of the input
      alignment.
    """

    _label = 'Picard processing'

    PICARD_COMMAND = 'picard {}'
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
        self._appendCommand(self.PICARD_COMMAND.format(args))

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

        self._appendCommand(self.PICARD_COMMAND.format(args))

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

        # ReorderSam expects a sequence dictionary, not the FASTA itself.
        referenceDict = os.path.splitext(referenceFasta)[0] + '.dict'
        if not os.path.isfile(referenceDict):
            referenceDict = self._getExtraPath('reference.dict')
            dictArgs = 'CreateSequenceDictionary R={} O={}'.format(
                self._quote(referenceFasta), self._quote(referenceDict)
            )
            Plugin.runCondaCommand(self, dictArgs, RNASEQ_DIC, 'picard')
            self._appendCommand(self.PICARD_COMMAND.format(dictArgs))

        if not os.path.isfile(referenceDict) or os.path.getsize(referenceDict) == 0:
            raise RuntimeError('Reference sequence dictionary was not created: {}'.format(referenceDict))

        args = (
            'ReorderSam '
            'I={inputBam} '
            'O={outputBam} '
            'SD={referenceDict} '
            'CREATE_INDEX=false'
        ).format(
            inputBam=self._quote(inputBam),
            outputBam=self._quote(outputBam),
            referenceDict=self._quote(referenceDict)
        )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'picard'
        )

        self._validateBamFile(outputBam)
        self._appendCommand(self.PICARD_COMMAND.format(args))

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