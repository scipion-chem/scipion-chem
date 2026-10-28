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

from pyworkflow.protocol.params import BooleanParam, PointerParam
from pwem.protocols import EMProtocol

from pwchem.utils.utilsRNA import GATKPreparationMixin


class ProtPicard(GATKPreparationMixin, EMProtocol):
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
        ``sample`` is used. Leading and trailing whitespace is removed.
        The resulting sample name must contain printable ASCII characters,
        as required for SAM read-group values. Existing read groups are
        replaced by a single group; library and platform values are the
        fixed values shown above, rather than inferred from sequencing data.

    MarkDuplicates
        Identifies reads that are likely to represent PCR or optical
        duplicates.

        Duplicate reads are marked but are not removed from the BAM file.
        A non-empty Picard duplication metrics file is also generated.
        When read-group replacement is disabled, the input header must already
        declare read groups with unique, non-empty IDs. This does not certify
        that every alignment record carries the correct RG tag.

    ReorderSam
        Reorders the BAM according to the sequence dictionary of the
        reference genome associated with the input ``AlignmentFile``.

        This operation requires a valid reference FASTA stored in the
        alignment metadata. Reference auxiliary files are prepared locally:
        usable FASTA indexes and dictionaries are copied; missing, outdated or
        invalid ones are generated with samtools and Picard. Sequence names,
        lengths and order must agree between the local index and dictionary.
        The input BAM may use a different sequence order, but its sequences
        must exist in the reference with matching lengths.

        In the RNA-seq workflow, this operation can be used after the GATK
        BAM-processing steps.


    ---------------------------------------------------------------------
    Execution order
    ---------------------------------------------------------------------

    When several operations are selected, they are always executed in the
    following order::

        Input validation and optional reference preparation
            -> AddOrReplaceReadGroups
            -> MarkDuplicates
            -> ReorderSam
            -> Final BAM indexing and output creation

    Operations that are not selected are skipped while preserving this
    execution order.

    At least one Picard operation must be selected. AddOrReplaceReadGroups
    requests coordinate sorting. When that operation is disabled, the input
    BAM header must already declare coordinate order, so the final result
    can be indexed. Sorting is not inferred from the input object metadata.


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
        file is associated with the output ``AlignmentFile``. Empty or
        invalid indexes stop execution rather than being published.

        When supported by the ``AlignmentFile`` object, the output is
        marked as indexed and sorted and the executed Picard commands are
        stored as command metadata. The sample name is also updated when
        read groups are replaced, including use of the ``sample`` fallback.
        Reference metadata are preserved, including the scientific name,
        assembly and original reference FASTA path. Local preparation does
        not modify the source reference or its auxiliary files.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - An input ``AlignmentFile`` has been provided.
    - The alignment format is BAM, and the input file exists and is not empty.
    - At least one Picard operation has been selected.
    - A non-empty reference FASTA is available when ``ReorderSam`` is enabled.
    - At least one thread is available for final BAM indexing.
    - An explicitly supplied sample name can be written into SAM read groups.

    During preparation, samtools quickcheck checks the input BAM header and
    end-of-file marker, and the extracted SAM header is inspected for sequence
    declarations, sort order and duplicate read-group IDs. When reordering is
    enabled, the reference must be uncompressed and its sequence lengths must
    match the BAM. The order may differ, since changing it is ReorderSam's job.

    These checks reuse utilsRNA routines. They do not scan every alignment,
    verify reference MD5s, establish biological sample identity or replace
    Picard validation. MarkDuplicates retains VALIDATION_STRINGENCY=SILENT
    to match the reference workflow; that setting is not full BAM validation.


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
        Record of the Picard and samtools commands executed by the protocol.

    input.header.sam
        Input BAM header used to check sequence and read-group declarations.

    reference/reference.fa, reference/reference.fa.fai, reference/reference.dict
        Local reference link, FASTA index and sequence dictionary used for
        ReorderSam. Existing source auxiliary files are copied when usable.
        Index checks inspect signatures, timestamps and dictionary structure;
        they do not certify the integrity of every index entry.


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
      Previous BAMs, indexes and metrics at the generated output paths are
      removed before each corresponding operation, preventing stale results
      from hiding a missing output during repeated execution.

    - The final BAM is always indexed with samtools.

    - The output ``AlignmentFile`` preserves the metadata of the input
      alignment.
    """

    CREATE_INDEX_OPTION = 'CREATE_INDEX=false'

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
        form.addParallelSection(threads=4, mpi=1)

    # -------------------------------------------------------
    # Steps
    # -------------------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareInputsStep)
        if self.addReadGroups.get():
            self._insertFunctionStep(self.addReadGroupsStep)
        if self.markDuplicates.get():
            self._insertFunctionStep(self.markDuplicatesStep)
        if self.reorderBam.get():
            self._insertFunctionStep(self.reorderBamStep)
        self._insertFunctionStep(self.indexFinalBamStep)
        self._insertFunctionStep(self.createOutputStep)

    def prepareInputsStep(self):
        errors = self._validate()
        if errors:
            raise RuntimeError('Invalid Picard inputs: ' + '; '.join(errors))
        os.makedirs(self._getExtraPath(), exist_ok=True)
        with open(self._getCommandFile(), 'w') as commandFile:
            commandFile.write('')
        bam = self.inputAlignment.get().getFileName()
        self.runTool('samtools', ['quickcheck', '-v', bam])
        header = self._getExtraPath('input.header.sam')
        self.runTool('samtools', ['view', '-H', '-o', header, bam])
        contigs, groups, sortOrder = self.readSamHeader(header)
        if not self.addReadGroups.get():
            self._checkExistingReadGroups(groups, sortOrder)
        if self.reorderBam.get():
            reference = self.prepareReference(self.getReferenceFasta(self.inputAlignment.get()))
            for name, length in contigs.items():
                if reference.get(name) != length:
                    raise RuntimeError('Input BAM sequence {} (length {}) is absent from the reference '
                                       'or has a different length.'.format(name, length))

    def _checkExistingReadGroups(self, groups, sortOrder):
        if sortOrder != 'coordinate':
            raise RuntimeError('Input BAM must declare coordinate sorting when '
                               'AddOrReplaceReadGroups is disabled.')
        ids = [group.get('ID', '').strip() for group in groups]
        if any(not groupId for groupId in ids) or len(ids) != len(set(ids)):
            raise RuntimeError('Input BAM contains empty or duplicate read-group IDs.')
        if self.markDuplicates.get() and not groups:
            raise RuntimeError('MarkDuplicates requires read groups. Enable AddOrReplaceReadGroups '
                               'or supply a BAM with existing read groups.')

    def ensureSequenceDictionary(self, referenceFasta):
        """Use Picard for reference preparation without requiring GATK."""
        dictionary = self.getSequenceDictionaryPath(referenceFasta)
        if self.isCurrentIndex(dictionary, referenceFasta) and self.validSequenceDictionary(dictionary):
            return
        self.removeIndex(dictionary)
        self.runTool('picard', ['CreateSequenceDictionary', 'R=' + referenceFasta, 'O=' + dictionary])
        if not self.validSequenceDictionary(dictionary):
            raise RuntimeError('Picard did not create a valid sequence dictionary: {}'.format(dictionary))

    # -------------------------------------------------------
    # Picard operations
    # -------------------------------------------------------

    def _runPicard(self, arguments, outputBam):
        self.clearBamOutput(outputBam)
        self.runTool('picard', arguments)
        self.validateBamFile(outputBam)

    def addReadGroupsStep(self):
        alignment = self.inputAlignment.get()
        self._runPicard([
            'AddOrReplaceReadGroups', 'I=' + alignment.getFileName(),
            'O=' + self._getAddReadGroupsBam(), 'SO=coordinate', self.CREATE_INDEX_OPTION,
            'RGID=id', 'RGLB=library', 'RGPL=ILLUMINA', 'RGPU=machine',
            'RGSM=' + self._getInputSampleName(alignment)
        ], self._getAddReadGroupsBam())

    def markDuplicatesStep(self):
        metrics = self._getMarkDuplicatesMetricsFile()
        self.removeIndex(metrics)
        self._runPicard([
            'MarkDuplicates', 'I=' + self._getInputForMarkDuplicates(),
            'O=' + self._getMarkDuplicatesBam(), 'M=' + metrics,
            'REMOVE_DUPLICATES=false', self.CREATE_INDEX_OPTION, 'VALIDATION_STRINGENCY=SILENT'
        ], self._getMarkDuplicatesBam())
        self.checkNonempty(metrics, 'Picard duplication metrics')

    def reorderBamStep(self):
        dictionary = self.getSequenceDictionaryPath(self.getPreparedReference())
        self.checkNonempty(dictionary, 'Reference sequence dictionary')
        self._runPicard([
            'ReorderSam', 'I=' + self._getInputForReorderSam(),
            'O=' + self._getReorderedBam(), 'SD=' + dictionary, self.CREATE_INDEX_OPTION
        ], self._getReorderedBam())

    # -------------------------------------------------------
    # Output
    # -------------------------------------------------------

    def indexFinalBamStep(self):
        self.ensureBamIndex(self._getFinalBam())

    def createOutputStep(self):
        finalBam = self._getFinalBam()
        finalBai = finalBam + '.bai'
        self.validateBamFile(finalBam)
        if not self.isCurrentIndex(finalBai, finalBam) or not self.hasIndexSignature(finalBai, b'BAI\x01'):
            raise RuntimeError('Final BAM index is missing, outdated or invalid: {}'.format(finalBai))
        alignment = self.inputAlignment.get()
        output = alignment.clone()
        output.setFileName(finalBam)
        output.setIndexFile(finalBai)
        output.setIsIndexed(True)
        output.setIsSorted(True)
        output.setCommand(self.readCommandHistory())
        if self.addReadGroups.get():
            output.setSampleName(self._getInputSampleName(alignment))
        self._defineOutputs(outputAlignment=output)
        self._defineSourceRelation(self.inputAlignment, output)

    # -------------------------------------------------------
    # BAM flow and paths
    # -------------------------------------------------------

    def _getInputForMarkDuplicates(self):
        return (self._getAddReadGroupsBam() if self.addReadGroups.get()
                else self.inputAlignment.get().getFileName())

    def _getInputForReorderSam(self):
        return self._getMarkDuplicatesBam() if self.markDuplicates.get() else self._getInputForMarkDuplicates()

    def _getFinalBam(self):
        if self.reorderBam.get():
            return self._getReorderedBam()
        if self.markDuplicates.get():
            return self._getMarkDuplicatesBam()
        if self.addReadGroups.get():
            return self._getAddReadGroupsBam()
        raise RuntimeError('Select at least one Picard operation.')

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

    @staticmethod
    def _getInputSampleName(alignment):
        return (alignment.getSampleName() or '').strip() or 'sample'

    # -------------------------------------------------------
    # Validation and reporting
    # -------------------------------------------------------

    def _validate(self):
        errors = []
        alignment = self.inputAlignment.get()
        if alignment is None:
            return ['An input AlignmentFile object is required.']
        try:
            self.validateBamFile(alignment.getFileName())
        except RuntimeError as error:
            errors.append(str(error))
        if alignment.getFormat() != 'BAM':
            errors.append('Input alignment must be BAM.')
        if not self._getSelectedOperationNames():
            errors.append('Select at least one Picard operation.')
        if self.reorderBam.get():
            try:
                self.checkNonempty(self.getReferenceFasta(alignment), 'Reference FASTA')
            except RuntimeError as error:
                errors.append(str(error))
        if self.addReadGroups.get():
            sample = self._getInputSampleName(alignment)
            if any(not 32 <= ord(character) <= 126 for character in sample):
                errors.append('Sample name must contain only printable ASCII characters for SAM read groups.')
        if self.numberOfThreads.get() < 1:
            errors.append('At least one thread is required.')
        return errors

    def _getSelectedOperationNames(self):
        return [name for selected, name in (
            (self.addReadGroups.get(), 'AddOrReplaceReadGroups'),
            (self.markDuplicates.get(), 'MarkDuplicates'),
            (self.reorderBam.get(), 'ReorderSam')
        ) if selected]

    def _summary(self):
        operations = ' -> '.join(self._getSelectedOperationNames())
        if not hasattr(self, 'outputAlignment'):
            return ['Selected Picard operations: {}'.format(operations)]
        summary = ['Picard operations: {}'.format(operations),
                   'Output BAM: {}'.format(self.outputAlignment.getFileName())]
        if self.markDuplicates.get():
            summary.append('Duplication metrics: {}'.format(self._getMarkDuplicatesMetricsFile()))
        return summary

    def _methods(self):
        if not hasattr(self, 'outputAlignment'):
            return []
        return ['The input BAM alignment was processed with Picard using {}. '
                'The resulting BAM file was indexed with samtools.'
                .format(' followed by '.join(self._getSelectedOperationNames()))]
