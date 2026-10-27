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
from pyworkflow.protocol.params import EnumParam, FileParam, PointerParam, StringParam

from pwchem.objects import RNASeqExpression
from pwchem.utils.utilsRNA import GATKPreparationMixin, checkGtfAnnotation


class ProtCufflinks(GATKPreparationMixin, EMProtocol):
    """Assemble transcripts with Cufflinks or compare expression with Cuffdiff.

    This protocol implements the expression branch of the reference RNA-seq
    script. It runs independently of Mutect2, using prepared BAM alignments
    and a GTF annotation. Variant calling results are not expression inputs.

    Modes
    -----
    Cufflinks (one sample)
        Runs the equivalent of ``cufflinks BAM -g GTF --num-threads 4 -o DIR``.
        The lowercase ``-g`` option preserves the script's guided assembly;
        it is not replaced by reference-only quantification with ``-G``.
        The assembled transcripts and gene/transcript FPKM tables are retained.

    Cuffdiff (two samples)
        Runs the equivalent of ``cuffdiff GTF BAM1 BAM2 -L tumor,normal -o DIR``.
        The first and second BAM are distinct inputs. Labels default to
        ``tumor`` and ``normal`` and can be changed for other comparisons,
        including cell-line/allogeneic samples. Labels do not establish a
        biological relationship between samples. This initial implementation
        accepts one BAM per condition, matching the supplied script; it does
        not pool replicates or run cuffmerge/cuffquant.

    Input parameters
    ----------------
    inputAlignment : AlignmentFile
        First coordinate-sorted BAM, normally produced by STAR, Picard and
        GATK preprocessing. The reference script uses its final reordered
        BAM. Preprocessing is not repeated here. The input must have a SAM
        sequence dictionary and a header declaring coordinate sorting.

    inputComparison : AlignmentFile, optional
        Second BAM for Cuffdiff. It must use the same sequence dictionary,
        including sequence lengths and order. Structured species/assembly
        metadata must agree when supplied. If both BAMs identify reference
        FASTA paths, they must resolve to the same file; aliases are accepted.
        These checks cannot certify biological pairing or reference MD5s.

    annotationFile : file, optional
        Uncompressed GTF2 annotation. When empty, the GTF path is taken from
        the first alignment's reference metadata. An explicit annotation
        overrides that path and is used for both conditions. A matching
        assembly and compatible transcript identifiers remain essential for
        later variant-to-protein processing.

    firstLabel, secondLabel : str
        Cuffdiff condition labels in BAM input order. Labels must be distinct,
        non-empty identifiers containing letters, digits, dots, underscores
        or hyphens. Defaults reproduce ``-L tumor,normal``. They are stored
        separately from the biological sample names in the output metadata.

    libraryType : enum
        ``fr-unstranded`` by default, with ``fr-firststrand`` and
        ``fr-secondstrand`` available for stranded libraries. This is the
        library protocol, not STAR's intron-motif option or the orientation
        of an individual read. No strandedness is inferred from the BAM.

    Execution and validation
    ------------------------
    Input BAMs and the annotation are linked into the protocol directory.
    samtools quickcheck and extracted headers check basic BAM integrity and
    coordinate order. A BAM index is not required by this stage and is not
    generated merely for expression analysis. The GTF is read in streaming
    mode to check nine-column records, coordinates and transcript features.
    At least one transcript feature must overlap a BAM sequence; features
    on extra annotation contigs are reported and left for Cufflinks to handle.
    Cufflinks retains its own checks for alignment tags and transcript models.

    Commands use the selected thread count and disable the tool's update
    check. Other statistical and assembly settings retain Cufflinks defaults,
    as in the reference script. Analysis runs in a fresh temporary directory;
    required output headers are checked before replacing the previous local
    result directory. Missing outputs cannot be hidden by an earlier run.
    Input files and their indexes are not modified.

    Output
    ------
    outputExpression : RNASeqExpression
        Primary file: ``isoforms.fpkm_tracking``. The object also records
        ``genes.fpkm_tracking``, annotation/reference metadata, sample names,
        condition labels and executed commands. Cufflinks adds the assembled
        ``transcripts.gtf``; Cuffdiff adds ``isoform_exp.diff`` and
        ``gene_exp.diff``. Other files generated by the tool remain in the
        same output directory. Tracking tables with a valid header and no
        data rows are accepted; expression or significance is not fabricated.

    Both modes retain the native table format. A one-sample Cufflinks table
    has one FPKM column, whereas Cuffdiff has condition-specific FPKM columns.
    A downstream peptide protocol must select the appropriate condition and
    handle these schemas explicitly. This protocol does not convert either
    table to amino-acid sequences or filter differential results to significant
    genes. Statuses such as NOTEST/LOWDATA remain in the native output.

    Reproducing the reference workflow
    ---------------------------------
    Run each sample through FASTQ import/FastQC, fastp with minimum length 40,
    STAR, Picard read groups/duplicate marking, GATK SplitNCigarReads/BQSR,
    and Picard ReorderSam. Connect the final BAM here and to Mutect2 as two
    branches. For a single tumor or separate cell/allogeneic samples without
    a normal, use separate Cufflinks runs. For a two-sample expression
    comparison, choose Cuffdiff and provide both BAMs explicitly.

    The script's allogeneic branch labels a comparison as tumor/normal even
    when its second BAM is the allogeneic sample. This protocol exposes the
    comparison explicitly rather than inferring that sample's role. A shared
    normal can be supplied to separate comparisons when that is the intended
    experiment. Differential analysis here reproduces the one-BAM-per-condition
    script and does not add biological replication or certify its conclusions.

    Requirements
    ------------
    Cufflinks/Cuffdiff 2.2.1 and samtools must be available in the RNA-seq
    environment configured by the plugin. The annotation must be decompressed
    before import. The complete RNA-seq workflow still requires separate VCF
    intersection, VEP annotation and peptide-generation stages.
    """

    _label = 'Cufflinks / Cuffdiff'
    CUFFLINKS = 0
    CUFFDIFF = 1
    CUFFDIFF_CONDITION = 'mode == {}'.format(CUFFDIFF)
    MODES = ['Cufflinks (one sample)', 'Cuffdiff (two samples)']
    LIBRARY_TYPES = ['fr-unstranded', 'fr-firststrand', 'fr-secondstrand']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('mode', EnumParam, default=self.CUFFLINKS, choices=self.MODES,
                      label='Analysis mode:')
        form.addParam('inputAlignment', PointerParam, pointerClass='AlignmentFile',
                      label='First alignment:', help='Prepared coordinate-sorted BAM.')
        form.addParam('inputComparison', PointerParam, pointerClass='AlignmentFile',
                      condition=self.CUFFDIFF_CONDITION, label='Comparison alignment:',
                      help='Second BAM for Cuffdiff; one BAM per condition.')
        form.addParam('annotationFile', FileParam, default='', allowsNull=True,
                      label='Annotation GTF:',
                      help='Optional override. Otherwise use the first alignment reference GTF.')
        form.addSection(label='Expression options')
        form.addParam('firstLabel', StringParam, default='tumor', condition=self.CUFFDIFF_CONDITION,
                      label='First condition label:')
        form.addParam('secondLabel', StringParam, default='normal', condition=self.CUFFDIFF_CONDITION,
                      label='Second condition label:')
        form.addParam('libraryType', EnumParam, default=0, choices=self.LIBRARY_TYPES,
                      label='Library type:')
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareInputsStep)
        self._insertFunctionStep(self.expressionStep)
        self._insertFunctionStep(self.createOutputStep)

    def _getAlignments(self):
        alignments = [self.inputAlignment.get()]
        if self.mode.get() == self.CUFFDIFF:
            alignments.append(self.inputComparison.get())
        return alignments

    def _getAnnotation(self):
        override = self.annotationFile.get()
        alignment = self.inputAlignment.get()
        path = override or (alignment.getReferenceGtf() if alignment is not None else None)
        return os.path.abspath(path) if path else None

    def _getLabels(self):
        return [(self.firstLabel.get() or '').strip(), (self.secondLabel.get() or '').strip()]

    def _getProgram(self):
        return 'cuffdiff' if self.mode.get() == self.CUFFDIFF else 'cufflinks'

    def _getOutputDir(self):
        return self._getExtraPath(self._getProgram())

    def _getCommandFile(self):
        return self._getExtraPath('cufflinks_commands.txt')

    def _getAlignmentPath(self, index):
        return self._getExtraPath('inputs', 'sample{}.bam'.format(index + 1))

    def _getAnnotationPath(self):
        return self._getExtraPath('inputs', 'annotation.gtf')

    def prepareInputsStep(self):
        errors = self._validate()
        if errors:
            raise RuntimeError('Invalid expression inputs: ' + '; '.join(errors))
        os.makedirs(self._getExtraPath(), exist_ok=True)
        with open(self._getCommandFile(), 'w') as commandFile:
            commandFile.write('')
        firstContigs = None
        for index, alignment in enumerate(self._getAlignments()):
            bam = self._getAlignmentPath(index)
            self.stageFile(alignment.getFileName(), bam)
            self.runTool('samtools', ['quickcheck', '-v', bam])
            header = bam + '.header.sam'
            self.runTool('samtools', ['view', '-H', '-o', header, bam])
            contigs, _, order = self.readSamHeader(header)
            if order != 'coordinate':
                raise RuntimeError('Sample {} BAM must declare coordinate sorting.'.format(index + 1))
            if firstContigs is None:
                firstContigs = contigs
            elif list(firstContigs.items()) != list(contigs.items()):
                raise RuntimeError('Cuffdiff BAMs must use the same sequence names, lengths and order.')
        annotation = self._getAnnotationPath()
        self.stageFile(self._getAnnotation(), annotation)
        self._checkAnnotation(annotation, firstContigs)

    def _checkAnnotation(self, path, contigs):
        extraContigs = checkGtfAnnotation(path, contigs)
        if extraContigs:
            self.warning('GTF contains {} sequence(s) absent from the BAM dictionary.'.format(len(extraContigs)))

    def expressionStep(self):
        program = self._getProgram()
        os.makedirs(self._getExtraPath(), exist_ok=True)
        with tempfile.TemporaryDirectory(prefix=program + '_', dir=self._getExtraPath()) as runDir:
            arguments = ['-o', runDir, '--num-threads', self.numberOfThreads.get(),
                         '--library-type', self.LIBRARY_TYPES[self.libraryType.get()], '--no-update-check']
            if self.mode.get() == self.CUFFDIFF:
                arguments.extend(['-L', ','.join(self._getLabels()), self._getAnnotationPath(),
                                  self._getAlignmentPath(0), self._getAlignmentPath(1)])
            else:
                arguments.extend(['-g', self._getAnnotationPath(), self._getAlignmentPath(0)])
            self.runTool(program, arguments)
            self._checkOutputs(runDir)
            target = self._getOutputDir()
            if os.path.isdir(target):
                shutil.rmtree(target)
            os.replace(runDir, target)

    @staticmethod
    def _checkTable(path, required):
        GATKPreparationMixin.checkNonempty(path, 'Expression table')
        with open(path) as handle:
            columns = handle.readline().rstrip('\r\n').split('\t')
        if len(columns) != len(set(columns)) or not set(required).issubset(columns):
            raise RuntimeError('Unexpected expression table header: {}'.format(path))
        return columns

    def _checkOutputs(self, directory):
        for filename in ('isoforms.fpkm_tracking', 'genes.fpkm_tracking'):
            columns = self._checkTable(os.path.join(directory, filename), ['tracking_id', 'gene_id', 'locus'])
            if self.mode.get() == self.CUFFLINKS:
                valid = 'FPKM' in columns
            else:
                valid = len([column for column in columns if column.endswith('_FPKM')]) == 2
            if not valid:
                raise RuntimeError('Missing expression columns in {}.'.format(filename))
        if self.mode.get() == self.CUFFDIFF:
            required = ['test_id', 'sample_1', 'sample_2', 'status', 'p_value', 'q_value', 'significant']
            for filename in ('isoform_exp.diff', 'gene_exp.diff'):
                self._checkTable(os.path.join(directory, filename), required)
        elif not os.path.isfile(os.path.join(directory, 'transcripts.gtf')):
            raise RuntimeError('Cufflinks did not produce transcripts.gtf.')

    def createOutputStep(self):
        directory = self._getOutputDir()
        self._checkOutputs(directory)
        first = self.inputAlignment.get()
        output = RNASeqExpression(filename=os.path.join(directory, 'isoforms.fpkm_tracking'))
        output.setGenesFile(os.path.join(directory, 'genes.fpkm_tracking'))
        output.setProgram(self._getProgram())
        output.setSampleName(first.getSampleName())
        output.setReferenceGtf(self._getAnnotationPath())
        output.setReferenceFasta(first.getReferenceFasta())
        output.setScientificName(first.getReferenceScientificName())
        output.setAssembly(first.getReferenceAssembly())
        output.setCommand(self.readCommandHistory())
        if self.mode.get() == self.CUFFDIFF:
            output.setComparisonSample(self.inputComparison.get().getSampleName())
            firstLabel, secondLabel = self._getLabels()
            output.setSampleLabel(firstLabel)
            output.setComparisonLabel(secondLabel)
            output.setIsoformDiffFile(os.path.join(directory, 'isoform_exp.diff'))
            output.setGeneDiffFile(os.path.join(directory, 'gene_exp.diff'))
        else:
            output.setSampleLabel(first.getSampleName())
            output.setTranscriptsFile(os.path.join(directory, 'transcripts.gtf'))
        self._defineOutputs(outputExpression=output)
        self._defineSourceRelation(self.inputAlignment, output)
        if self.mode.get() == self.CUFFDIFF:
            self._defineSourceRelation(self.inputComparison, output)

    def _validate(self):
        errors = []
        if self.mode.get() not in (self.CUFFLINKS, self.CUFFDIFF):
            errors.append('Select a valid expression mode.')
        alignments = self._getAlignments()
        errors.extend(self._validateAlignments(alignments))
        errors.extend(self._validateAnnotation())
        if self.mode.get() == self.CUFFDIFF:
            errors.extend(self._validateCuffdiff(alignments))
        if self.numberOfThreads.get() < 1:
            errors.append('At least one thread is required.')
        if self.libraryType.get() not in range(len(self.LIBRARY_TYPES)):
            errors.append('Select a valid library type.')
        return errors

    def _validateAlignments(self, alignments):
        errors = []
        for index, alignment in enumerate(alignments, start=1):
            if alignment is None:
                errors.append('Sample {} AlignmentFile is required.'.format(index))
                continue
            if alignment.getFormat() != 'BAM':
                errors.append('Sample {} alignment must be BAM.'.format(index))
            try:
                self.checkNonempty(alignment.getFileName(), 'Sample {} BAM'.format(index))
                if alignment.getReferenceFasta():
                    self.checkNonempty(alignment.getReferenceFasta(), 'Sample {} reference'.format(index))
            except RuntimeError as error:
                errors.append(str(error))
        return errors

    def _validateAnnotation(self):
        errors = []
        try:
            annotation = self._getAnnotation()
            self.checkNonempty(annotation, 'Annotation GTF')
            with open(annotation, 'rb') as handle:
                if annotation.lower().endswith('.gz') or handle.read(2) == b'\x1f\x8b':
                    errors.append('Annotation GTF must be uncompressed.')
        except (RuntimeError, OSError) as error:
            errors.append(str(error))
        return errors

    def _validateCuffdiff(self, alignments):
        errors = []
        labels = self._getLabels()
        if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', label) for label in labels) or labels[0] == labels[1]:
            errors.append('Cuffdiff labels must be distinct, non-empty identifiers.')
        if ',' in self._getExtraPath('inputs'):
            errors.append('Cuffdiff working paths cannot contain commas (reserved for replicate lists).')
        if all(alignment is not None for alignment in alignments):
            errors.extend(self._validateComparisonReferences(alignments))
        return errors

    @staticmethod
    def _validateComparisonReferences(alignments):
        errors = []
        paths = [alignment.getFileName() for alignment in alignments]
        if all(path and os.path.isfile(path) for path in paths) and os.path.samefile(*paths):
            errors.append('Cuffdiff requires different BAM files, including aliases.')
        refs = [alignment.getReferenceFasta() for alignment in alignments]
        if all(path and os.path.isfile(path) for path in refs) and not os.path.samefile(*refs):
            errors.append('Cuffdiff BAMs must identify the same reference FASTA file.')
        for getter in ('getReferenceScientificName', 'getReferenceAssembly'):
            values = [getattr(alignment, getter)() for alignment in alignments]
            if all(values) and values[0] != values[1]:
                errors.append('Cuffdiff reference metadata disagree: {}.'.format(getter))
        return errors

    def _summary(self):
        summary = ['Expression analysis: {}'.format(self._getProgram())]
        if hasattr(self, 'outputExpression'):
            summary.append('Transcript expression: {}'.format(self.outputExpression.getFileName()))
            if self.mode.get() == self.CUFFDIFF:
                summary.append('Conditions: {} -> {}'.format(*self._getLabels()))
        return summary

    def _methods(self):
        if not hasattr(self, 'outputExpression'):
            return []
        return ['Transcript and gene expression were analysed with {} using the supplied GTF '
                'and {} coordinate-sorted BAM alignment(s). Native FPKM tables were retained.'
                .format(self._getProgram(), len(self._getAlignments()))]
