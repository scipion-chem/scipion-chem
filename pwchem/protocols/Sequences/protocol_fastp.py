# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
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

import json
import os
import re
import shlex

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import BooleanParam, PointerParam, IntParam

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import FastqFile
from pwchem.utils.utilsRNA import getFastqStats, runFastqc


class ProtFastpFilter(EMProtocol):
    """
      Filter and preprocess FASTQ reads using fastp.

      This protocol performs quality filtering and optional adapter trimming
      on single-end or paired-end FASTQ datasets using fastp. It accepts a
      ``FastqFile`` object as input and generates a new ``FastqFile`` containing
      the filtered reads for downstream RNA-seq analyses.

      Filtering can be configured according to read length, base quality,
      percentage of unqualified bases, number of ambiguous bases (N), and
      average read quality. Adapter trimming can also be enabled or disabled.

      After filtering, the protocol calculates the number of reads and the
      mean read length from the generated FASTQ file(s). For paired-end data,
      both reads are processed together and the protocol verifies that the
      resulting R1 and R2 files contain the same number of reads.

      fastp HTML and JSON reports are generated automatically. Optionally,
      FastQC can be executed on the filtered reads to provide an additional
      quality-control assessment.


      ---------------------------------------------------------------------
      Input parameters
      ---------------------------------------------------------------------

      inputFastq : FastqFile
          Input FASTQ dataset.

          Both single-end and paired-end datasets are supported. Paired-end
          inputs must contain both R1 and R2 FASTQ files.

      runFastqc : bool
          If True, FastQC is executed on the filtered FASTQ file(s).

          FastQC is performed after fastp filtering.

      lengthRequired : int
          Minimum read length required after filtering. The default of 40
          matches the reference RNA-seq script.

          Reads shorter than this value are discarded.

          Corresponds to the fastp option::

              -l / --length_required

      qualifiedQualityPhred : int
          Minimum Phred quality score for a base to be considered qualified.

          Corresponds to the fastp option::

              -q / --qualified_quality_phred

      unqualifiedPercentLimit : int
          Maximum percentage of unqualified bases allowed in a read.

          Reads exceeding this percentage are discarded.

          Corresponds to the fastp option::

              -u / --unqualified_percent_limit

      nBaseLimit : int
          Maximum number of ambiguous N bases allowed in a read.

          Reads containing more N bases than this value are discarded.

          Corresponds to the fastp option::

              -n / --n_base_limit

      averageQual : int
          Minimum average quality score required for a read.

          Reads with an average quality below this value are discarded.

          Corresponds to the fastp option::

              -e / --average_qual

      disableAdapterTrimming : bool
          If True, automatic adapter trimming performed by fastp is disabled.

          Corresponds to the fastp option::

              -A / --disable_adapter_trimming


      ---------------------------------------------------------------------
      Output
      ---------------------------------------------------------------------

      outputFastq : FastqFile
          Filtered FASTQ dataset.

          The output object contains:

          - Filtered FASTQ file path(s).
          - Sample name.
          - Sequencing type (single-end or paired-end).
          - Number of reads after filtering.
          - Mean read length after filtering.
          - fastp HTML report.
          - fastp JSON report.
          - FastQC HTML report(s), when requested and reads remain.

          The output has a source relation to the input FastqFile. For paired-end
          data, the stored read count is the number of pairs (records per file),
          while the fastp JSON total_reads count includes both mates.

          Filtered FASTQ files are stored uncompressed.


      ---------------------------------------------------------------------
      Workflow
      ---------------------------------------------------------------------

      The protocol performs the following steps:

      1. Retrieve the input FASTQ file(s) and sequencing metadata.

      2. Execute fastp using the selected filtering parameters.

      3. Generate the filtered FASTQ file(s) together with the fastp HTML
         and JSON reports.

      4. Calculate the number of reads and mean read length from the filtered
         FASTQ file(s).

      5. For paired-end datasets, verify that R1 and R2 contain the same
         number of reads.

      6. Check that the FASTQ record count agrees with the fastp JSON report.
         If requested and reads remain, execute FastQC on the filtered reads.

      7. Create the output ``FastqFile`` containing the filtered files,
         statistics and quality-control reports.


      ---------------------------------------------------------------------
      Notes
      ---------------------------------------------------------------------

      - Both single-end and paired-end FASTQ datasets are supported.

      - Input files must contain FASTQ quality scores.

      - The minimum read length is not limited by the stored input mean read
        length: variable-length datasets can contain longer reads. A value
        above the mean may therefore be valid, but can discard many reads.

      - Filtered FASTQ files are written in uncompressed FASTQ format.

      - Output filenames are generated from a sanitized version of the input
        sample name. Directory separators and special characters are replaced
        with underscores; the original name is preserved in output metadata.

      - If the input sample does not have a sample name, ``sample`` is used.

      - For paired-end datasets, R1 and R2 are filtered together.

      - The protocol verifies that filtered paired-end files contain the
        same number of reads.

      - For paired-end datasets, the stored read length corresponds to the
        mean of the R1 and R2 mean read lengths.

      - FastQC, when enabled, is executed on the filtered reads rather than
        on the original input files.

      - fastp HTML and JSON reports are required and must be non-empty. The
        JSON report must be readable and its retained-read count must match
        the FASTQ outputs. This does not validate every FASTQ record field.

      - If no reads pass filtering, the empty FASTQ outputs and fastp reports
        are retained, the output count is zero and optional FastQC is skipped
        with a warning. Review filtering thresholds before downstream alignment.

      - Input FASTQ files must exist, be non-empty and, for paired-end data,
        refer to different files even when symbolic links are used.

      - Parameters follow fastp 0.24.0 limits: qualified and average quality
        0-93, unqualified bases 0-100 percent, N bases 0-50, minimum length
        zero or greater, and 1-16 threads.

      - fastp HTML and JSON reports are always generated. Previous files at
        the generated output paths are removed before execution, so an old
        report or FASTQ cannot hide a missing result on a repeated run.

      - Generated fastp and FastQC reports are stored as attributes of the
        output ``FastqFile`` object.


      ---------------------------------------------------------------------
      Requirements
      ---------------------------------------------------------------------

      The following programs must be available through the RNA-seq
      environment configured by the Scipion-Chem plugin:

      - fastp
      - FastQC, when ``runFastqc`` is enabled
      """

    _label = 'fastp filter'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputFastq', PointerParam,
                      pointerClass='FastqFile',
                      label='Input FASTQ: ',
                      help='Input FASTQ dataset to be filtered.')

        form.addParam('runFastqc', BooleanParam,
                      default=False,
                      label='Run FastQC: ',
                      help='Execute FastQC on the filtered FASTQ files.')

        form.addSection(label='Fastp parameters')

        form.addParam('lengthRequired', IntParam,
                      default=40,
                      label='Minimum read length: ',
                      help='Reads shorter than this value will be discarded. '
                           'Corresponds to fastp option -l / --length_required.')

        form.addParam('qualifiedQualityPhred', IntParam,
                      default=15,
                      label='Qualified quality phred: ',
                      help='Minimum base quality considered qualified. '
                           'Corresponds to fastp option -q / --qualified_quality_phred.')

        form.addParam('unqualifiedPercentLimit', IntParam,
                      default=40,
                      label='Unqualified percent limit: ',
                      help='Maximum percentage of unqualified bases allowed per read. '
                           'Corresponds to fastp option -u / --unqualified_percent_limit.')

        form.addParam('nBaseLimit', IntParam,
                      default=5,
                      label='N base limit: ',
                      help='Reads with more N bases than this value will be discarded. '
                           'Corresponds to fastp option -n / --n_base_limit.')

        form.addParam('averageQual', IntParam,
                      default=0,
                      label='Average quality: ',
                      help='Reads with average quality below this value will be discarded. '
                           'Corresponds to fastp option -e / --average_qual.')

        form.addParam('disableAdapterTrimming', BooleanParam,
                      default=False,
                      label='Disable adapter trimming: ',
                      help='Disable fastp adapter trimming. '
                           'Corresponds to fastp option -A / --disable_adapter_trimming.')

        form.addParallelSection(threads=4, mpi=0)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.filterStep)

    def filterStep(self):
        errors = self._validate()
        if errors:
            raise RuntimeError('Invalid fastp inputs: ' + '; '.join(errors))
        inputFastq = self.inputFastq.get()

        outputFastq = self._runFastp(inputFastq)

        if self.runFastqc.get() and outputFastq.getNumReads() > 0:
            htmlFiles = self._runFastqc(outputFastq)

            if outputFastq.isPaired():
                outputFastq.setFastqcHtmlR1(htmlFiles[0])
                outputFastq.setFastqcHtmlR2(htmlFiles[1])
            else:
                outputFastq.setFastqcHtml(htmlFiles[0])

        if self.runFastqc.get() and outputFastq.getNumReads() == 0:
            self.warning('No reads passed fastp filtering. FastQC was skipped; review the fastp reports.')

        self._defineOutputs(outputFastq=outputFastq)
        self._defineSourceRelation(self.inputFastq, outputFastq)

    def _runFastp(self, inputFastq):
        isPaired = inputFastq.isPaired()
        sampleName = inputFastq.getSampleName() or 'sample'
        fileStem = re.sub(r'[^A-Za-z0-9_.-]+', '_', sampleName).strip('.') or 'sample'
        outDir = self._getExtraPath('fastp')
        os.makedirs(outDir, exist_ok=True)

        out1 = os.path.join(outDir, f'{fileStem}_filtered_R1.fastq')
        out2 = os.path.join(outDir, f'{fileStem}_filtered_R2.fastq') if isPaired else None
        reportHtml = os.path.join(outDir, f'{fileStem}_fastp.html')
        reportJson = os.path.join(outDir, f'{fileStem}_fastp.json')

        args = [
            '-i', inputFastq.getFileName(), '-o', out1,
            '-w', self.numberOfThreads.get(), '-l', self.lengthRequired.get(),
            '-q', self.qualifiedQualityPhred.get(), '-u', self.unqualifiedPercentLimit.get(),
            '-n', self.nBaseLimit.get(), '-e', self.averageQual.get(),
            '-h', reportHtml, '-j', reportJson
        ]
        if isPaired:
            args.extend(['-I', inputFastq.getFileName2(), '-O', out2])
        if self.disableAdapterTrimming.get():
            args.append('-A')

        inputFiles = [inputFastq.getFileName()]
        if isPaired:
            inputFiles.append(inputFastq.getFileName2())
        for path in (out1, out2, reportHtml, reportJson):
            if path and any(os.path.abspath(path) == os.path.abspath(source)
                            or (os.path.exists(path) and os.path.samefile(path, source)) for source in inputFiles):
                raise RuntimeError('fastp output would overwrite an input file: {}'.format(path))
        for path in (out1, out2, reportHtml, reportJson):
            if path and os.path.lexists(path):
                os.remove(path)
        arguments = ' '.join(shlex.quote(str(value)) for value in args)
        Plugin.runCondaCommand(self, arguments, RNASEQ_DIC, 'fastp')

        outputFiles = [out1, out2] if isPaired else [out1]
        for path in outputFiles:
            if not os.path.isfile(path):
                raise RuntimeError(f'fastp did not generate the expected FASTQ: {path}')
        for path in (reportHtml, reportJson):
            if not os.path.isfile(path) or os.path.getsize(path) == 0:
                raise RuntimeError(f'fastp report is missing or empty: {path}')

        numReads, readLength = getFastqStats(out1)
        if isPaired:
            numReads2, readLength2 = getFastqStats(out2)
            if numReads != numReads2:
                raise RuntimeError('Filtered paired FASTQ files have different numbers of reads: '
                                   f'R1={numReads} R2={numReads2}')
            readLength = int(round((readLength + readLength2) / 2))

        try:
            with open(reportJson) as handle:
                report = json.load(handle)
            reportedReads = report['summary']['after_filtering']['total_reads']
            if type(reportedReads) is not int or reportedReads < 0:
                raise ValueError('Invalid total_reads value')
        except (OSError, ValueError, KeyError, TypeError) as error:
            raise RuntimeError(f'Invalid fastp JSON report: {reportJson}') from error
        expectedReads = numReads * (2 if isPaired else 1)
        if reportedReads != expectedReads:
            raise RuntimeError('Filtered FASTQ read count does not match the fastp report: '
                               f'files={expectedReads}, report={reportedReads}')

        outputFastq = FastqFile()
        outputFastq.setFileName(out1)
        outputFastq.setIsPaired(isPaired)
        outputFastq.setIsCompressed(False)
        outputFastq.setFormat('FASTQ')
        outputFastq.setHasQuality(True)
        outputFastq.setReadLength(readLength)
        outputFastq.setNumReads(numReads)
        outputFastq.setSampleName(sampleName)
        if isPaired:
            outputFastq.setFileName2(out2)
        outputFastq.setFastpHtml(reportHtml)
        outputFastq.setFastpJson(reportJson)
        return outputFastq

    def _runFastqc(self, outputFastq):
        if not outputFastq.supportsFastQC():
            raise RuntimeError(
                'FastQC cannot be executed because the input file does not '
                'contain quality scores.'
            )

        fastqFiles = [outputFastq.getFileName()]

        if outputFastq.isPaired():
            fastqFiles.append(outputFastq.getFileName2())

        return runFastqc(self, fastqFiles)

    def _validate(self):
        errors = []
        inputFastq = self.inputFastq.get()
        if inputFastq is None:
            return ['An input FASTQ object is required.']

        files = [('Read 1', inputFastq.getFileName())]
        if inputFastq.isPaired():
            files.append(('Read 2', inputFastq.getFileName2()))
        for label, path in files:
            if not path or not os.path.isfile(path):
                errors.append(f'{label} FASTQ file is missing: {path}')
            elif os.path.getsize(path) == 0:
                errors.append(f'{label} FASTQ file is empty: {path}')
        if (len(files) == 2 and all(path and os.path.isfile(path) for _, path in files)
                and os.path.samefile(files[0][1], files[1][1])):
            errors.append('Paired-end reads must use different FASTQ files.')
        if inputFastq.getFormat() != 'FASTQ' or not inputFastq.hasQuality():
            errors.append('Input must be FASTQ with quality scores.')

        if self.lengthRequired.get() < 0:
            errors.append('Minimum read length must be greater than or equal to 0.')
        for value, label, maximum in (
            (self.qualifiedQualityPhred.get(), 'Qualified quality phred', 93),
            (self.unqualifiedPercentLimit.get(), 'Unqualified percent limit', 100),
            (self.nBaseLimit.get(), 'N base limit', 50),
            (self.averageQual.get(), 'Average quality', 93)
        ):
            if not 0 <= value <= maximum:
                errors.append(f'{label} must be between 0 and {maximum}.')
        if not 1 <= self.numberOfThreads.get() <= 16:
            errors.append('fastp supports between 1 and 16 threads.')
        return errors

    def _summary(self):
        summary = []

        inputFastq = self.inputFastq.get()

        summary.append('FASTQ filtering with fastp')
        summary.append(f'Minimum read length: {self.lengthRequired.get()}')
        summary.append(f'Threads: {self.numberOfThreads.get()}')
        summary.append(f'Qualified quality phred: {self.qualifiedQualityPhred.get()}')
        summary.append(f'Unqualified percent limit: {self.unqualifiedPercentLimit.get()}')
        summary.append(f'N base limit: {self.nBaseLimit.get()}')
        summary.append(f'Average quality: {self.averageQual.get()}')
        summary.append(
            f'Adapter trimming: '
            f'{"disabled" if self.disableAdapterTrimming.get() else "enabled"}'
        )

        if inputFastq:
            summary.append(f'Input: {inputFastq.getFileName()}')
            if (inputFastq.getReadLength() or 0) > 0:
                summary.append(
                    f'Input read length: {inputFastq.getReadLength()} bp'
                )
            if (inputFastq.getNumReads() or 0) > 0:
                summary.append(
                    f'Input number of reads: {inputFastq.getNumReads()}'
                )
            summary.append(
                f'Sequencing type: '
                f'{"paired-end" if inputFastq.isPaired() else "single-end"}'
            )

        if hasattr(self, 'outputFastq'):
            outputFastq = self.outputFastq

            if outputFastq.getReadLength() > 0:
                summary.append(
                    f'Output read length: {outputFastq.getReadLength()} bp'
                )

            summary.append(f'Output number of reads: {outputFastq.getNumReads()}')

        if not self.runFastqc.get():
            fastqcStatus = 'not requested'
        elif not hasattr(self, 'outputFastq'):
            fastqcStatus = 'pending'
        elif self.outputFastq.getNumReads() == 0:
            fastqcStatus = 'skipped (no reads passed filtering)'
        else:
            fastqcStatus = 'executed'
        summary.append(f'FastQC: {fastqcStatus}')

        return summary

    def _methods(self):
        if not hasattr(self, 'outputFastq'):
            return []
        methods = []

        methods.append(
            'Reads were filtered using fastp with the following parameters: '
            f'minimum read length {self.lengthRequired.get()} bp, '
            f'{self.numberOfThreads.get()} threads, '
            f'qualified quality phred {self.qualifiedQualityPhred.get()}, '
            f'unqualified percent limit {self.unqualifiedPercentLimit.get()}%, '
            f'N base limit {self.nBaseLimit.get()}, '
            f'minimum average quality {self.averageQual.get()}, '
            f'and adapter trimming '
            f'{"disabled" if self.disableAdapterTrimming.get() else "enabled"}.'
        )

        outputFastq = self.outputFastq
        methods.append(
            f'The filtered FASTQ dataset contains '
            f'{outputFastq.getNumReads()} '
            f'{"read pairs" if outputFastq.isPaired() else "reads"} with an average read '
            f'length of {outputFastq.getReadLength()} bp.'
        )

        if self.runFastqc.get() and self.outputFastq.getNumReads() > 0:
            methods.append(
                'FastQC quality control was performed on the filtered FASTQ files.'
            )

        return methods
