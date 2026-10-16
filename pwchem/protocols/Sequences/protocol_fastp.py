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

import os

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
          Minimum read length required after filtering.

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
          - FastQC HTML report(s), when requested.

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

      6. If requested, execute FastQC on the filtered reads.

      7. Create the output ``FastqFile`` containing the filtered files,
         statistics and quality-control reports.


      ---------------------------------------------------------------------
      Notes
      ---------------------------------------------------------------------

      - Both single-end and paired-end FASTQ datasets are supported.

      - Input files must contain FASTQ quality scores.

      - The minimum read length cannot exceed the input read length when
        this information is available.

      - Filtered FASTQ files are written in uncompressed FASTQ format.

      - Output filenames are generated using the input sample name.

      - If the input sample does not have a sample name, ``sample`` is used.

      - For paired-end datasets, R1 and R2 are filtered together.

      - The protocol verifies that filtered paired-end files contain the
        same number of reads.

      - For paired-end datasets, the stored read length corresponds to the
        mean of the R1 and R2 mean read lengths.

      - FastQC, when enabled, is executed on the filtered reads rather than
        on the original input files.

      - fastp HTML and JSON reports are always generated.

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
                      default=15,
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
                      default=10,
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
        inputFastq = self.inputFastq.get()

        outputFastq = self._runFastp(inputFastq)

        if self.runFastqc.get():
            htmlFiles = self._runFastqc(outputFastq)

            if outputFastq.isPaired():
                outputFastq.setFastqcHtmlR1(htmlFiles[0])
                outputFastq.setFastqcHtmlR2(htmlFiles[1])
            else:
                outputFastq.setFastqcHtml(htmlFiles[0])

        self._defineOutputs(outputFastq=outputFastq)

    def _runFastp(self, inputFastq):
        fn1 = inputFastq.getFileName()
        isPaired = inputFastq.isPaired()

        lengthRequired = self.lengthRequired.get()
        threads = self.numberOfThreads.get()
        qualifiedQualityPhred = self.qualifiedQualityPhred.get()
        unqualifiedPercentLimit = self.unqualifiedPercentLimit.get()
        nBaseLimit = self.nBaseLimit.get()
        averageQual = self.averageQual.get()
        extraArgs = ''

        if self.disableAdapterTrimming.get():
            extraArgs += '-A '

        outDir = self._getExtraPath('fastp')
        os.makedirs(outDir, exist_ok=True)

        sampleName = inputFastq.getSampleName()

        if not sampleName:
            sampleName = 'sample'

        out1 = os.path.join(
            outDir,
            f'{sampleName}_filtered_R1.fastq'
        )

        reportHtml = os.path.join(
            outDir,
            f'{sampleName}_fastp.html'
        )

        reportJson = os.path.join(
            outDir,
            f'{sampleName}_fastp.json'
        )

        arguments = (
            f'-i "{fn1}" '
            f'-o "{out1}" '
            f'-w {threads} '
            f'-l {lengthRequired} '
            f'-q {qualifiedQualityPhred} '
            f'-u {unqualifiedPercentLimit} '
            f'-n {nBaseLimit} '
            f'-e {averageQual} '
            f'{extraArgs}'
            f'-h "{reportHtml}" '
            f'-j "{reportJson}"'
        )

        if isPaired:
            fn2 = inputFastq.getFileName2()

            out2 = os.path.join(
                outDir,
                f'{sampleName}_filtered_R2.fastq'
            )

            arguments = (
                f'-i "{fn1}" '
                f'-I "{fn2}" '
                f'-o "{out1}" '
                f'-O "{out2}" '
                f'-w {threads} '
                f'-l {lengthRequired} '
                f'-q {qualifiedQualityPhred} '
                f'-u {unqualifiedPercentLimit} '
                f'-n {nBaseLimit} '
                f'-e {averageQual} '
                f'{extraArgs}'
                f'-h "{reportHtml}" '
                f'-j "{reportJson}"'
            )

        Plugin.runCondaCommand(self, arguments, RNASEQ_DIC, 'fastp')

        if not os.path.exists(out1):
            raise RuntimeError(
                f'fastp did not generate the expected output FASTQ: {out1}'
            )

        if isPaired and not os.path.exists(out2):
            raise RuntimeError(
                f'fastp did not generate the expected read 2 FASTQ: {out2}'
            )

        if not os.path.exists(reportHtml):
            raise RuntimeError(
                f'fastp did not generate the expected HTML report: {reportHtml}'
            )

        if not os.path.exists(reportJson):
            raise RuntimeError(
                f'fastp did not generate the expected JSON report: {reportJson}'
            )

        numReads, readLength = getFastqStats(out1)

        if isPaired:
            numReads2, readLength2 = getFastqStats(out2)

            if numReads != numReads2:
                raise RuntimeError(
                    'Filtered paired FASTQ files have different number of reads: '
                    'R1={} R2={}'.format(numReads, numReads2)
                )

            readLength = int(round((readLength + readLength2) / 2))

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
            errors.append('An input FASTQ object is required.')
            return errors

        if inputFastq.isPaired() and not inputFastq.hasFileName2():
            errors.append(
                'Input FASTQ is marked as paired-end but read 2 is missing.'
            )

        if not inputFastq.hasQuality():
            errors.append(
                'Input file does not contain quality scores and cannot be '
                'processed as FASTQ.'
            )

        if self.lengthRequired.get() < 0:
            errors.append('Minimum read length must be greater than or equal to 0.')

        if inputFastq.getReadLength() > 0:
            if self.lengthRequired.get() > inputFastq.getReadLength():
                errors.append(
                    f'Minimum read length ({self.lengthRequired.get()}) is larger '
                    f'than the input read length ({inputFastq.getReadLength()}).'
                )

        if not 0 <= self.qualifiedQualityPhred.get() <= 36:
            errors.append('Qualified quality phred must be between 0 and 36.')

        if not 0 <= self.unqualifiedPercentLimit.get() <= 100:
            errors.append('Unqualified percent limit must be between 0 and 100.')

        if self.nBaseLimit.get() < 0:
            errors.append('N base limit must be greater than or equal to 0.')

        if self.averageQual.get() < 0:
            errors.append('Average quality must be greater than or equal to 0.')

        return errors

    def _summary(self):
        summary = []

        inputFastq = self.inputFastq.get()

        summary.append('FASTQ filtering performed with fastp')
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
            if inputFastq.getReadLength() > 0:
                summary.append(
                    f'Input read length: {inputFastq.getReadLength()} bp'
                )
            if inputFastq.getNumReads() > 0:
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

            if outputFastq.getNumReads() > 0:
                summary.append(
                    f'Output number of reads: {outputFastq.getNumReads()}'
                )

        summary.append(
            f'FastQC: {"executed" if self.runFastqc.get() else "not executed"}'
        )

        return summary

    def _methods(self):
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

        if hasattr(self, 'outputFastq'):
            outputFastq = self.outputFastq

            methods.append(
                f'The filtered FASTQ dataset contains '
                f'{outputFastq.getNumReads()} reads with an average read '
                f'length of {outputFastq.getReadLength()} bp.'
            )

        if self.runFastqc.get():
            methods.append(
                'FastQC quality control was performed on the filtered FASTQ files.'
            )

        return methods