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
import gzip

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import BooleanParam, FileParam, StringParam

from pwchem import Plugin
from pwchem.constants import OPENBABEL_DIC
from pwchem.objects import FastqFile


class ProtImportFastq(EMProtocol):
    """
    Import FASTQ files from RNA sequencing experiments.

    This protocol imports single-end or paired-end FASTQ files into Scipion
    and generates a FastqFile object for downstream processing.

    Basic FASTQ metadata, including sequencing type, number of reads and mean
    read length, are extracted from the input FASTQ file(s) and stored in the
    output object for downstream analysis.

    Optionally, FastQC can be executed to assess read quality and generate
    HTML reports.

    ----------------------------
    Input parameters
    ----------------------------

    sampleName : str, optional
        Identifier assigned to the sample.

    isPaired : bool
        If True, the input is treated as paired-end data and two FASTQ files
        are required.

    inputFastq1 : file
        FASTQ file corresponding to read 1 or single-end data.

    inputFastq2 : file, optional
        FASTQ file corresponding to read 2. Only required for paired-end data.

    runFastqc : bool
        If True, FastQC is executed and quality control reports are generated.

    ----------------------------
    Output
    ----------------------------

    outputFastq : FastqFile
        Imported FASTQ dataset ready for downstream analysis.

        The output object contains:

        - FASTQ file path(s).
        - Sample name.
        - Sequencing type (single-end or paired-end).
        - Number of reads.
        - Mean read length.
        - FastQC HTML report, if generated.
        - FastQC HTML reports for read 1 and read 2, if paired-end.

    ----------------------------
    Notes
    ----------------------------

    - FastQC must be installed through the plugin binaries.
    - Number of reads and mean read length are calculated directly from the
      input FASTQ file(s).
    - For paired-end datasets, read counts are validated to ensure that read 1
      and read 2 contain the same number of reads.
    - The reported read length corresponds to the mean read length. For paired-end
      datasets, the mean value from read 1 and read 2 is stored.
    - Generated HTML reports are stored as attributes of the output FastqFile
      object.
    """

    _label = 'import fastq'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('sampleName', StringParam,
                      label='Sample name: ',
                      allowsNull=True,
                      help='Name used to identify the sample. If empty, '
                           'it will be inferred from the FASTQ file name.')

        form.addParam('isPaired', BooleanParam,
                      default=False,
                      label='Paired-end: ',
                      help='Select Yes if paired-end sequencing.')

        form.addParam('inputFastq1', FileParam,
                      label='FASTQ file (read 1 / single-end): ',
                      allowsNull=False)

        form.addParam('inputFastq2', FileParam,
                      condition='isPaired',
                      label='FASTQ file (read 2): ',
                      allowsNull=False)

        form.addParam('runFastqc', BooleanParam,
                      default=False,
                      label='Run FastQC: ',
                      help='Execute FastQC and generate HTML reports.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.importStep)

    def importStep(self):
        fn1 = self.inputFastq1.get()
        numReads, readLength = self._getFastqStats(fn1)

        fastq = FastqFile()
        fastq.setFileName(fn1)
        fastq.setIsPaired(self.isPaired.get())
        fastq.setIsCompressed(fn1.endswith('.gz'))
        fastq.setFormat('FASTQ')
        fastq.setHasQuality(True)
        fastq.setNumReads(numReads)
        fastq.setReadLength(readLength)

        sample = self.sampleName.get()
        sampleName = sample.strip() if sample and sample.strip() else \
            self._getDefaultSampleName(fn1)
        fastq.setSampleName(sampleName)

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            numReads2, readLength2 = self._getFastqStats(fn2)

            if numReads != numReads2:
                raise RuntimeError(
                    'Paired FASTQ files have different number of reads: '
                    'R1={} R2={}'.format(numReads, numReads2)
                )

            fastq.setFileName2(fn2)
            fastq.setIsCompressed(fn1.endswith('.gz') and fn2.endswith('.gz'))

            readLength = int(round((readLength + readLength2) / 2))
            fastq.setReadLength(readLength)

        if self.runFastqc.get():
            htmlFiles = self._runFastqc()

            if self.isPaired.get():
                if len(htmlFiles) < 2:
                    raise RuntimeError(
                        'Expected two FastQC HTML reports for paired-end data.'
                    )

                fastq.setFastqcHtmlR1(htmlFiles[0])
                fastq.setFastqcHtmlR2(htmlFiles[1])
            else:
                fastq.setFastqcHtml(htmlFiles[0])

        self._defineOutputs(outputFastq=fastq)

    def _openFastq(self, fn):
        if fn.endswith('.gz'):
            return gzip.open(fn, 'rt')

        return open(fn, 'r')

    def _getFastqStats(self, fn):
        """
        Calculate number of reads and mean read length from a FASTQ file.
        """
        numReads = 0
        totalLength = 0

        with self._openFastq(fn) as f:
            while True:
                header = f.readline()

                if not header:
                    break

                seq = f.readline().strip()
                f.readline()
                f.readline()

                numReads += 1
                totalLength += len(seq)

        readLength = int(round(totalLength / numReads)) if numReads else 0

        return numReads, readLength

    def _runFastqc(self):
        fn1 = self.inputFastq1.get()
        outDir = self._getExtraPath('fastqc')
        os.makedirs(outDir, exist_ok=True)

        arguments = f'-o "{outDir}" "{fn1}"'

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            arguments += f' "{fn2}"'

        Plugin.runCondaCommand(self, arguments, OPENBABEL_DIC, 'fastqc')

        htmlFiles = self._getFastqcHtmlFiles(outDir)

        if not htmlFiles:
            raise RuntimeError(
                'FastQC finished but no HTML report was generated.'
            )

        return htmlFiles

    def _getDefaultSampleName(self, fn):
        sampleName = os.path.basename(fn)

        for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
            if sampleName.endswith(ext):
                sampleName = sampleName[:-len(ext)]
                break

        for suffix in ['_R1', '_R2', '_1', '_2', '.R1', '.R2', '.1', '.2']:
            if sampleName.endswith(suffix):
                sampleName = sampleName[:-len(suffix)]
                break

        return sampleName

    def _getFastqcHtmlFiles(self, outDir):
        return sorted([
            os.path.join(outDir, f)
            for f in os.listdir(outDir)
            if f.endswith('.html')
        ])

    def _validate(self):
        errors = []
        validExt = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

        fn1 = self.inputFastq1.get()

        if fn1 and not fn1.endswith(validExt):
            errors.append('Read 1 must be a FASTQ file.')

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()

            if not fn2:
                errors.append('Read 2 is required for paired-end data.')
            elif not fn2.endswith(validExt):
                errors.append('Read 2 must be a FASTQ file.')
            elif fn1 == fn2:
                errors.append('Read 1 and read 2 must be different files.')

        return errors

    def _summary(self):
        summary = []

        sample = self.sampleName.get()
        fn1 = self.inputFastq1.get()
        sampleName = sample.strip() if sample and sample.strip() else \
            self._getDefaultSampleName(fn1)

        numReads, readLength = self._getFastqStats(fn1)

        summary.append('Format: FASTQ')
        summary.append('Quality scores: yes')
        summary.append(f'Sample name: {sampleName}')
        summary.append(f'Read 1: {fn1}')
        summary.append(f'Number of reads: {numReads}')

        if readLength > 0:
            summary.append(f'Mean read length: {readLength} bp')

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            numReads2, readLength2 = self._getFastqStats(fn2)

            summary.append(f'Read 2: {fn2}')
            summary.append(f'Number of reads R2: {numReads2}')

            if readLength2 > 0:
                summary.append(f'Mean read length R2: {readLength2} bp')

            summary.append('Sequencing type: paired-end')
        else:
            summary.append('Sequencing type: single-end')

        summary.append(
            f'FastQC: {"executed" if self.runFastqc.get() else "not executed"}'
        )

        return summary

    def _methods(self):
        methods = []

        if self.isPaired.get():
            methods.append('A paired-end FASTQ dataset was imported.')
        else:
            methods.append('A single-end FASTQ dataset was imported.')

        sample = self.sampleName.get()
        fn1 = self.inputFastq1.get()
        sampleName = sample.strip() if sample and sample.strip() else \
            self._getDefaultSampleName(fn1)

        methods.append(f'Sample name: {sampleName}.')

        if self.runFastqc.get():
            methods.append('FastQC quality control was performed.')

        return methods