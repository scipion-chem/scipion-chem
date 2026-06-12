# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
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

    Optionally, FastQC can be executed to assess read quality and produce
    HTML reports.

    ----------------------------
    Input parameters
    ----------------------------

    sampleName : str, optional
        Identifier assigned to the sample.

    isPaired : bool
        If True, the input is treated as paired-end data and two FASTQ files are required.

    inputFastq1 : file
        FASTQ file corresponding to read 1 (or single-end data).

    inputFastq2 : file, optional
        FASTQ file corresponding to read 2 (only required if isPaired=True).

    runFastqc : bool
        If True, FastQC is executed and quality control reports are generated.

    ----------------------------
    Output
    ----------------------------

    outputFastq : FastqFile
        Imported FASTQ dataset ready for downstream analysis.

    FastQC HTML reports are stored as attributes of the outputFastq object.

    ----------------------------
    Notes
    ----------------------------

    - FastQC must be installed through the plugin binaries.
    - Generated HTML reports are stored in the protocol's extra directory.
    """

    _label = 'import fastq'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('sampleName', StringParam,
                      label='Sample name',
                      allowsNull=True,
                      help='Name used to identify the sample. If empty, it will be inferred from the FASTQ file name.')

        form.addParam('isPaired', BooleanParam,
                      default=False,
                      label='Paired-end',
                      help='Select Yes if paired-end sequencing.')

        form.addParam('inputFastq1', FileParam,
                      label='FASTQ file (read 1 / single-end)',
                      allowsNull=False)

        form.addParam('inputFastq2', FileParam,
                      condition='isPaired',
                      label='FASTQ file (read 2)',
                      allowsNull=False)

        form.addParam('runFastqc', BooleanParam,
                      default=False,
                      label='Run FastQC:',
                      help='Execute FastQC and generate HTML reports.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.importStep)
        if self.runFastqc.get():
            self._insertFunctionStep(self.fastqcStep)

    def importStep(self):
        fn1 = self.inputFastq1.get()

        fastq = FastqFile()
        fastq.setFileName(fn1)
        fastq.setIsPaired(self.isPaired.get())
        fastq.setIsCompressed(fn1.endswith('.gz'))
        fastq.setFormat('FASTQ')
        fastq.setHasQuality(True)

        sample = self.sampleName.get()
        sampleName = sample.strip() if sample and sample.strip() else self._getDefaultSampleName(fn1)
        fastq.setSampleName(sampleName)

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            fastq.setFileName2(fn2)
            fastq.setIsCompressed(fn1.endswith('.gz') and fn2.endswith('.gz'))

        self._defineOutputs(outputFastq=fastq)

    def fastqcStep(self):
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
            raise RuntimeError('FastQC finished but no HTML report was generated.')

        fastq = self.outputFastq
        if self.isPaired.get():
            if len(htmlFiles) < 2:
                raise RuntimeError('Expected two FastQC HTML reports for paired-end data.')

            fastq.setFastqcHtmlR1(htmlFiles[0])
            fastq.setFastqcHtmlR2(htmlFiles[1])
        else:
            fastq.setFastqcHtml(htmlFiles[0])

        self._store(fastq)

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
        sampleName = sample.strip() if sample and sample.strip() else self._getDefaultSampleName(fn1)

        summary.append('Format: FASTQ')
        summary.append('Quality scores: yes')


        summary.append(f'Sample name: {sampleName}')
        summary.append(f'Read 1: {fn1}')

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            summary.append(f'Read 2: {fn2}')
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
        sampleName = sample.strip() if sample and sample.strip() else self._getDefaultSampleName(fn1)

        methods.append(f'Sample name: {sampleName}.')

        if self.runFastqc.get():
            methods.append('FastQC quality control was performed.')

        return methods