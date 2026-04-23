# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perezliens@cnb.csic.es)
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

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import BooleanParam, FileParam, StringParam

from pwchem.objects import FastqFile


class ProtImportFastq(EMProtocol):
    _label = 'import fastq'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('sampleName', StringParam,
                      label='Sample name',
                      allowsNull=True)

        form.addParam('isPaired', BooleanParam,
                      default=False,
                      label='Paired-end')

        form.addParam('inputFastq1', FileParam,
                      label='FASTQ read 1 / single-end',
                      allowsNull=False)

        form.addParam('inputFastq2', FileParam,
                      condition='isPaired',
                      label='FASTQ read 2',
                      allowsNull=False)

    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        fn1 = self.inputFastq1.get()

        fastq = FastqFile()
        fastq.setFileName(fn1)
        fastq.setIsPaired(self.isPaired.get())

        if self.sampleName.get():
            fastq.setSampleName(self.sampleName.get())

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            fastq.setFileName2(fn2)
            fastq.setIsCompressed(fn1.endswith('.gz') and fn2.endswith('.gz'))
        else:
            fastq.setIsCompressed(fn1.endswith('.gz'))

        self._defineOutputs(outputFastq=fastq)

    def _validate(self):
        errors = []
        valid_ext = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

        fn1 = self.inputFastq1.get()
        if fn1 and not fn1.endswith(valid_ext):
            errors.append('Read 1 must be a FASTQ file.')

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            if not fn2:
                errors.append('Read 2 is required for paired-end data.')
            elif not fn2.endswith(valid_ext):
                errors.append('Read 2 must be a FASTQ file.')
            elif fn1 == fn2:
                errors.append('Read 1 and read 2 must be different files.')

        return errors

    def _summary(self):
        summary = []

        sample = self.sampleName.get()
        fn1 = self.inputFastq1.get()

        summary.append('Format: FASTQ')
        summary.append('Quality scores: yes')

        if sample:
            summary.append('Sample name: {}'.format(sample))

        summary.append('Read 1: {}'.format(fn1))

        if self.isPaired.get():
            fn2 = self.inputFastq2.get()
            summary.append('Read 2: {}'.format(fn2))
            summary.append('Sequencing type: paired-end')
            summary.append('Compressed: {}'.format(
                'yes' if fn1.endswith('.gz') and fn2.endswith('.gz') else 'no'))
        else:
            summary.append('Sequencing type: single-end')
            summary.append('Compressed: {}'.format(
                'yes' if fn1.endswith('.gz') else 'no'))

        return summary

    def _methods(self):
        methods = []

        if self.isPaired.get():
            methods.append('Imported a paired-end FASTQ dataset.')
        else:
            methods.append('Imported a single-end FASTQ dataset.')

        if self.sampleName.get():
            methods.append('Sample name: {}.'.format(self.sampleName.get()))

        return methods