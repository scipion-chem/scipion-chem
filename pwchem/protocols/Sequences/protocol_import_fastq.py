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
from pyworkflow.protocol.params import BooleanParam, FileParam, StringParam
from pwchem.objects import FastqFile
from pwchem.utils.utilsRNA import getFastqStats, runFastqc


class ProtImportFastq(EMProtocol):
    """
    Import FASTQ files from RNA sequencing experiments.

    This protocol imports single-end or paired-end FASTQ files into
    Scipion and creates a ``FastqFile`` object that can be used as input
    for downstream RNA-seq protocols.

    Basic sequencing metadata are calculated directly from the input
    FASTQ file(s), including the number of reads and mean read length.
    Optionally, FastQC can be executed to generate quality-control reports.

    The original FASTQ files are referenced by the output object and are
    not copied or modified by the protocol.


    ---------------------------------------------------------------------
    Input parameters
    ---------------------------------------------------------------------

    sampleName : str, optional
        Name used to identify the sample.

        If no name is provided, it is inferred automatically from the
        first FASTQ filename. Common FASTQ extensions and read suffixes
        such as ``_R1``, ``_R2``, ``_1`` and ``_2`` are removed when
        generating the default sample name.

    isPaired : bool
        Defines whether the sequencing dataset is single-end or paired-end.

        If True, both read 1 and read 2 FASTQ files must be provided.

    inputFastq1 : file
        FASTQ file containing read 1 for paired-end data or the reads for
        single-end data.

        Supported extensions are::

            .fastq
            .fq
            .fastq.gz
            .fq.gz

    inputFastq2 : file, optional
        FASTQ file containing read 2.

        This parameter is required only when ``isPaired`` is enabled.

        Supported extensions are the same as for ``inputFastq1``.

    runFastqc : bool
        If True, FastQC is executed on the imported FASTQ file(s).

        For single-end data, one FastQC report is generated.

        For paired-end data, independent FastQC reports are generated for
        read 1 and read 2.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputFastq : FastqFile
        Imported FASTQ dataset.

        The output object contains:

        - FASTQ file path(s).
        - Sample name.
        - Sequencing type (single-end or paired-end).
        - FASTQ format information.
        - Presence of quality scores.
        - Compression status.
        - Number of reads.
        - Mean read length.
        - FastQC HTML report(s), when requested.

        The output references the original input FASTQ files; the sequence
        data are not copied or modified during import.


    ---------------------------------------------------------------------
    FASTQ metadata
    ---------------------------------------------------------------------

    Number of reads
        The number of sequencing reads is calculated directly from the
        input FASTQ file.

        For paired-end datasets, read 1 and read 2 must contain the same
        number of reads.

    Read length
        The mean read length is calculated directly from the FASTQ data.

        For single-end datasets, the mean read length of the input file is
        stored.

        For paired-end datasets, the mean read lengths of read 1 and read 2
        are calculated independently and their mean value is stored in the
        output ``FastqFile``.

    Compression
        Single-end data are marked as compressed when the input filename
        ends in ``.gz``.

        Paired-end data are marked as compressed only when both input
        FASTQ files are compressed.


    ---------------------------------------------------------------------
    Workflow
    ---------------------------------------------------------------------

    The protocol performs the following steps:

    1. Read the input FASTQ file and calculate the number of reads and
       mean read length.

    2. Determine the sample name from the user-provided value or infer it
       from the FASTQ filename.

    3. For paired-end datasets, calculate the statistics of read 2 and
       verify that both files contain the same number of reads.

    4. If requested, execute FastQC on the input FASTQ file(s).

    5. Create a ``FastqFile`` containing the input paths, sequencing
       metadata and available quality-control reports.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - A read 1 FASTQ file has been provided.
    - Input files use a supported FASTQ extension.
    - Read 2 is provided for paired-end datasets.
    - Read 1 and read 2 are different files for paired-end datasets.

    During import, paired-end datasets are additionally checked to ensure
    that read 1 and read 2 contain the same number of reads.


    ---------------------------------------------------------------------
    Requirements
    ---------------------------------------------------------------------

    FastQC must be available through the RNA-seq environment configured by
    the Scipion-Chem plugin when ``runFastqc`` is enabled.


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Both single-end and paired-end FASTQ datasets are supported.

    - Compressed and uncompressed FASTQ files are supported.

    - Supported extensions are ``.fastq``, ``.fq``, ``.fastq.gz`` and
      ``.fq.gz``.

    - The input FASTQ files are not copied or modified.

    - If no sample name is provided, it is inferred from the first FASTQ
      filename.

    - Number of reads and mean read length are calculated directly from
      the FASTQ file(s).

    - Paired-end FASTQ files must contain the same number of reads.

    - For paired-end datasets, the stored read length corresponds to the
      mean of the R1 and R2 mean read lengths.

    - FastQC is optional and does not modify the imported sequencing data.

    - Generated FastQC HTML reports are stored as attributes of the output
      ``FastqFile`` object.
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
        numReads, readLength = getFastqStats(fn1)

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
            numReads2, readLength2 = getFastqStats(fn2)

            if numReads != numReads2:
                raise RuntimeError(
                    'Paired FASTQ files have different number of reads: '
                    'R1={} R2={}'.format(numReads, numReads2)
                )

            fastq.setFileName2(fn2)
            fastq.setIsCompressed(
                fn1.endswith('.gz') and fn2.endswith('.gz')
            )

            readLength = int(round((readLength + readLength2) / 2))
            fastq.setReadLength(readLength)

        if self.runFastqc.get():
            htmlFiles = self._runFastqc()

            if self.isPaired.get():
                fastq.setFastqcHtmlR1(htmlFiles[0])
                fastq.setFastqcHtmlR2(htmlFiles[1])
            else:
                fastq.setFastqcHtml(htmlFiles[0])

        self._defineOutputs(outputFastq=fastq)

    def _runFastqc(self):
        fastqFiles = [self.inputFastq1.get()]

        if self.isPaired.get():
            fastqFiles.append(self.inputFastq2.get())

        return runFastqc(self, fastqFiles)

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

    def _validate(self):
        errors = []
        validExt = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')

        fn1 = self.inputFastq1.get()

        if not fn1:
            errors.append('Read 1 is required.')
        elif not fn1.endswith(validExt):
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

        if not hasattr(self, 'outputFastq'):
            return ['Output FASTQ not available yet.']

        fastq = self.outputFastq

        summary.append('Format: FASTQ')
        summary.append('Quality scores: yes')
        summary.append(f'Sample name: {fastq.getSampleName()}')
        summary.append(f'Read 1: {fastq.getFileName()}')
        summary.append(f'Number of reads: {fastq.getNumReads()}')

        if fastq.getReadLength() > 0:
            summary.append(
                f'Mean read length: {fastq.getReadLength()} bp'
            )

        if fastq.isPaired():
            summary.append(f'Read 2: {fastq.getFileName2()}')
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