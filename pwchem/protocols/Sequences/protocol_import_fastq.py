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
    not copied or modified by the protocol. Paths stored in the output are
    absolute, so downstream protocols can access them from their working directories.


    ---------------------------------------------------------------------
    Input parameters
    ---------------------------------------------------------------------

    sampleName : str, optional
        Name used to identify the sample.

        If no name is provided, it is inferred automatically from the
        first FASTQ filename. Common FASTQ extensions and read suffixes
        such as ``_R1``, ``_R2``, ``_1`` and ``_2`` are removed when
        generating the default sample name. If this removes the entire
        name, ``sample`` is used as a fallback.

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
        number of reads. The stored count is the number of read pairs,
        rather than the sum of reads in both files.

    Read length
        The mean read length is calculated directly from the FASTQ data.

        For single-end datasets, the mean read length of the input file is
        stored.

        For paired-end datasets, the mean read lengths of read 1 and read 2
        are calculated independently and their mean value is stored in the
        output ``FastqFile``. Each file mean and the final pair mean are
        rounded to the nearest integer using Python rounding.

    Compression
        Single-end data are marked as compressed when the input filename
        ends in ``.gz``.

        Paired-end data are marked as compressed only when both input
        FASTQ files are compressed. Both mates must use the same compression
        mode, because downstream STAR alignment uses one reader for both files.
        Compression is inferred from the extension; the file contents are read
        during import to detect unreadable or invalid gzip data.


    ---------------------------------------------------------------------
    Workflow
    ---------------------------------------------------------------------

    The protocol performs the following steps:

    1. Validate input paths, extensions and paired-end file compatibility.

    2. Read each input FASTQ in streaming mode and calculate its number
       of reads and mean read length. For paired-end data, verify matching
       counts and calculate the mean of the two rounded read lengths.

    3. Determine the sample name from the user-provided value or infer it
       from the first FASTQ filename, using ``sample`` if no name remains.

    4. If requested, execute FastQC on the input FASTQ file(s).

    5. Create a ``FastqFile`` containing the input paths, sequencing
       metadata and available quality-control reports.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - A read 1 FASTQ file has been provided.
    - Input files use a supported FASTQ extension.
    - Each input is an existing, non-empty regular file.
    - Read 2 is provided for paired-end datasets.
    - Read 1 and read 2 are different files, including symbolic and hard links.
    - Paired-end inputs use the same compression mode.

    During import, paired-end datasets are additionally checked to ensure
    that read 1 and read 2 contain the same number of reads. FASTQ records
    are read in streaming mode and must have four lines, with an @ header,
    a + separator and matching sequence/quality lengths. Empty read sets,
    truncated records and unreadable compressed data stop execution before
    an output is published. Wrapped FASTQ records are not supported.

    Matching counts do not establish that read identifiers correspond or
    that the mates are in the same order; the input pair must already be
    synchronized. When FastQC is enabled, input filenames must produce
    distinct report names, even when located in separate directories.


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

    def _getInputFiles(self):
        files = [self.inputFastq1.get()]
        if self.isPaired.get():
            files.append(self.inputFastq2.get())
        return [os.path.abspath(path) if path else '' for path in files]

    def importStep(self):
        errors = self._validate()
        if errors:
            raise RuntimeError('Invalid FASTQ inputs: ' + '; '.join(errors))

        files = self._getInputFiles()
        stats = [getFastqStats(path) for path in files]
        for path, (numReads, _) in zip(files, stats):
            if numReads == 0:
                raise RuntimeError('FASTQ contains no reads: {}'.format(path))
        numReads, readLength = stats[0]
        if len(stats) == 2:
            numReads2, readLength2 = stats[1]
            if numReads != numReads2:
                raise RuntimeError('Paired FASTQ files have different number of reads: '
                                   'R1={} R2={}'.format(numReads, numReads2))
            readLength = int(round((readLength + readLength2) / 2))

        sampleName = (self.sampleName.get() or '').strip() or self._getDefaultSampleName(files[0])
        fastq = FastqFile()
        fastq.setFileName(files[0])
        fastq.setSampleName(sampleName)
        fastq.setIsPaired(self.isPaired.get())
        fastq.setIsCompressed(files[0].endswith('.gz'))
        fastq.setFormat('FASTQ')
        fastq.setHasQuality(True)
        fastq.setNumReads(numReads)
        fastq.setReadLength(readLength)
        if len(files) == 2:
            fastq.setFileName2(files[1])

        if self.runFastqc.get():
            htmlFiles = self._runFastqc()
            if self.isPaired.get():
                fastq.setFastqcHtmlR1(htmlFiles[0])
                fastq.setFastqcHtmlR2(htmlFiles[1])
            else:
                fastq.setFastqcHtml(htmlFiles[0])

        self._defineOutputs(outputFastq=fastq)

    def _runFastqc(self):
        return runFastqc(self, self._getInputFiles())

    @staticmethod
    def _getDefaultSampleName(fn):
        sampleName = os.path.basename(fn)
        for extension in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if sampleName.endswith(extension):
                sampleName = sampleName[:-len(extension)]
                break
        for suffix in ('_R1', '_R2', '_1', '_2', '.R1', '.R2', '.1', '.2'):
            if sampleName.endswith(suffix):
                sampleName = sampleName[:-len(suffix)]
                break
        return sampleName or 'sample'

    def _validate(self):
        errors = []
        files = self._getInputFiles()
        validExtensions = ('.fastq', '.fq', '.fastq.gz', '.fq.gz')
        for number, path in enumerate(files, start=1):
            label = 'Read {}'.format(number)
            if not path:
                errors.append('{} is required{}.'.format(
                    label, ' for paired-end data' if number == 2 else ''))
                continue
            if not path.endswith(validExtensions):
                errors.append('{} must be a FASTQ file.'.format(label))
            if not os.path.isfile(path):
                errors.append('{} file does not exist: {}'.format(label, path))
            elif os.path.getsize(path) == 0:
                errors.append('{} file is empty: {}'.format(label, path))

        if len(files) == 2 and all(files):
            if all(os.path.isfile(path) for path in files) and os.path.samefile(*files):
                errors.append('Read 1 and read 2 must be different files.')
            if files[0].endswith('.gz') != files[1].endswith('.gz'):
                errors.append('Paired FASTQ files must use the same compression mode.')
        return errors

    def _summary(self):
        if not hasattr(self, 'outputFastq'):
            return ['Output FASTQ not available yet.']

        fastq = self.outputFastq

        summary = [
            'Format: FASTQ',
            'Quality scores: yes',
            f'Sample name: {fastq.getSampleName()}',
            f'Read 1: {fastq.getFileName()}',
            f'Number of reads: {fastq.getNumReads()}'
        ]

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
        if not hasattr(self, 'outputFastq'):
            return []
        fastq = self.outputFastq
        sequencingType = 'paired-end' if fastq.isPaired() else 'single-end'
        methods = [
            'A {} FASTQ dataset was imported.'.format(sequencingType),
            'Sample name: {}.'.format(fastq.getSampleName())
        ]
        if self.runFastqc.get():
            methods.append('FastQC quality control was performed.')
        return methods
