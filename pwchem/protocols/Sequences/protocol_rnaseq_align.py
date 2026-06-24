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

from pyworkflow.protocol.params import PointerParam,EnumParam,IntParam,FileParam,BooleanParam
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import AlignmentFile


class ProtRNASeqAlignment(EMProtocol):
    """Align FASTQ reads to a reference genome using STAR or HISAT2."""

    ALIGN_STAR = 0
    ALIGN_HISAT2 = 1

    REF_UPLOAD = 0
    REF_DOWNLOAD = 1

    SPECIES_HUMAN = 0
    SPECIES_MOUSE = 1

    _label = 'RNA-seq alignment'

    def _defineParams(self, form):
        form.addSection(label='Input FASTQ')

        form.addParam(
            'inputFastq',
            PointerParam,
            pointerClass='FastqFile',
            label='Input FASTQ: '
        )

        form.addSection(label='Alignment')

        form.addParam(
            'aligner',
            EnumParam,
            choices=['STAR', 'HISAT2'],
            default=self.ALIGN_STAR,
            label='Aligner: '
        )

        form.addParam(
            'threads',
            IntParam,
            default=4,
            label='Threads: '
        )

        form.addParam(
            'keepIntermediateFiles',
            BooleanParam,
            default=False,
            label='Keep intermediate files: '
        )

        form.addSection(label='Reference genome: ')

        form.addParam(
            'referenceMode',
            EnumParam,
            choices=['Upload reference genome', 'Download by species'],
            default=self.REF_UPLOAD,
            label='Reference genome source: '
        )

        form.addParam(
            'genomeFasta',
            FileParam,
            label='Genome FASTA: ',
            condition='referenceMode == %d' % self.REF_UPLOAD
        )

        form.addParam(
            'annotationGtf',
            FileParam,
            label='Annotation GTF: ',
            allowsNull=True,
            condition='referenceMode == %d' % self.REF_UPLOAD
        )

        form.addParam(
            'species',
            EnumParam,
            choices=['Homo sapiens', 'Mus musculus'],
            default=self.SPECIES_HUMAN,
            label='Species: ',
            condition='referenceMode == %d' % self.REF_DOWNLOAD
        )

        form.addParam(
            'downloadAnnotation',
            BooleanParam,
            default=True,
            label='Download annotation GTF: ',
            condition='referenceMode == %d' % self.REF_DOWNLOAD
        )

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareReferenceStep)
        self._insertFunctionStep(self.buildIndexStep)
        self._insertFunctionStep(self.alignStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------------------------------------------------------
    # Reference preparation
    # -------------------------------------------------------------------------

    def prepareReferenceStep(self):
        if self.referenceMode.get() == self.REF_UPLOAD:
            self._prepareUploadedReference()
        else:
            self._downloadReference()

    def _prepareUploadedReference(self):
        fastaPath = self.genomeFasta.get()
        gtfPath = self.annotationGtf.get() if self.annotationGtf.hasValue() else None

        self._fastaFile = os.path.abspath(fastaPath)
        self._gtfFile = os.path.abspath(gtfPath) if gtfPath else None
        self._referenceName = 'custom'

    def _downloadReference(self):
        referenceDir = self._getExtraPath('reference')
        os.makedirs(referenceDir, exist_ok=True)

        fastaUrl, gtfUrl, referenceName = self._getReferenceUrls()

        fastaGz = os.path.join(referenceDir, os.path.basename(fastaUrl))
        fastaFile = fastaGz.replace('.gz', '')

        if not os.path.exists(fastaFile):
            Plugin.runCondaCommand(
                self,
                '-L -o {} {}'.format(fastaGz, fastaUrl),
                RNASEQ_DIC,
                'curl'
            )

            Plugin.runCondaCommand(
                self,
                '-f {}'.format(fastaGz),
                RNASEQ_DIC,
                'gunzip'
            )

        self._fastaFile = fastaFile
        self._gtfFile = None
        self._referenceName = referenceName

        if self.downloadAnnotation.get() and gtfUrl:
            gtfGz = os.path.join(referenceDir, os.path.basename(gtfUrl))
            gtfFile = gtfGz.replace('.gz', '')

            if not os.path.exists(gtfFile):
                Plugin.runCondaCommand(
                    self,
                    '-L -o {} {}'.format(gtfGz, gtfUrl),
                    RNASEQ_DIC,
                    'curl'
                )

                Plugin.runCondaCommand(
                    self,
                    '-f {}'.format(gtfGz),
                    RNASEQ_DIC,
                    'gunzip'
                )

            self._gtfFile = gtfFile

    def _getReferenceUrls(self):
        if self.species.get() == self.SPECIES_HUMAN:
            return (
                'https://ftp.ensembl.org/pub/current_fasta/'
                'homo_sapiens/dna/'
                'Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz',
                'https://ftp.ensembl.org/pub/current_gtf/'
                'homo_sapiens/'
                'Homo_sapiens.GRCh38.116.gtf.gz',
                'Homo sapiens GRCh38'
            )

        return (
            'https://ftp.ensembl.org/pub/current_fasta/'
            'mus_musculus/dna/'
            'Mus_musculus.GRCm39.dna.primary_assembly.fa.gz',
            'https://ftp.ensembl.org/pub/current_gtf/'
            'mus_musculus/'
            'Mus_musculus.GRCm39.116.gtf.gz',
            'Mus musculus GRCm39'
        )

    # -------------------------------------------------------------------------
    # Index generation
    # -------------------------------------------------------------------------

    def buildIndexStep(self):
        if self.aligner.get() == self.ALIGN_STAR:
            self._buildStarIndex()
        else:
            self._buildHisat2Index()

    def _buildStarIndex(self):
        indexDir = self._getExtraPath('star_index')
        os.makedirs(indexDir, exist_ok=True)

        genomeParameters = os.path.join(indexDir, 'genomeParameters.txt')

        if os.path.exists(genomeParameters):
            self._indexPath = indexDir
            return

        args = [
            '--runMode genomeGenerate',
            '--genomeDir {}'.format(indexDir),
            '--genomeFastaFiles {}'.format(self._fastaFile),
            '--runThreadN {}'.format(self.threads.get())
        ]

        if self._gtfFile:
            args.append('--sjdbGTFfile {}'.format(self._gtfFile))

        Plugin.runCondaCommand(
            self,
            ' '.join(args),
            RNASEQ_DIC,
            'STAR'
        )

        self._indexPath = indexDir

    def _buildHisat2Index(self):
        indexDir = self._getExtraPath('hisat2_index')
        os.makedirs(indexDir, exist_ok=True)

        indexPrefix = os.path.join(indexDir, 'genome')
        indexCheckFile = indexPrefix + '.1.ht2'

        if os.path.exists(indexCheckFile):
            self._indexPath = indexPrefix
            return

        Plugin.runCondaCommand(
            self,
            '{} {}'.format(self._fastaFile, indexPrefix),
            RNASEQ_DIC,
            'hisat2-build'
        )

        self._indexPath = indexPrefix

    # -------------------------------------------------------------------------
    # Alignment
    # -------------------------------------------------------------------------

    def alignStep(self):
        if self.aligner.get() == self.ALIGN_STAR:
            self._alignWithStar()
        else:
            self._alignWithHisat2()

    def _alignWithStar(self):
        fastqObj = self.inputFastq.get()

        outPrefix = self._getExtraPath('star_')

        readFiles = [fastqObj.getFileName()]

        if fastqObj.isPaired():
            readFiles.append(fastqObj.getFileName2())

        args = [
            '--genomeDir {}'.format(self._indexPath),
            '--readFilesIn {}'.format(' '.join(readFiles)),
            '--runThreadN {}'.format(self.threads.get()),
            '--outSAMtype BAM SortedByCoordinate',
            '--outFileNamePrefix {}'.format(outPrefix)
        ]

        if fastqObj.isCompressed():
            args.append('--readFilesCommand zcat')

        Plugin.runCondaCommand(
            self,
            ' '.join(args),
            RNASEQ_DIC,
            'STAR'
        )

        bamFile = self._getExtraPath('star_Aligned.sortedByCoord.out.bam')
        baiFile = bamFile + '.bai'

        Plugin.runCondaCommand(
            self,
            'index {}'.format(bamFile),
            RNASEQ_DIC,
            'samtools'
        )

        self._alignmentFile = bamFile
        self._alignmentFormat = 'BAM'
        self._alignerName = 'STAR'
        self._logFile = self._getExtraPath('star_Log.final.out')
        self._isSorted = True
        self._isIndexed = True
        self._indexFile = baiFile

    def _alignWithHisat2(self):
        fastqObj = self.inputFastq.get()

        bamFile = self._getExtraPath('hisat2_alignment.sorted.bam')
        baiFile = bamFile + '.bai'
        logFile = self._getExtraPath('hisat2_alignment.log')

        if fastqObj.isPaired():
            args = (
                '-x {} -1 {} -2 {} -p {} 2> {} | '
                'samtools sort -@ {} -o {}'
            ).format(
                self._indexPath,
                fastqObj.getFileName(),
                fastqObj.getFileName2(),
                self.threads.get(),
                logFile,
                self.threads.get(),
                bamFile
            )
        else:
            args = (
                '-x {} -U {} -p {} 2> {} | '
                'samtools sort -@ {} -o {}'
            ).format(
                self._indexPath,
                fastqObj.getFileName(),
                self.threads.get(),
                logFile,
                self.threads.get(),
                bamFile
            )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'hisat2'
        )

        Plugin.runCondaCommand(
            self,
            'index {}'.format(bamFile),
            RNASEQ_DIC,
            'samtools'
        )

        self._alignmentFile = bamFile
        self._alignmentFormat = 'BAM'
        self._alignerName = 'HISAT2'
        self._logFile = logFile
        self._isSorted = True
        self._isIndexed = True
        self._indexFile = baiFile

    # -------------------------------------------------------------------------
    # Output
    # -------------------------------------------------------------------------

    def createOutputStep(self):
        fastqObj = self.inputFastq.get()

        outputAlignment = AlignmentFile(filename=self._alignmentFile)
        outputAlignment.setFormat(self._alignmentFormat)
        outputAlignment.setSampleName(fastqObj.getSampleName())
        outputAlignment.setAligner(self._alignerName)
        outputAlignment.setReferenceGenome(self._referenceName)
        outputAlignment.setIsSorted(self._isSorted)
        outputAlignment.setIsIndexed(self._isIndexed)
        outputAlignment.setIndexFile(self._indexFile)
        outputAlignment.setLogFile(self._logFile)

        self._defineOutputs(outputAlignment=outputAlignment)
        self._defineSourceRelation(self.inputFastq, outputAlignment)

    # -------------------------------------------------------------------------
    # Validation
    # -------------------------------------------------------------------------

    def _validate(self):
        errors = []

        if self.referenceMode.get() == self.REF_UPLOAD:
            fastaPath = self.genomeFasta.get()

            if not fastaPath or not os.path.exists(fastaPath):
                errors.append('A valid genome FASTA file is required.')

            if self.annotationGtf.hasValue():
                gtfPath = self.annotationGtf.get()

                if gtfPath and not os.path.exists(gtfPath):
                    errors.append('The provided GTF file does not exist.')

        fastqObj = self.inputFastq.get()

        if fastqObj.isPaired() and not fastqObj.hasFileName2():
            errors.append('Paired-end input requires read 2 file.')

        return errors

    # -------------------------------------------------------------------------
    # Summary / methods
    # -------------------------------------------------------------------------

    def _summary(self):
        if not hasattr(self, '_alignmentFile'):
            return []

        return [
            'Aligner: {}'.format(self._alignerName),
            'Reference: {}'.format(self._referenceName),
            'Output BAM: {}'.format(self._alignmentFile),
            'Index BAI: {}'.format(self._indexFile),
            'Sorted: {}'.format(self._isSorted),
            'Indexed: {}'.format(self._isIndexed)
        ]

    def _methods(self):
        return [
            'RNA-seq reads were aligned to the selected reference genome using '
            '{}. The resulting alignment was stored as a sorted and indexed BAM '
            'file.'.format(self._alignerName)
        ]