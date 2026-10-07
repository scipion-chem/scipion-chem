# ***************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
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
# ***************************************************************************
import gzip
import json
import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import (
    ProtImportFastq,
    ProtImportGenomes,
    ProtRNASeqAlignment,
    ProtPicard
)


N_TEST_READS = 1000


class TestPicard(BaseTest):
    """Test Picard processing of RNA-seq BAM alignments."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')

        cls.fastqFile = cls.dataset.getFile(
            'mouse/SRR1552445.fastq.gz'
        )

        cls.genomeFile = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa'
        )

        cls.gtfFile = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.116.gtf'
        )

    # ---------------------------------------------------------------------
    # Input preparation
    # ---------------------------------------------------------------------

    @staticmethod
    def _openTextFile(fileName):
        """Open plain-text or gzip-compressed files in text mode."""
        if fileName.endswith('.gz'):
            return gzip.open(fileName, 'rt')

        return open(fileName, 'r')

    def _createSmallFastq(self):
        """Create a reduced FASTQ containing the first N_TEST_READS reads."""
        outputFastq = self.getOutputPath(
            'mouse_picard_test.fastq'
        )
        maxLines = N_TEST_READS * 4

        with self._openTextFile(self.fastqFile) as inputFile, \
                open(outputFastq, 'w') as outputFile:

            for lineNumber, line in enumerate(inputFile):
                if lineNumber >= maxLines:
                    break

                outputFile.write(line)

        self.assertTrue(
            os.path.isfile(outputFastq),
            'Reduced FASTQ was not created.'
        )

        self.assertGreater(
            os.path.getsize(outputFastq),
            0,
            'Reduced FASTQ is empty.'
        )

        return outputFastq

    def _createChromosomeReference(self):
        """Create a reduced reference containing one chromosome."""
        outputFasta = self.getOutputPath(
            'mouse_picard_test.fa'
        )
        outputGtf = self.getOutputPath(
            'mouse_picard_test.gtf'
        )

        chromosome = None

        with self._openTextFile(self.genomeFile) as inputFile, \
                open(outputFasta, 'w') as outputFile:

            for line in inputFile:

                if line.startswith('>'):
                    currentChromosome = line[1:].split()[0]

                    if chromosome is None:
                        chromosome = currentChromosome

                    elif currentChromosome != chromosome:
                        break

                outputFile.write(line)

        self.assertIsNotNone(
            chromosome,
            'No sequence was found in the reference FASTA.'
        )

        featureCount = 0

        with self._openTextFile(self.gtfFile) as inputFile, \
                open(outputGtf, 'w') as outputFile:

            for line in inputFile:

                if line.startswith('#'):
                    outputFile.write(line)
                    continue

                fields = line.rstrip().split('\t')

                if len(fields) < 9:
                    continue

                if fields[0] == chromosome:
                    outputFile.write(line)
                    featureCount += 1

        self.assertTrue(
            os.path.isfile(outputFasta),
            'Reduced FASTA was not created.'
        )

        self.assertGreater(
            os.path.getsize(outputFasta),
            0,
            'Reduced FASTA is empty.'
        )

        self.assertGreater(
            featureCount,
            0,
            'Reduced GTF does not contain annotations.'
        )

        return outputFasta, outputGtf

    def _createInputAlignment(self):
        """Create a small STAR alignment to use as Picard input."""
        smallFastq = self._createSmallFastq()

        referenceFasta, referenceGtf = (
            self._createChromosomeReference()
        )

        sampleName = 'mouse_picard_test'

        # -------------------------------------------------------------
        # Import FASTQ
        # -------------------------------------------------------------

        importFastqProtocol = self.newProtocol(
            ProtImportFastq,
            objLabel='Import Picard test FASTQ',
            inputFastq1=smallFastq,
            sampleName=sampleName,
            isPaired=False,
            runFastqc=False
        )

        self.launchProtocol(
            importFastqProtocol
        )

        self.assertTrue(
            hasattr(importFastqProtocol, 'outputFastq'),
            'ProtImportFastq did not produce outputFastq.'
        )

        # -------------------------------------------------------------
        # Import reduced reference genome
        # -------------------------------------------------------------

        genomesData = [
            {
                'scientificName': 'Mus musculus',
                'assembly': 'GRCm39',
                'release': '116',
                'fastaFile': referenceFasta,
                'gtfFile': referenceGtf
            }
        ]

        importGenomeProtocol = self.newProtocol(
            ProtImportGenomes,
            objLabel='Import Picard test genome',
            genomesData=json.dumps(genomesData)
        )

        self.launchProtocol(
            importGenomeProtocol
        )

        self.assertTrue(
            hasattr(importGenomeProtocol, 'referenceGenomes'),
            'ProtImportGenomes did not produce referenceGenomes.'
        )

        self.assertEqual(
            importGenomeProtocol.referenceGenomes.getSize(),
            1
        )

        # -------------------------------------------------------------
        # Create STAR alignment
        # -------------------------------------------------------------

        alignmentProtocol = self.newProtocol(
            ProtRNASeqAlignment,
            objLabel='Create Picard test alignment',
            genomeIndex=0,
            aligner=ProtRNASeqAlignment.ALIGN_STAR,
            rnaStrandness=0,
            keepIntermediateFiles=False,
            numberOfThreads=4,
            numberOfMpi=1
        )

        alignmentProtocol.inputFastq.set(
            importFastqProtocol.outputFastq
        )

        alignmentProtocol.inputGenomes.set(
            importGenomeProtocol.referenceGenomes
        )

        self.launchProtocol(
            alignmentProtocol
        )

        self.assertTrue(
            hasattr(alignmentProtocol, 'outputAlignment'),
            'ProtRNASeqAlignment did not produce outputAlignment.'
        )

        alignment = alignmentProtocol.outputAlignment

        self.assertTrue(
            os.path.isfile(alignment.getFileName()),
            'Input BAM for Picard was not created.'
        )

        self.assertGreater(
            os.path.getsize(alignment.getFileName()),
            0,
            'Input BAM for Picard is empty.'
        )

        return alignment

    # ---------------------------------------------------------------------
    # Output checks
    # ---------------------------------------------------------------------

    def _checkPicardOutput(
            self,
            protocol,
            expectedSampleName
    ):
        """Check the Picard output AlignmentFile and BAM index."""
        self.assertTrue(
            hasattr(protocol, 'outputAlignment'),
            'ProtPicard did not produce outputAlignment.'
        )

        output = protocol.outputAlignment

        bamFile = output.getFileName()

        self.assertTrue(
            os.path.isfile(bamFile),
            'Picard BAM was not created: {}'.format(
                bamFile
            )
        )

        self.assertGreater(
            os.path.getsize(bamFile),
            0,
            'Picard BAM is empty: {}'.format(
                bamFile
            )
        )

        baiFile = output.getIndexFile()

        self.assertTrue(
            os.path.isfile(baiFile),
            'Picard BAM index was not created: {}'.format(
                baiFile
            )
        )

        self.assertGreater(
            os.path.getsize(baiFile),
            0,
            'Picard BAM index is empty: {}'.format(
                baiFile
            )
        )

        self.assertTrue(
            output.isIndexed()
        )

        self.assertTrue(
            output.isSorted()
        )

        self.assertEqual(
            output.getSampleName(),
            expectedSampleName
        )

        return output

    # ---------------------------------------------------------------------
    # Tests
    # ---------------------------------------------------------------------

    def testPicardDefaultProcessing(self):
        """Test AddOrReplaceReadGroups followed by MarkDuplicates."""
        print(
            '\nPicard: AddOrReplaceReadGroups + MarkDuplicates'
        )

        inputAlignment = self._createInputAlignment()

        protocol = self.newProtocol(
            ProtPicard,
            objLabel='Picard default processing',
            addReadGroups=True,
            markDuplicates=True,
            reorderBam=False
        )

        protocol.inputAlignment.set(
            inputAlignment
        )

        self.launchProtocol(
            protocol
        )

        output = self._checkPicardOutput(
            protocol,
            inputAlignment.getSampleName()
        )

        metricsFile = (
            protocol._getMarkDuplicatesMetricsFile()
        )

        self.assertTrue(
            os.path.isfile(metricsFile),
            'Picard duplication metrics file was not created.'
        )

        self.assertGreater(
            os.path.getsize(metricsFile),
            0,
            'Picard duplication metrics file is empty.'
        )

        self.assertEqual(
            os.path.abspath(
                output.getReferenceFasta()
            ),
            os.path.abspath(
                inputAlignment.getReferenceFasta()
            )
        )

    def testPicardReorderBam(self):
        """Test Picard ReorderSam using the alignment reference FASTA."""
        print(
            '\nPicard: ReorderSam'
        )

        inputAlignment = self._createInputAlignment()

        protocol = self.newProtocol(
            ProtPicard,
            objLabel='Picard ReorderSam',
            addReadGroups=False,
            markDuplicates=False,
            reorderBam=True
        )

        protocol.inputAlignment.set(
            inputAlignment
        )

        self.launchProtocol(
            protocol
        )

        output = self._checkPicardOutput(
            protocol,
            inputAlignment.getSampleName()
        )

        self.assertEqual(
            os.path.abspath(
                output.getReferenceFasta()
            ),
            os.path.abspath(
                inputAlignment.getReferenceFasta()
            )
        )

        self.assertTrue(
            os.path.isfile(
                protocol._getReorderedBam()
            ),
            'Reordered BAM was not created.'
        )