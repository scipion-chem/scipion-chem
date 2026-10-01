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
import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import ProtImportFastq, ProtRNASeqAlignment
N_TEST_READS = 1000


class TestRNASeqAlignment(BaseTest):
    """Test RNA-seq alignment with STAR and HISAT2."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')

        # Mouse
        cls.mouseFastq = cls.dataset.getFile(
            'mouse/SRR1552445.fastq.gz'
        )
        cls.mouseGenomeFile = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa'
        )
        cls.mouseGtfFile = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.116.gtf'
        )

        # Human
        cls.humanFastq = cls.dataset.getFile(
            'human/SRR390728_1.fastq.gz'
        )
        cls.humanGenomeFile = cls.dataset.getFile(
            'human/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
        )
        cls.humanGtfFile = cls.dataset.getFile(
            'human/Homo_sapiens.GRCh38.116.gtf'
        )

    # ---------------------------------------------------------------------
    # File helpers
    # ---------------------------------------------------------------------

    @staticmethod
    def _openTextFile(fileName):
        """Open plain-text or gzip-compressed files in text mode."""
        if fileName.endswith('.gz'):
            return gzip.open(fileName, 'rt')

        return open(fileName, 'r')

    def _createSmallFastq(self, fastqFile, species):
        """
        Create a reduced FASTQ containing the first N_TEST_READS reads.

        FASTQ records contain four lines each, so the output contains at
        most N_TEST_READS * 4 lines.
        """
        outputFastq = self.getOutputPath(
            '{}_test.fastq'.format(species)
        )

        maxLines = N_TEST_READS * 4

        with self._openTextFile(fastqFile) as inputFile, \
                open(outputFastq, 'w') as outputFile:

            for lineNumber, line in enumerate(inputFile):
                if lineNumber >= maxLines:
                    break

                outputFile.write(line)

        if not os.path.isfile(outputFastq):
            raise RuntimeError(
                'Reduced FASTQ was not created: {}'.format(
                    outputFastq
                )
            )

        if os.path.getsize(outputFastq) == 0:
            raise RuntimeError(
                'Reduced FASTQ is empty: {}'.format(
                    outputFastq
                )
            )

        return outputFastq

    def _createChromosomeReference(
            self,
            fastaFile,
            gtfFile,
            species
    ):
        """
        Create a reduced reference containing one chromosome.

        The first sequence found in the FASTA file is selected. The
        corresponding sequence is written to a new FASTA file and the
        GTF annotation is filtered to retain only entries belonging to
        that sequence.
        """
        outputFasta = self.getOutputPath(
            '{}_test.fa'.format(species)
        )
        outputGtf = self.getOutputPath(
            '{}_test.gtf'.format(species)
        )

        chromosome = None

        # -------------------------------------------------------------
        # Reduce FASTA
        # -------------------------------------------------------------

        with self._openTextFile(fastaFile) as inputFile, \
                open(outputFasta, 'w') as outputFile:

            for line in inputFile:

                if line.startswith('>'):
                    currentChromosome = line[1:].split()[0]

                    if chromosome is None:
                        chromosome = currentChromosome

                    elif currentChromosome != chromosome:
                        break

                outputFile.write(line)

        if chromosome is None:
            raise RuntimeError(
                'No sequence was found in FASTA file: {}'.format(
                    fastaFile
                )
            )

        if not os.path.isfile(outputFasta):
            raise RuntimeError(
                'Reduced FASTA was not created: {}'.format(
                    outputFasta
                )
            )

        if os.path.getsize(outputFasta) == 0:
            raise RuntimeError(
                'Reduced FASTA is empty: {}'.format(
                    outputFasta
                )
            )

        # -------------------------------------------------------------
        # Reduce GTF
        # -------------------------------------------------------------

        featureCount = 0

        with self._openTextFile(gtfFile) as inputFile, \
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

        if not os.path.isfile(outputGtf):
            raise RuntimeError(
                'Reduced GTF was not created: {}'.format(
                    outputGtf
                )
            )

        if featureCount == 0:
            raise RuntimeError(
                'No GTF annotations were found for sequence "{}".'
                .format(chromosome)
            )

        print(
            'Using {} sequence "{}" as reduced reference.'
            .format(species, chromosome)
        )

        return outputFasta, outputGtf

    # ---------------------------------------------------------------------
    # FASTQ import
    # ---------------------------------------------------------------------

    def _importFastq(self, fastqFile, sampleName):
        """Import a reduced single-end FASTQ file."""
        protocol = self.newProtocol(
            ProtImportFastq,
            objLabel='Import {}'.format(sampleName),
            inputFastq1=fastqFile,
            sampleName=sampleName,
            isPaired=False,
            runFastqc=False
        )

        self.launchProtocol(protocol)

        self.assertTrue(
            hasattr(protocol, 'outputFastq'),
            'ProtImportFastq did not produce outputFastq.'
        )

        fastqObj = protocol.outputFastq

        self.assertTrue(
            os.path.isfile(fastqObj.getFileName()),
            'Imported FASTQ file does not exist: {}'.format(
                fastqObj.getFileName()
            )
        )

        self.assertEqual(
            fastqObj.getNumReads(),
            N_TEST_READS
        )

        return fastqObj
    # ---------------------------------------------------------------------
    # Alignment checks
    # ---------------------------------------------------------------------

    def _checkAlignment(
            self,
            protocol,
            alignerName,
            sampleName,
            referenceFasta,
            referenceGtf
    ):
        """Check the alignment output and its basic properties."""

        self.assertTrue(
            hasattr(protocol, 'outputAlignment'),
            '{} did not produce outputAlignment.'.format(
                alignerName
            )
        )

        output = protocol.outputAlignment

        bamFile = output.getFileName()

        self.assertTrue(
            os.path.isfile(bamFile),
            'BAM file was not created: {}'.format(bamFile)
        )

        self.assertGreater(
            os.path.getsize(bamFile),
            0,
            'BAM file is empty: {}'.format(bamFile)
        )

        baiFile = output.getIndexFile()

        self.assertTrue(
            os.path.isfile(baiFile),
            'BAI file was not created: {}'.format(baiFile)
        )

        self.assertGreater(
            os.path.getsize(baiFile),
            0,
            'BAI file is empty: {}'.format(baiFile)
        )

        self.assertEqual(
            output.getFormat(),
            'BAM'
        )

        self.assertEqual(
            output.getSampleName(),
            sampleName
        )

        self.assertEqual(
            output.getAligner(),
            alignerName
        )

        self.assertTrue(
            output.isSorted()
        )

        self.assertTrue(
            output.isIndexed()
        )

        self.assertEqual(
            os.path.abspath(output.getReferenceFasta()),
            os.path.abspath(referenceFasta)
        )

        if referenceGtf:
            self.assertEqual(
                os.path.abspath(output.getReferenceGtf()),
                os.path.abspath(referenceGtf)
            )

        self._checkStatistics(output)

    def _checkStatistics(self, output):
        """Check basic alignment statistics."""
        totalReads = output.getTotalReads()
        mappedReads = output.getMappedReads()
        unmappedReads = output.getUnmappedReads()
        mappingRate = output.getMappingRate()

        self.assertGreater(
            totalReads,
            0
        )

        self.assertGreaterEqual(
            mappedReads,
            0
        )

        self.assertGreaterEqual(
            unmappedReads,
            0
        )

        self.assertGreaterEqual(
            mappingRate,
            0.0
        )

        self.assertLessEqual(
            mappingRate,
            100.0
        )

        self.assertEqual(
            mappedReads + unmappedReads,
            totalReads
        )

    # ---------------------------------------------------------------------
    # Species test
    # ---------------------------------------------------------------------

    def _testSpecies(
            self,
            species,
            fastqFile,
            fastaFile,
            gtfFile,
            scientificName,
            assembly,
            release
    ):
        """
        Test STAR and HISAT2 using the same reduced dataset.

        The reference is reduced to one chromosome and the FASTQ to
        N_TEST_READS reads. Both reduced inputs are generated only once
        and reused for STAR and HISAT2.
        """

        # -------------------------------------------------------------
        # Reduce reference once
        # -------------------------------------------------------------

        smallFasta, smallGtf = self._createChromosomeReference(
            fastaFile,
            gtfFile,
            species
        )

        # -------------------------------------------------------------
        # Reduce FASTQ once
        # -------------------------------------------------------------

        smallFastq = self._createSmallFastq(
            fastqFile,
            species
        )

        # -------------------------------------------------------------
        # Import FASTQ once
        # -------------------------------------------------------------

        sampleName = '{}_single'.format(species)

        fastqObj = self._importFastq(
            smallFastq,
            sampleName
        )

        # -------------------------------------------------------------
        # Run both aligners
        # -------------------------------------------------------------

        aligners = (
            (
                ProtRNASeqAlignment.ALIGN_STAR,
                'STAR'
            ),
            (
                ProtRNASeqAlignment.ALIGN_HISAT2,
                'HISAT2'
            )
        )

        for aligner, alignerName in aligners:

            print(
                '\nRNA-seq alignment: {} {} single-end'
                .format(
                    alignerName,
                    species
                )
            )

            protocol = self.newProtocol(
                ProtRNASeqAlignment,
                objLabel='{} {} alignment'.format(
                    alignerName,
                    species
                ),
                referenceSource=(
                    ProtRNASeqAlignment.REFERENCE_FROM_FILES
                ),
                manualFasta=smallFasta,
                manualGtf=smallGtf,
                manualReferenceName='{} {} {}'.format(
                    scientificName,
                    assembly,
                    release
                ),
                manualReferenceSource='Test data',
                aligner=aligner,
                rnaStrandness=0,
                keepIntermediateFiles=False,
                numberOfThreads=4,
                numberOfMpi=1
            )

            protocol.inputFastq.set(fastqObj)

            self.launchProtocol(protocol)

            self._checkAlignment(
                protocol,
                alignerName,
                sampleName,
                smallFasta,
                smallGtf
            )

    # ---------------------------------------------------------------------
    # Test
    # ---------------------------------------------------------------------

    def testRNASeqAlignment(self):
        """Test STAR and HISAT2 with mouse and human RNA-seq data."""

        self._testSpecies(
            species='mouse',
            fastqFile=self.mouseFastq,
            fastaFile=self.mouseGenomeFile,
            gtfFile=self.mouseGtfFile,
            scientificName='Mus musculus',
            assembly='GRCm39',
            release='test'
        )

        self._testSpecies(
            species='human',
            fastqFile=self.humanFastq,
            fastaFile=self.humanGenomeFile,
            gtfFile=self.humanGtfFile,
            scientificName='Homo sapiens',
            assembly='GRCh38',
            release='test'
        )