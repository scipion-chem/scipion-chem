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
    ProtRNASeqAlignment
)


N_TEST_READS = 1000


class TestRNASeqAlignment(BaseTest):
    """Test RNA-seq alignment with STAR and HISAT2."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')
        # Mouse
        cls.mouseFastq = cls.dataset.getFile('mouseFastq')
        cls.mouseGenomeFile = cls.dataset.getFile('mouseGenomeFile')
        cls.mouseGtfFile = cls.dataset.getFile('mouseGtfFile')

        # Human
        cls.humanFastq = cls.dataset.getFile('humanFastqR1')
        cls.humanGenomeFile = cls.dataset.getFile('humanGenomeFile')
        cls.humanGtfFile = cls.dataset.getFile('humanGtfFile')

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
    # Genome import
    # ---------------------------------------------------------------------

    def _importGenome(
            self,
            fastaFile,
            gtfFile,
            scientificName,
            assembly,
            release
    ):
        """Import a reduced reference genome as a SetOfGenomes."""

        genomesData = [
            {
                'scientificName': scientificName,
                'assembly': assembly,
                'release': release,
                'fastaFile': fastaFile,
                'gtfFile': gtfFile
            }
        ]

        protocol = self.newProtocol(
            ProtImportGenomes,
            objLabel='Import {} {}'.format(
                scientificName,
                assembly
            ),
            genomesData=json.dumps(genomesData)
        )

        self.launchProtocol(protocol)

        self.assertTrue(
            hasattr(protocol, 'referenceGenomes'),
            'ProtImportGenomes did not produce referenceGenomes.'
        )

        genomeSet = protocol.referenceGenomes

        self.assertEqual(
            genomeSet.getSize(),
            1
        )

        return genomeSet

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
        Test STAR and HISAT2 using the same test dataset.

        The reference genome already contains one chromosome. The FASTQ
        is reduced to N_TEST_READS reads and reused for STAR and HISAT2.
        """
        # -------------------------------------------------------------
        # Import reference once
        # -------------------------------------------------------------

        genomeSet = self._importGenome(
            fastaFile,
            gtfFile,
            scientificName,
            assembly,
            release
        )

        genome = genomeSet.getFirstItem()

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
                genomeIndex=0,
                aligner=aligner,
                rnaStrandness=0,
                keepIntermediateFiles=False,
                numberOfThreads=4,
                numberOfMpi=1
            )

            protocol.inputFastq.set(
                fastqObj
            )

            protocol.inputGenomes.set(
                genomeSet
            )

            self.launchProtocol(protocol)

            self._checkAlignment(
                protocol,
                alignerName,
                sampleName,
                genome.getFastaFile(),
                genome.getGtfFile()
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