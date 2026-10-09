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
# ***************************************************************************

import gzip
import json
import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import (
    ProtImportFastq,
    ProtImportGenomes,
    ProtImportVCF,
    ProtRNASeqAlignment,
    ProtPicard,
    ProtGATK
)


N_TEST_READS = 1000


class TestGATK(BaseTest):
    """Test GATK processing of RNA-seq BAM alignments."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')

        cls.fastqFile = cls.dataset.getFile(
            'mouseFastq'
        )

        cls.genomeFile = cls.dataset.getFile(
            'mouseGenomeFile'
        )

        cls.gtfFile = cls.dataset.getFile(
            'mouseGtfFile'
        )

        cls.vcfFile = cls.dataset.getFile(
            'mouseVcfFile'
        )

        cls.vcfIndex = cls.dataset.getFile(
            'mouseVcfIndex'
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
            'mouse_gatk_test.fastq'
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

    def _createInputAlignment(self):
        """Create a Picard-processed STAR alignment for GATK."""
        smallFastq = self._createSmallFastq()

        sampleName = 'mouse_gatk_test'

        # -------------------------------------------------------------
        # Import FASTQ
        # -------------------------------------------------------------

        importFastqProtocol = self.newProtocol(
            ProtImportFastq,
            objLabel='Import GATK test FASTQ',
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
        # Import reference genome
        # -------------------------------------------------------------

        genomesData = json.dumps([
            {
                'scientificName': 'Mus musculus',
                'assembly': 'GRCm39',
                'release': '116',
                'fastaFile': self.genomeFile,
                'gtfFile': self.gtfFile
            }
        ])

        importGenomeProtocol = self.newProtocol(
            ProtImportGenomes,
            objLabel='Import GATK test genome',
            genomesData=genomesData
        )

        self.launchProtocol(
            importGenomeProtocol
        )

        self.assertTrue(
            hasattr(importGenomeProtocol, 'referenceGenomes'),
            'ProtImportGenomes did not produce referenceGenomes.'
        )

        # -------------------------------------------------------------
        # STAR alignment
        # -------------------------------------------------------------

        alignmentProtocol = self.newProtocol(
            ProtRNASeqAlignment,
            objLabel='Create GATK test alignment',
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

        # -------------------------------------------------------------
        # Picard preprocessing
        # -------------------------------------------------------------

        picardProtocol = self.newProtocol(
            ProtPicard,
            objLabel='Prepare BAM for GATK',
            addReadGroups=True,
            markDuplicates=True,
            reorderBam=False
        )

        picardProtocol.inputAlignment.set(
            alignmentProtocol.outputAlignment
        )

        self.launchProtocol(
            picardProtocol
        )

        self.assertTrue(
            hasattr(picardProtocol, 'outputAlignment'),
            'ProtPicard did not produce outputAlignment.'
        )

        bamFile = (
            picardProtocol.outputAlignment.getFileName()
        )

        self.assertTrue(
            os.path.isfile(bamFile),
            'Picard BAM was not created.'
        )

        self.assertGreater(
            os.path.getsize(bamFile),
            0,
            'Picard BAM is empty.'
        )

        return picardProtocol.outputAlignment

    def _importKnownSites(self):
        """Import the mouse known-sites VCF used by BaseRecalibrator."""

        vcfsData = json.dumps([
            {
                'scientificName': 'Mus musculus',
                'assembly': 'GRCm39',
                'vcfFile': self.vcfFile,
                'indexFile': self.vcfIndex
            }
        ])

        importVCFProtocol = self.newProtocol(
            ProtImportVCF,
            objLabel='Import GATK known-sites VCF',
            vcfsData=vcfsData
        )

        self.launchProtocol(
            importVCFProtocol
        )

        self.assertTrue(
            hasattr(importVCFProtocol, 'outputVCFs'),
            'ProtImportVCF did not produce outputVCFs.'
        )

        self.assertEqual(
            importVCFProtocol.outputVCFs.getSize(),
            1,
            'ProtImportVCF should produce one VCF.'
        )

        return importVCFProtocol.outputVCFs

    # ---------------------------------------------------------------------
    # Output checks
    # ---------------------------------------------------------------------

    def _checkGATKOutput(self, protocol):
        """Check the GATK output BAM, index and recalibration table."""

        self.assertTrue(
            hasattr(protocol, 'outputAlignment'),
            'ProtGATK did not produce outputAlignment.'
        )

        output = protocol.outputAlignment
        bamFile = output.getFileName()

        self.assertTrue(
            os.path.isfile(bamFile),
            'GATK BAM was not created: {}'.format(
                bamFile
            )
        )

        self.assertGreater(
            os.path.getsize(bamFile),
            0,
            'GATK BAM is empty: {}'.format(
                bamFile
            )
        )

        self.assertTrue(
            os.path.isfile(bamFile + '.bai'),
            'GATK BAM index was not created.'
        )

        recalibrationTable = (
            protocol._getRecalibrationTable()
        )

        self.assertTrue(
            os.path.isfile(recalibrationTable),
            'GATK recalibration table was not created.'
        )

        self.assertGreater(
            os.path.getsize(recalibrationTable),
            0,
            'GATK recalibration table is empty.'
        )

    # ---------------------------------------------------------------------
    # Tests
    # ---------------------------------------------------------------------

    def testGATKProcessing(self):
        """
        Test the complete RNA-seq GATK preprocessing workflow:

        SplitNCigarReads -> BaseRecalibrator -> ApplyBQSR
        """

        inputAlignment = (
            self._createInputAlignment()
        )

        knownSites = (
            self._importKnownSites()
        )

        gatkProtocol = self.newProtocol(
            ProtGATK,
            objLabel='GATK RNA-seq processing',
            splitNCigarReads=True,
            baseRecalibrator=True,
            vcfIndex=0,
            applyBQSR=True,
            numberOfThreads=4,
            numberOfMpi=1
        )

        gatkProtocol.inputAlignment.set(
            inputAlignment
        )

        gatkProtocol.knownSites.set(
            knownSites
        )

        self.launchProtocol(
            gatkProtocol
        )

        self._checkGATKOutput(
            gatkProtocol
        )