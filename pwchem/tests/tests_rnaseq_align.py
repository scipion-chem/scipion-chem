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

import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwchem.protocols import ProtImportFastq, ProtRNASeqAlignment
from pwchem.utils import assertHandle


def assertOutputExists(test, protocol, output):
    assertHandle(
        test.assertIsNotNone,
        output,
        cwd=protocol.getWorkingDir()
    )


def assertAlignmentFiles(test, protocol, output):
    """Check that BAM and BAI files exist and are not empty."""

    assertHandle(
        test.assertTrue,
        os.path.exists(output.getFileName()),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreater,
        os.path.getsize(output.getFileName()),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        output.hasIndexFile(),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        os.path.exists(output.getIndexFile()),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreater,
        os.path.getsize(output.getIndexFile()),
        0,
        cwd=protocol.getWorkingDir()
    )


def assertAlignmentProperties(
        test,
        protocol,
        output,
        aligner,
        sampleName,
        fastaFile,
        gtfFile):
    """Check common AlignmentFile properties."""

    assertHandle(
        test.assertEqual,
        output.getFormat(),
        'BAM',
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        output.getSampleName(),
        sampleName,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        output.getAligner(),
        aligner,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        output.isSorted(),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        output.isIndexed(),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        output.getReferenceFasta(),
        os.path.abspath(fastaFile),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        output.hasReferenceGtf(),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        output.getReferenceGtf(),
        os.path.abspath(gtfFile),
        cwd=protocol.getWorkingDir()
    )


def assertAlignmentStatistics(test, protocol, output):
    """Check that alignment statistics are valid and internally consistent."""

    assertHandle(
        test.assertGreater,
        output.getTotalReads(),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreaterEqual,
        output.getMappedReads(),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreaterEqual,
        output.getUnmappedReads(),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreaterEqual,
        output.getMappingRate(),
        0.0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertLessEqual,
        output.getMappingRate(),
        100.0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        output.getMappedReads() + output.getUnmappedReads(),
        output.getTotalReads(),
        cwd=protocol.getWorkingDir()
    )


class TestRNASeqAlignment(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')

        cls.mouseFastq = cls.dataset.getFile(
            'mouseFastq'
        )

        cls.mouseGenome = cls.dataset.getFile(
            'mouseGenomeFile'
        )

        cls.mouseGtf = cls.dataset.getFile(
            'mouseGtfFile'
        )

    def _importMouseFastq(self):
        """Import the mouse single-end FASTQ as a FastqFile object."""

        protImport = self.newProtocol(
            ProtImportFastq,
            sampleName='mouse_single',
            isPaired=False,
            inputFastq1=self.mouseFastq,
            runFastqc=False
        )

        self.launchProtocol(protImport)

        outputFastq = getattr(
            protImport,
            'outputFastq',
            None
        )

        assertOutputExists(
            self,
            protImport,
            outputFastq
        )

        return outputFastq

    def _runAlignment(self, inputFastq, aligner):
        """Run RNA-seq alignment using the selected aligner."""

        protAlignment = self.newProtocol(
            ProtRNASeqAlignment,
            inputFastq=inputFastq,
            referenceSource=(
                ProtRNASeqAlignment.REFERENCE_FROM_FILES
            ),
            manualFasta=self.mouseGenome,
            manualGtf=self.mouseGtf,
            manualReferenceName='Mus musculus GRCm39',
            manualReferenceSource='Ensembl',
            aligner=aligner,
            rnaStrandness=0,
            keepIntermediateFiles=False,
            numberOfThreads=4,
            numberOfMpi=1
        )

        self.launchProtocol(protAlignment)

        return protAlignment

    def _assertAlignment(
            self,
            protocol,
            expectedAligner):

        output = getattr(
            protocol,
            'outputAlignment',
            None
        )

        assertOutputExists(
            self,
            protocol,
            output
        )

        assertAlignmentFiles(
            self,
            protocol,
            output
        )

        assertAlignmentProperties(
            self,
            protocol,
            output,
            expectedAligner,
            'mouse_single',
            self.mouseGenome,
            self.mouseGtf
        )

        assertAlignmentStatistics(
            self,
            protocol,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getReferenceName(),
            'Mus musculus GRCm39',
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            output.getReferenceSource(),
            'Ensembl',
            cwd=protocol.getWorkingDir()
        )

    def testStarAlignmentSingleEnd(self):
        print("\nRNA-seq alignment: STAR single-end")

        inputFastq = self._importMouseFastq()

        protAlignment = self._runAlignment(
            inputFastq,
            ProtRNASeqAlignment.ALIGN_STAR
        )

        self._assertAlignment(
            protAlignment,
            'STAR'
        )

    def testHisat2AlignmentSingleEnd(self):
        print("\nRNA-seq alignment: HISAT2 single-end")

        inputFastq = self._importMouseFastq()

        protAlignment = self._runAlignment(
            inputFastq,
            ProtRNASeqAlignment.ALIGN_HISAT2
        )

        self._assertAlignment(
            protAlignment,
            'HISAT2'
        )