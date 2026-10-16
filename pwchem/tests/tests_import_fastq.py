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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import ProtImportFastq
from pwchem.utils import assertHandle


def assertOutputExists(test, protocol, outputFastq):
    assertHandle(
        test.assertIsNotNone,
        outputFastq,
        cwd=protocol.getWorkingDir()
    )


def assertFastqProperties(test, protocol, outputFastq, sampleName):
    assertHandle(
        test.assertEqual,
        outputFastq.getSampleName(),
        sampleName,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreater,
        outputFastq.getNumReads(),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreater,
        outputFastq.getReadLength(),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        outputFastq.getFormat(),
        'FASTQ',
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        outputFastq.hasQuality(),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        outputFastq.isCompressed(),
        cwd=protocol.getWorkingDir()
    )


class TestImportFastq(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet('genomics')

        cls.mouseFastq = cls.dataset.getFile(
            'mouseFastq'
        )

        cls.humanFastqR1 = cls.dataset.getFile(
            'humanFastqR1'
        )

        cls.humanFastqR2 = cls.dataset.getFile(
            'humanFastqR2'
        )

    def testImportSingleFastqWithFastqc(self):
        print("\nImport FASTQ: single-end with FastQC")

        prot = self.newProtocol(
            ProtImportFastq,
            sampleName='mouse_single',
            isPaired=False,
            inputFastq1=self.mouseFastq,
            runFastqc=True
        )

        self.launchProtocol(prot)

        outputFastq = getattr(
            prot,
            'outputFastq',
            None
        )

        assertOutputExists(
            self,
            prot,
            outputFastq
        )

        assertFastqProperties(
            self,
            prot,
            outputFastq,
            'mouse_single'
        )

        assertHandle(
            self.assertFalse,
            outputFastq.isPaired(),
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getFileName(),
            self.mouseFastq,
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastqcHtml(),
            cwd=prot.getWorkingDir()
        )

    def testImportPairedFastqWithFastqc(self):
        print("\nImport FASTQ: paired-end with FastQC")

        prot = self.newProtocol(
            ProtImportFastq,
            sampleName='human_paired',
            isPaired=True,
            inputFastq1=self.humanFastqR1,
            inputFastq2=self.humanFastqR2,
            runFastqc=True
        )

        self.launchProtocol(prot)

        outputFastq = getattr(
            prot,
            'outputFastq',
            None
        )

        assertOutputExists(
            self,
            prot,
            outputFastq
        )

        assertFastqProperties(
            self,
            prot,
            outputFastq,
            'human_paired'
        )

        assertHandle(
            self.assertTrue,
            outputFastq.isPaired(),
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getFileName(),
            self.humanFastqR1,
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getFileName2(),
            self.humanFastqR2,
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastqcHtmlR1(),
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastqcHtmlR2(),
            cwd=prot.getWorkingDir()
        )