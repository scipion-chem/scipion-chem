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

from pwchem.protocols import ProtImportFastq, ProtFastpFilter
from pwchem.utils import assertHandle


def assertOutputExists(test, protocol, outputFastq):
    assertHandle(
        test.assertIsNotNone,
        outputFastq,
        cwd=protocol.getWorkingDir()
    )


class TestFastpFilter(BaseTest):

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

    def _importFastq(self, sampleName, isPaired):
        kwargs = {
            'sampleName': sampleName,
            'isPaired': isPaired,
            'runFastqc': False
        }

        if isPaired:
            kwargs['inputFastq1'] = self.humanFastqR1
            kwargs['inputFastq2'] = self.humanFastqR2
        else:
            kwargs['inputFastq1'] = self.mouseFastq

        protImport = self.newProtocol(
            ProtImportFastq,
            **kwargs
        )

        self.launchProtocol(
            protImport
        )

        return protImport.outputFastq

    def _runFastp(self, inputFastq):
        protFastp = self.newProtocol(
            ProtFastpFilter,
            inputFastq=inputFastq,
            runFastqc=True,
            lengthRequired=15,
            qualifiedQualityPhred=15,
            unqualifiedPercentLimit=40,
            nBaseLimit=10,
            averageQual=0,
            disableAdapterTrimming=True
        )

        protFastp.numberOfThreads.set(1)

        self.launchProtocol(
            protFastp
        )

        return protFastp

    def _testFastpFilter(self, isPaired):
        sampleName = (
            'human_paired'
            if isPaired
            else 'mouse_single'
        )

        inputFastq = self._importFastq(
            sampleName=sampleName,
            isPaired=isPaired
        )

        protFastp = self._runFastp(
            inputFastq
        )

        outputFastq = getattr(
            protFastp,
            'outputFastq',
            None
        )

        assertOutputExists(
            self,
            protFastp,
            outputFastq
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getSampleName(),
            sampleName,
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertGreater,
            outputFastq.getNumReads(),
            0,
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertGreater,
            outputFastq.getReadLength(),
            0,
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getFormat(),
            'FASTQ',
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasQuality(),
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertFalse,
            outputFastq.isCompressed(),
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.isPaired(),
            isPaired,
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastpHtml(),
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastpJson(),
            cwd=protFastp.getWorkingDir()
        )

        if isPaired:
            assertHandle(
                self.assertTrue,
                outputFastq.hasFileName2(),
                cwd=protFastp.getWorkingDir()
            )

            assertHandle(
                self.assertTrue,
                outputFastq.hasFastqcHtmlR1(),
                cwd=protFastp.getWorkingDir()
            )

            assertHandle(
                self.assertTrue,
                outputFastq.hasFastqcHtmlR2(),
                cwd=protFastp.getWorkingDir()
            )

        else:
            assertHandle(
                self.assertTrue,
                outputFastq.hasFastqcHtml(),
                cwd=protFastp.getWorkingDir()
            )

    def testFastpFilterSingleEnd(self):
        print("\nFASTP filter: single-end")

        self._testFastpFilter(
            isPaired=False
        )

    def testFastpFilterPairedEnd(self):
        print("\nFASTP filter: paired-end")

        self._testFastpFilter(
            isPaired=True
        )