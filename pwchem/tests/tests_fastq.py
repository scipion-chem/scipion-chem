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

# Scipion em imports
from pyworkflow.tests import BaseTest, setupTestProject

# Scipion chem imports
from pwchem.protocols import ProtImportFastq, ProtFastpFilter
from pwchem.utils import assertHandle

FASTQ_R1_CONTENT = (
    '@READ_1\n'
    'ACGTACGTACGTACGTACGT\n'
    '+\n'
    'IIIIIIIIIIIIIIIIIIII\n'
    '@READ_2\n'
    'TTTTCCCCAAAAGGGGTTTT\n'
    '+\n'
    'IIIIIIIIIIIIIIIIIIII\n'
)

FASTQ_R2_CONTENT = (
    '@READ_1\n'
    'TGCATGCATGCATGCATGCA\n'
    '+\n'
    'IIIIIIIIIIIIIIIIIIII\n'
    '@READ_2\n'
    'AAAACCCCGGGGTTTTAAAA\n'
    '+\n'
    'IIIIIIIIIIIIIIIIIIII\n'
)

def write_fastq_files(cls, prefix):
    cls.fastqR1 = cls.proj.getTmpPath(f'{prefix}_R1.fastq')
    cls.fastqR2 = cls.proj.getTmpPath(f'{prefix}_R2.fastq')

    with open(cls.fastqR1, 'w') as f:
        f.write(FASTQ_R1_CONTENT)

    with open(cls.fastqR2, 'w') as f:
        f.write(FASTQ_R2_CONTENT)


def assert_output_exists(test, protocol, outputFastq):
    assertHandle(
        test.assertIsNotNone,
        outputFastq,
        cwd=protocol.getWorkingDir()
    )


def assert_fastq_stats(test, protocol, outputFastq, sampleName,
                      numReads=2, readLength=20):
    assertHandle(
        test.assertEqual,
        outputFastq.getSampleName(),
        sampleName,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        outputFastq.getNumReads(),
        numReads,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        outputFastq.getReadLength(),
        readLength,
        cwd=protocol.getWorkingDir()
    )


class TestImportFastq(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        write_fastq_files(cls, 'sample')

    def testImportSingleFastqWithFastqc(self):
        print("\nImport FASTQ: single-end with FastQC")

        prot = self.newProtocol(
            ProtImportFastq,
            sampleName='test_single',
            isPaired=False,
            inputFastq1=self.fastqR1,
            runFastqc=True
        )

        self.launchProtocol(prot)

        outputFastq = getattr(prot, 'outputFastq', None)

        assert_output_exists(self, prot, outputFastq)
        assert_fastq_stats(self, prot, outputFastq, 'test_single')

        assertHandle(
            self.assertFalse,
            outputFastq.isPaired(),
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
            sampleName='test_paired',
            isPaired=True,
            inputFastq1=self.fastqR1,
            inputFastq2=self.fastqR2,
            runFastqc=True
        )

        self.launchProtocol(prot)

        outputFastq = getattr(prot, 'outputFastq', None)

        assert_output_exists(self, prot, outputFastq)
        assert_fastq_stats(self, prot, outputFastq, 'test_paired')

        assertHandle(
            self.assertTrue,
            outputFastq.isPaired(),
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            outputFastq.getFileName2(),
            self.fastqR2,
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


class TestFastpFilter(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        write_fastq_files(cls, 'fastp_sample')

    def _importFastq(self, sampleName, isPaired):
        kwargs = {
            'sampleName': sampleName,
            'isPaired': isPaired,
            'inputFastq1': self.fastqR1,
            'runFastqc': False
        }

        if isPaired:
            kwargs['inputFastq2'] = self.fastqR2

        protImport = self.newProtocol(ProtImportFastq, **kwargs)
        self.launchProtocol(protImport)

        return protImport

    def _runFastp(self, inputFastq):
        protFastp = self.newProtocol(
            ProtFastpFilter,
            inputFastq=inputFastq,
            runFastqc=False,
            lengthRequired=1,
            threads=1,
            qualifiedQualityPhred=15,
            unqualifiedPercentLimit=40,
            nBaseLimit=10,
            averageQual=0,
            disableAdapterTrimming=True
        )

        self.launchProtocol(protFastp)

        return protFastp

    def _assertFastpReports(self, protocol, outputFastq):
        assertHandle(
            self.assertTrue,
            outputFastq.hasFastpHtml(),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFastpJson(),
            cwd=protocol.getWorkingDir()
        )

    def testFastpFilterSingleEnd(self):
        print("\nFASTP filter: single-end")

        protImport = self._importFastq(
            sampleName='fastp_single',
            isPaired=False
        )

        protFastp = self._runFastp(protImport.outputFastq)
        outputFastq = getattr(protFastp, 'outputFastq', None)

        assert_output_exists(self, protFastp, outputFastq)
        assert_fastq_stats(self, protFastp, outputFastq, 'fastp_single')

        assertHandle(
            self.assertFalse,
            outputFastq.isPaired(),
            cwd=protFastp.getWorkingDir()
        )

        self._assertFastpReports(protFastp, outputFastq)

    def testFastpFilterPairedEnd(self):
        print("\nFASTP filter: paired-end")

        protImport = self._importFastq(
            sampleName='fastp_paired',
            isPaired=True
        )

        protFastp = self._runFastp(protImport.outputFastq)
        outputFastq = getattr(protFastp, 'outputFastq', None)

        assert_output_exists(self, protFastp, outputFastq)
        assert_fastq_stats(self, protFastp, outputFastq, 'fastp_paired')

        assertHandle(
            self.assertTrue,
            outputFastq.isPaired(),
            cwd=protFastp.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            outputFastq.hasFileName2(),
            cwd=protFastp.getWorkingDir()
        )

        self._assertFastpReports(protFastp, outputFastq)