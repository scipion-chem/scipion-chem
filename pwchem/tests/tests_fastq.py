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
from pwchem.protocols import ProtImportFastq
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

class TestImportFastq(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.fastqR1 = cls.proj.getTmpPath('sample_R1.fastq')
        cls.fastqR2 = cls.proj.getTmpPath('sample_R2.fastq')

        with open(cls.fastqR1, 'w') as f:
            f.write(FASTQ_R1_CONTENT)

        with open(cls.fastqR2, 'w') as f:
            f.write(FASTQ_R2_CONTENT)

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

        assertHandle(
            self.assertIsNotNone,
            outputFastq,
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

        assertHandle(
            self.assertIsNotNone,
            outputFastq,
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