# **************************************************************************
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
# **************************************************************************

import os

from pyworkflow.tests import BaseTest, setupTestProject

from pwchem.protocols import ProtDownloadVCF
from pwchem.utils import assertHandle
from pwchem.utils.utilsRNA import assertOutputExists


def assertVcfFile(test, protocol, vcf):
    """Check that the downloaded VCF file exists and is not empty."""
    fileName = vcf.getFileName()

    assertHandle(
        test.assertTrue,
        bool(fileName),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        os.path.exists(fileName),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertGreater,
        os.path.getsize(fileName),
        0,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        fileName.endswith(('.vcf', '.vcf.gz')),
        cwd=protocol.getWorkingDir()
    )


class TestDownloadVCF(BaseTest):
    """Test known-variant VCF downloads from Ensembl, NCBI and EVA."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def testEnsemblVCF(self):
        """Download a known-variant VCF from Ensembl."""
        print("\nKnown variants VCF: Ensembl")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_ENSEMBL,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Saccharomyces cerevisiae',
            assemblies='Latest',
            releases='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        assertVcfFile(
            self,
            prot,
            vcf
        )

    def testNcbiHumanVCF(self):
        """Download the human dbSNP VCF from NCBI."""
        print("\nKnown variants VCF: NCBI dbSNP")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Homo sapiens',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        assertVcfFile(
            self,
            prot,
            vcf
        )

    def testNcbiEvaFallback(self):
        """
        Check the automatic NCBI -> EVA fallback for a non-human species.

        When NCBI is selected for a non-human species, dbSNP does not
        provide the current VCF and the protocol delegates the download
        to EVA.
        """
        print("\nKnown variants VCF: NCBI -> EVA fallback")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Oryza sativa',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        assertVcfFile(
            self,
            prot,
            vcf
        )

    def testMultipleEvaFallbackVCFs(self):
        """
        Download VCFs for multiple non-human species through EVA.
        """
        print("\nKnown variants VCF: multiple NCBI -> EVA fallbacks")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Bos taurus;Oryza sativa',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            2,
            cwd=prot.getWorkingDir()
        )

        for vcf in output:
            assertVcfFile(
                self,
                prot,
                vcf
            )