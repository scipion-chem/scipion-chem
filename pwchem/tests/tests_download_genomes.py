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

from pyworkflow.tests import BaseTest, setupTestProject

from pwchem.protocols import ProtDownloadGenomes
from pwchem.utils import assertHandle
from pwchem.utils.utilsRNA import (
    assertOutputExists,
    assertGenomeFiles
)


def assertGenomeMetadata(test, protocol, genome, source):
    """Check the basic metadata stored in a downloaded Genome object."""
    assertHandle(
        test.assertTrue,
        bool(genome.getScientificName()),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        bool(genome.getAssembly()),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertTrue,
        bool(genome.getRelease()),
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        genome.getSource(),
        source,
        cwd=protocol.getWorkingDir()
    )


def assertResolvedGenome(test, protocol, genome):
    """Check that Latest values were resolved."""
    assertHandle(
        test.assertNotEqual,
        genome.getAssembly().lower(),
        'latest',
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertNotEqual,
        genome.getRelease().lower(),
        'latest',
        cwd=protocol.getWorkingDir()
    )


def assertNcbiAccession(test, protocol, genome):
    """Check that the NCBI release is a versioned assembly accession."""
    assertHandle(
        test.assertRegex,
        genome.getRelease(),
        r'^(GCF|GCA)_\d+\.\d+$',
        cwd=protocol.getWorkingDir()
    )


class TestDownloadGenomes(BaseTest):
    """Test reference genome downloads from Ensembl and NCBI."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def testEnsemblCommonGenome(self):
        """Download one common genome from Ensembl."""
        print("\nReference genomes: Ensembl common genome")

        prot = self.newProtocol(
            ProtDownloadGenomes,
            source=ProtDownloadGenomes.SOURCE_ENSEMBL,
            genomeSelection=ProtDownloadGenomes.GENOME_COMMON,
            commonSpecies='Saccharomyces cerevisiae',
            assemblies='Latest',
            releases='Latest',
            downloadAnnotation=True,
            overwrite=False
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'referenceGenomes',
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

        genome = output.getFirstItem()

        assertHandle(
            self.assertEqual,
            genome.getScientificName(),
            'Saccharomyces cerevisiae',
            cwd=prot.getWorkingDir()
        )

        assertGenomeMetadata(
            self,
            prot,
            genome,
            'Ensembl'
        )

        assertGenomeFiles(
            self,
            prot,
            genome,
            expectGtf=True
        )

        assertResolvedGenome(
            self,
            prot,
            genome
        )

    def testNcbiCommonGenome(self):
        """Download one common genome from NCBI."""
        print("\nReference genomes: NCBI common genome")

        prot = self.newProtocol(
            ProtDownloadGenomes,
            source=ProtDownloadGenomes.SOURCE_NCBI,
            genomeSelection=ProtDownloadGenomes.GENOME_COMMON,
            commonSpecies='Homo sapiens',
            ncbiAssemblies='Latest',
            downloadAnnotation=True,
            overwrite=False
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'referenceGenomes',
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

        genome = output.getFirstItem()

        assertHandle(
            self.assertEqual,
            genome.getScientificName(),
            'Homo sapiens',
            cwd=prot.getWorkingDir()
        )

        assertGenomeMetadata(
            self,
            prot,
            genome,
            'NCBI'
        )

        assertGenomeFiles(
            self,
            prot,
            genome,
            expectGtf=True
        )

        assertResolvedGenome(
            self,
            prot,
            genome
        )

        assertNcbiAccession(
            self,
            prot,
            genome
        )

    def testMultipleEnsemblGenomes(self):
        """Download multiple common genomes from Ensembl."""
        print("\nReference genomes: multiple Ensembl genomes")

        prot = self.newProtocol(
            ProtDownloadGenomes,
            source=ProtDownloadGenomes.SOURCE_ENSEMBL,
            genomeSelection=ProtDownloadGenomes.GENOME_COMMON,
            commonSpecies=(
                'Saccharomyces cerevisiae;'
                'Caenorhabditis elegans'
            ),
            assemblies='Latest',
            releases='Latest',
            downloadAnnotation=True,
            overwrite=False
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'referenceGenomes',
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

        scientificNames = []

        for genome in output:
            scientificNames.append(
                genome.getScientificName()
            )

            assertGenomeMetadata(
                self,
                prot,
                genome,
                'Ensembl'
            )

            assertGenomeFiles(
                self,
                prot,
                genome,
                expectGtf=True
            )

            assertResolvedGenome(
                self,
                prot,
                genome
            )

        assertHandle(
            self.assertIn,
            'Saccharomyces cerevisiae',
            scientificNames,
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertIn,
            'Caenorhabditis elegans',
            scientificNames,
            cwd=prot.getWorkingDir()
        )

    def testNcbiCustomMultipleGenomes(self):
        """Download multiple custom genomes from NCBI."""
        print("\nReference genomes: multiple NCBI custom genomes")

        prot = self.newProtocol(
            ProtDownloadGenomes,
            source=ProtDownloadGenomes.SOURCE_NCBI,
            genomeSelection=ProtDownloadGenomes.GENOME_CUSTOM,
            customSpecies='Homo sapiens;Mus musculus',
            ncbiAssemblies='Latest',
            downloadAnnotation=False,
            overwrite=False
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'referenceGenomes',
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

        scientificNames = []

        for genome in output:
            scientificNames.append(
                genome.getScientificName()
            )

            assertGenomeMetadata(
                self,
                prot,
                genome,
                'NCBI'
            )

            assertGenomeFiles(
                self,
                prot,
                genome,
                expectGtf=False
            )

            assertResolvedGenome(
                self,
                prot,
                genome
            )

            assertNcbiAccession(
                self,
                prot,
                genome
            )

        assertHandle(
            self.assertIn,
            'Homo sapiens',
            scientificNames,
            cwd=prot.getWorkingDir()
        )

        assertHandle(
            self.assertIn,
            'Mus musculus',
            scientificNames,
            cwd=prot.getWorkingDir()
        )

    def testNcbiGenBankFallback(self):
        """Download a genome requiring the NCBI GenBank fallback."""
        print("\nReference genomes: NCBI GenBank fallback")

        prot = self.newProtocol(
            ProtDownloadGenomes,
            source=ProtDownloadGenomes.SOURCE_NCBI,
            genomeSelection=ProtDownloadGenomes.GENOME_CUSTOM,
            customSpecies='Oryza sativa',
            ncbiAssemblies='Latest',
            downloadAnnotation=False,
            overwrite=False
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'referenceGenomes',
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

        genome = output.getFirstItem()

        assertHandle(
            self.assertEqual,
            genome.getScientificName(),
            'Oryza sativa',
            cwd=prot.getWorkingDir()
        )

        assertGenomeMetadata(
            self,
            prot,
            genome,
            'NCBI'
        )

        assertGenomeFiles(
            self,
            prot,
            genome,
            expectGtf=False
        )

        assertResolvedGenome(
            self,
            prot,
            genome
        )

        assertNcbiAccession(
            self,
            prot,
            genome
        )

        assertHandle(
            self.assertTrue,
            genome.getRelease().startswith('GCA_'),
            cwd=prot.getWorkingDir()
        )