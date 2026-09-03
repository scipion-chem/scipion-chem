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

from pyworkflow.tests import BaseTest, setupTestProject
from pwchem.protocols import ProtReferenceGenomes
from pwchem.utils import assertHandle


def assertOutputExists(test, protocol, output):
    assertHandle(
        test.assertIsNotNone,
        output,
        cwd=protocol.getWorkingDir()
    )


def assertGenomeMetadata(test, protocol, genome, source):
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


def assertGenomeFiles(test, protocol, genome, expectGtf=True):
    assertHandle(
        test.assertTrue,
        os.path.exists(genome.getFastaFile()),
        cwd=protocol.getWorkingDir()
    )

    if expectGtf:
        assertHandle(
            test.assertTrue,
            genome.hasGtfFile(),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            test.assertTrue,
            os.path.exists(genome.getGtfFile()),
            cwd=protocol.getWorkingDir()
        )


def assertResolvedGenome(test, protocol, genome):
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
    assertHandle(
        test.assertRegex,
        genome.getRelease(),
        r'^GCF_\d+\.\d+$',
        cwd=protocol.getWorkingDir()
    )


class TestReferenceGenomes(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def testEnsemblCommonGenome(self):
        print("\nReference genomes: Ensembl common genome")

        prot = self.newProtocol(
            ProtReferenceGenomes,
            source=ProtReferenceGenomes.SOURCE_ENSEMBL,
            genomeSelection=ProtReferenceGenomes.GENOME_COMMON,
            commonGenomes='saccharomyces_cerevisiae',
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
        print("\nReference genomes: NCBI common genome")

        prot = self.newProtocol(
            ProtReferenceGenomes,
            source=ProtReferenceGenomes.SOURCE_NCBI,
            genomeSelection=ProtReferenceGenomes.GENOME_COMMON,
            commonGenomes='homo_sapiens',
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
        print("\nReference genomes: multiple Ensembl genomes")

        prot = self.newProtocol(
            ProtReferenceGenomes,
            source=ProtReferenceGenomes.SOURCE_ENSEMBL,
            genomeSelection=ProtReferenceGenomes.GENOME_COMMON,
            commonGenomes=(
                'saccharomyces_cerevisiae;'
                'caenorhabditis_elegans'
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
        print("\nReference genomes: multiple NCBI custom genomes")

        prot = self.newProtocol(
            ProtReferenceGenomes,
            source=ProtReferenceGenomes.SOURCE_NCBI,
            genomeSelection=ProtReferenceGenomes.GENOME_CUSTOM,
            customGenomes='Homo sapiens;Mus musculus',
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