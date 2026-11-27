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
import tempfile
import unittest
from unittest.mock import patch

from pyworkflow.tests import BaseTest, setupTestProject

from pwchem.protocols import ProtDownloadGenomes, ProtDownloadVCF
from pwchem.utils import assertHandle
from pwchem.utils.utilsRNA import assertOutputExists, assertGenomeFiles
from pwchem.utils.utilsSpecies import GenomicResourceMixin


class TestGenomicResourceHelpers(unittest.TestCase):
    """Offline regression checks for shared genome/VCF metadata helpers."""

    def testSharedMethodsAreInherited(self):
        for cls in (ProtDownloadGenomes, ProtDownloadVCF):
            self.assertTrue(issubclass(cls, GenomicResourceMixin))
            self.assertEqual(cls.splitParameterValues(' a; b,, '), ['a', 'b'])
            self.assertEqual(cls.expandParameterValues('Latest', 2, 'releases'), ['Latest', 'Latest'])
            self.assertEqual(cls.expandParameterValues('1;2', 2, 'releases'), ['1', '2'])
            with self.assertRaises(ValueError):
                cls.expandParameterValues('', 2, 'releases')
            with self.assertRaises(ValueError):
                cls.expandParameterValues('1;2;3', 2, 'releases')

    def testReleaseUsesProviderRequest(self):
        for cls in (ProtDownloadGenomes, ProtDownloadVCF):
            protocol = cls()
            with patch.object(protocol, '_requestJson', return_value={'releases': ['99', '116', '100']}) as request:
                self.assertEqual(protocol.getLatestEnsemblRelease(), 116)
                request.assert_called_once_with(protocol.ENSEMBL_REST_DATA_URL, 'latest Ensembl release')
            with patch.object(protocol, '_requestJson', return_value={}):
                with self.assertRaises(RuntimeError):
                    protocol.getLatestEnsemblRelease()

    def testLocalMetadataAndFileChecks(self):
        with tempfile.TemporaryDirectory() as directory:
            path = os.path.join(directory, 'metadata.jsonl')
            self.assertEqual(GenomicResourceMixin.readJsonLines(path), [])
            self.assertFalse(GenomicResourceMixin.isNonemptyFile(path))
            with open(path, 'w') as handle:
                handle.write('{"accession": "GCF_1.1"}\n\n{"accession": "GCF_2.1"}\n')
            self.assertEqual(GenomicResourceMixin.readJsonLines(path),
                             [{'accession': 'GCF_1.1'}, {'accession': 'GCF_2.1'}])
            self.assertTrue(GenomicResourceMixin.isNonemptyFile(path))
            self.assertEqual(GenomicResourceMixin.safeName('Homo sapiens/GRCh38'), 'Homo_sapiens_GRCh38')


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

        self.assertGenomeMetadata(
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

        self.assertResolvedGenome(
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

        self.assertGenomeMetadata(
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

        self.assertResolvedGenome(
            prot,
            genome
        )

        self.assertNcbiAccession(
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

            self.assertGenomeMetadata(
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

            self.assertResolvedGenome(
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

            self.assertGenomeMetadata(
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

            self.assertResolvedGenome(
                prot,
                genome
            )

            self.assertNcbiAccession(
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

        self.assertGenomeMetadata(
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

        self.assertResolvedGenome(
            prot,
            genome
        )

        self.assertNcbiAccession(
            prot,
            genome
        )

        assertHandle(
            self.assertTrue,
            genome.getRelease().startswith('GCA_'),
            cwd=prot.getWorkingDir()
        )

    def assertGenomeMetadata(self, protocol, genome, source):
        """Check the basic metadata stored in a downloaded Genome object."""
        assertHandle(
            self.assertTrue,
            bool(genome.getScientificName()),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            bool(genome.getAssembly()),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            bool(genome.getRelease()),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertEqual,
            genome.getSource(),
            source,
            cwd=protocol.getWorkingDir()
        )


    def assertResolvedGenome(self, protocol, genome):
        """Check that Latest values were resolved."""
        assertHandle(
            self.assertNotEqual,
            genome.getAssembly().lower(),
            'latest',
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertNotEqual,
            genome.getRelease().lower(),
            'latest',
            cwd=protocol.getWorkingDir()
        )


    def assertNcbiAccession(self, protocol, genome):
        """Check that the NCBI release is a versioned assembly accession."""
        assertHandle(
            self.assertRegex,
            genome.getRelease(),
            r'^(GCF|GCA)_\d+\.\d+$',
            cwd=protocol.getWorkingDir()
        )
