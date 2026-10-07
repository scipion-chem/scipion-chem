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

import json
import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import ProtImportGenomes
from pwchem.utils import assertHandle
from pwchem.utils.utilsRNA import (
    assertOutputExists,
    assertGenomeFiles
)


def assertGenomeMetadata(
        test,
        protocol,
        genome,
        scientificName,
        assembly,
        release
):
    """Check the metadata stored in an imported Genome object."""
    assertHandle(
        test.assertEqual,
        genome.getScientificName(),
        scientificName,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        genome.getAssembly(),
        assembly,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        genome.getRelease(),
        release,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        genome.getSource(),
        'Imported',
        cwd=protocol.getWorkingDir()
    )


def assertFilesInExtra(test, protocol, genome):
    """Check that imported genome files are stored in the extra directory."""
    extraPath = os.path.abspath(
        protocol._getExtraPath()
    )

    fastaPath = os.path.abspath(
        genome.getFastaFile()
    )

    gtfPath = os.path.abspath(
        genome.getGtfFile()
    )

    assertHandle(
        test.assertEqual,
        os.path.commonpath([
            extraPath,
            fastaPath
        ]),
        extraPath,
        cwd=protocol.getWorkingDir()
    )

    assertHandle(
        test.assertEqual,
        os.path.commonpath([
            extraPath,
            gtfPath
        ]),
        extraPath,
        cwd=protocol.getWorkingDir()
    )


class TestImportGenomes(BaseTest):
    """Test the import of reference genomes from local FASTA/GTF files."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet(
            'genomics'
        )

        cls.humanFasta = cls.dataset.getFile(
            'human/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
        )
        cls.humanGtf = cls.dataset.getFile(
            'human/Homo_sapiens.GRCh38.116.gtf'
        )

        cls.mouseFasta = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.dna.primary_assembly.fa'
        )
        cls.mouseGtf = cls.dataset.getFile(
            'mouse/Mus_musculus.GRCm39.116.gtf'
        )

    def testImportMultipleGenomes(self):
        """Import human and mouse reference genomes."""
        print(
            "\nImport genomes: multiple reference genomes"
        )

        genomesData = [
            {
                'scientificName': 'Homo sapiens',
                'assembly': 'GRCh38',
                'release': '116',
                'fastaFile': self.humanFasta,
                'gtfFile': self.humanGtf
            },
            {
                'scientificName': 'Mus musculus',
                'assembly': 'GRCm39',
                'release': '116',
                'fastaFile': self.mouseFasta,
                'gtfFile': self.mouseGtf
            }
        ]

        prot = self.newProtocol(
            ProtImportGenomes,
            genomesData=json.dumps(
                genomesData
            )
        )

        self.launchProtocol(
            prot
        )

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

        # Scipion sets may reuse the same object instance while iterating,
        # so Genome objects are validated directly during iteration.
        for genome in output:
            scientificName = genome.getScientificName()

            scientificNames.append(
                scientificName
            )

            if scientificName == 'Homo sapiens':
                assertGenomeMetadata(
                    self,
                    prot,
                    genome,
                    'Homo sapiens',
                    'GRCh38',
                    '116'
                )

                assertGenomeFiles(
                    self,
                    prot,
                    genome,
                    expectGtf=True
                )

                assertFilesInExtra(
                    self,
                    prot,
                    genome
                )

                assertHandle(
                    self.assertEqual,
                    os.path.basename(
                        genome.getFastaFile()
                    ),
                    'Homo_sapiens_GRCh38_genome.fa',
                    cwd=prot.getWorkingDir()
                )

                assertHandle(
                    self.assertEqual,
                    os.path.basename(
                        genome.getGtfFile()
                    ),
                    'Homo_sapiens_GRCh38_annotation.gtf',
                    cwd=prot.getWorkingDir()
                )

            elif scientificName == 'Mus musculus':
                assertGenomeMetadata(
                    self,
                    prot,
                    genome,
                    'Mus musculus',
                    'GRCm39',
                    '116'
                )

                assertGenomeFiles(
                    self,
                    prot,
                    genome,
                    expectGtf=True
                )

                assertFilesInExtra(
                    self,
                    prot,
                    genome
                )

                assertHandle(
                    self.assertEqual,
                    os.path.basename(
                        genome.getFastaFile()
                    ),
                    'Mus_musculus_GRCm39_genome.fa',
                    cwd=prot.getWorkingDir()
                )

                assertHandle(
                    self.assertEqual,
                    os.path.basename(
                        genome.getGtfFile()
                    ),
                    'Mus_musculus_GRCm39_annotation.gtf',
                    cwd=prot.getWorkingDir()
                )

            else:
                self.fail(
                    'Unexpected genome in output: {}'.format(
                        scientificName
                    )
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