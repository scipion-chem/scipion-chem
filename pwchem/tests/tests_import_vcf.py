
# ***************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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

from pwchem.protocols import ProtImportVCF


class TestImportVCF(BaseTest):
    """Test importing VCF files."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls.dataset = DataSet.getDataSet(
            'genomics'
        )

        cls.vcfFile = cls.dataset.getFile(
            'mouseVcfFile'
        )

        cls.vcfIndex = cls.dataset.getFile(
            'mouseVcfIndex'
        )

    def testImportVCF(self):
        """Import a mouse VCF file and its index."""

        vcfsData = json.dumps([
            {
                'scientificName': 'Mus musculus',
                'assembly': 'GRCm39',
                'vcfFile': self.vcfFile,
                'indexFile': self.vcfIndex
            }
        ])

        protocol = self.newProtocol(
            ProtImportVCF,
            objLabel='Import mouse VCF',
            vcfsData=vcfsData
        )

        self.launchProtocol(
            protocol
        )

        # -------------------------------------------------------------
        # Check protocol output
        # -------------------------------------------------------------

        self.assertTrue(
            hasattr(protocol, 'outputVCFs'),
            'ProtImportVCF did not produce outputVCFs.'
        )

        outputVCFs = protocol.outputVCFs

        self.assertEqual(
            outputVCFs.getSize(),
            1,
            'ProtImportVCF should produce exactly one VCF.'
        )

        # -------------------------------------------------------------
        # Check imported VCF
        # -------------------------------------------------------------

        vcf = next(
            iter(outputVCFs)
        )

        vcfFile = vcf.getFileName()

        self.assertTrue(
            os.path.isfile(vcfFile),
            'Imported VCF file does not exist: {}'.format(
                vcfFile
            )
        )

        self.assertGreater(
            os.path.getsize(vcfFile),
            0,
            'Imported VCF file is empty.'
        )

        # -------------------------------------------------------------
        # Check imported index
        # -------------------------------------------------------------

        indexFile = vcf.getIndexFile()

        self.assertTrue(
            os.path.isfile(indexFile),
            'Imported VCF index does not exist: {}'.format(
                indexFile
            )
        )

        self.assertGreater(
            os.path.getsize(indexFile),
            0,
            'Imported VCF index is empty.'
        )

        # -------------------------------------------------------------
        # Check metadata
        # -------------------------------------------------------------

        self.assertEqual(
            vcf.getScientificName(),
            'Mus musculus'
        )

        self.assertEqual(
            vcf.getAssembly(),
            'GRCm39'
        )

        self.assertEqual(
            vcf.getSource(),
            'Imported'
        )

        self.assertEqual(
            vcf.getDatabase(),
            'Local'
        )

        self.assertEqual(
            vcf.getVariantType(),
            'known'
        )

        self.assertTrue(
            vcf.getFileName().endswith('.vcf.gz')
        )