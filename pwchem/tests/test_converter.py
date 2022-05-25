# **************************************************************************
# *
# * Name:     test of protocol_converter.py
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import glob
from pathlib import Path

from pyworkflow.tests import *
from pwem.protocols.protocol_import import ProtImportPdb

from pwchem.protocols import ProtChemImportSmallMolecules as importSM

from pwchem.protocols import ConvertStructures as converter


class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        cls.lig_data = cls.dsLig.getFile('mix')

    @classmethod
    def _importPDB(cls):
        inputPdbData = 1  # file
        args = {'inputPdbData': 0,
                'pdbId': '4erf'
                }

        protocol = cls.newProtocol(ProtImportPdb, **args)
        cls.launchProtocol(protocol)
        pdb = protocol.outputPdb
        return pdb

    @classmethod
    def _importSmallM(cls, path):
        args = {'multiple': True,
                'filesPath': path,
                'filesPattern': '*'}

        protocol = cls.newProtocol(importSM, **args)
        cls.launchProtocol(protocol)
        setofSM = protocol.outputSmallMolecules
        return setofSM


class TestConverter(TestImportBase):


    def test_1(self):
        """ Covert a mix of small molecules file into pdb format
        """
        print("\nCovert a mix of small molecules file into pdb format \n")

        # Import SetOfSmallMolecules
        smallM = self._importSmallM(self.lig_data)

        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': smallM,
                'outputFormatSmall': 0,  # PDB
                }

        protocol = self.newProtocol(converter, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMolecules

        convert_file = glob.glob(protocol._getExtraPath("*"))


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==4,
                        "There was a problem with the import or conversion and the SetOfSmallMolecules is empty")

        files = ""
        for file in convert_file:
            if not file.endswith(".pdb"):
                files += "%s; "

        self.assertTrue(files == "",
                        "The conversion was incorrect and those files have a wrong format : %s" %files)


    def test_2(self):
        """ Covert a mix of small molecules file into smi format
        """
        print("\nCovert a mix of small molecules file into smi format \n")

        # Import SetOfSmallMolecules
        smallM = self._importSmallM(self.lig_data)

        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': smallM,
                'outputFormatSmall': 4,  # smiles or smi
                }

        protocol = self.newProtocol(converter, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMolecules

        convert_file = glob.glob(protocol._getExtraPath("*"))


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==4,
                        "There was a problem with the import or conversion and the SetOfSmallMolecules is empty")

        files = ""
        for file in convert_file:
            if not file.endswith(".smi"):
                files += "%s; "

        self.assertTrue(files == "",
                        "The conversion was incorrect and those files have a wrong format : %s" %files)


    def test_3(self):
        """ Covert a pdb protein file into mol2 format
        """
        print("\nCovert a pdb protein file into mol2 format \n")

        # Import PDB as Scipion object
        target = self._importPDB()

        args = {'inputType': 1,  # Protein structure
                'inputStructure': target,
                'outputFormatTarget': 2,  # smiles or smi
                }

        protocol = self.newProtocol(converter, **args)
        self.launchProtocol(protocol)
        prot = protocol.outputStructure


        self.assertIsNotNone(prot,
                             "There was a problem with the conversion and the new file not exist")

        prot_end = prot.getFileName()
        self.assertTrue(prot_end.endswith(".mol2"),
                        "The conversion was incorrect and this file have a wrong format : %s" %os.path.basename(prot_end))