# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_IMPORT_SMALLMOLECULES.PY
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

import os
from pyworkflow.tests import *
from pwchem.protocols import ProtChemImportSmallMolecules as PISmallM


class TestImportBase(BaseTest):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet("smallMolecules")

class TestImportSmallMolecules(TestImportBase):
    MULTIPLE_T = True
    MULTIPLE_F = False
    filesPattern = '*'

    def testImport_one_mol_smi(self):
        """ Import a single file of a small molecule provided by the user in smi format
        """
        print("\nImport Experiment:  1 molecule smi format")

        args = {'multiple': self.MULTIPLE_F,
                'filePath': self.dsModBuild.getFile(os.path.join("mix", "2000_smi.smi"))
                }

        prot1 = self.newProtocol(PISmallM, **args)
        self.launchProtocol(prot1)
        small_1 = prot1.outputSmallMolecules
        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==1,
                        "There was a problem with the import and the SetOfSmallMolecules is empty")

    def testImport_one_mol_pdb(self):
        """ Import a single file of a small molecule provided by the user in pdb format
        """
        print("\nImport Experiment:  1 molecule pdb format")

        args = {'multiple': self.MULTIPLE_F,
                'filePath': self.dsModBuild.getFile(os.path.join("mix", "2000_pdb.pdb"))
                }

        prot1 = self.newProtocol(PISmallM, **args)
        self.launchProtocol(prot1)
        small_1 = prot1.outputSmallMolecules
        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==1,
                        "There was a problem with the import and the SetOfSmallMolecules is empty")

    def testImport_one_mol_sdf(self):
        """ Import a single file of a small molecule provided by the user in sdf format
        """
        print("\nImport Experiment:  1 molecule sdf format")

        args = {'multiple': self.MULTIPLE_F,
                'filePath': self.dsModBuild.getFile(os.path.join("mix", "2000_sdf.sdf"))
                }

        prot1 = self.newProtocol(PISmallM, **args)
        self.launchProtocol(prot1)
        small_1 = prot1.outputSmallMolecules
        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==1,
                        "There was a problem with the import and the SetOfSmallMolecules is empty")


    def testImport_mols_mix(self):
        """Import several files of a small molecule provided by the user in different formats
        """
        print("\nImport Experiment:  4 molecules in several formats")

        args = {'multiple': self.MULTIPLE_T,
                'filesPath': self.dsModBuild.getFile(os.path.join("mix")),
                'filesPattern': self.filesPattern
                }

        prot1 = self.newProtocol(PISmallM, **args)
        self.launchProtocol(prot1)
        small_1 = prot1.outputSmallMolecules
        self.assertIsNotNone(small_1,
                             "There was a problem with the import")
        self.assertTrue(small_1.getSize()==4,
                        "There was a problem with the import and the SetOfSmallMolecules is empty")

