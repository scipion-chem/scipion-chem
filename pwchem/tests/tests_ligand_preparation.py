# ***************************************import numpy as np***********************************
# *
# * Name:     test of protocol_ligand_preparation.py
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

from pyworkflow.tests import *
from pwchem.protocols import *

from .tests_imports import TestImportBase


class TestOBLigandPreparation(TestImportBase):

    def test_1(self):
        """ Prepare a set of 4 ligands with conformer generation (genetic algotithm)
        """
        print("\n Prepare a set of 4 ligands with conformer generation (genetic algotithm) \n")

        inputMols = self.protImportSmallMols.outputSmallMolecules
        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': inputMols,
                'method_charges': 0,  # gasteiger
                'doConformers': True,
                "method_conf": 0,
                "number_conf": 10,
                "rmsd_cutoff": 0.375,
                }

        protocol = self.newProtocol(ProtChemOBabelPrepareLigands, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMolecules


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")

        for mol in small_1:
            self.assertTrue((mol.getFileName()).endswith(".mol2"),
                            "The format of first molecule is wrong. It must be in mol2 format")
            try:
                self.assertTrue((mol.getConformersFileName()).endswith("_conformers.mol2"),
                            "The format of conformers molecules is wrong. It must be in mol2 format")


            except:
                self.assertTrue((mol.getConformersFileName()).endswith("Not available"),
                                "Something was wrong in column of _ConformessFile")


    def test_2(self):
        """ Prepare a set of 4 ligands with conformer generation (confab algotithm)
        """
        print("\n Prepare a set of 4 ligands with conformer generation (confab algotithm) \n")

        inputMols = self.protImportSmallMols.outputSmallMolecules
        args = {'inputType': 0,  # SmallMolecules
                'inputSmallMols': inputMols,
                'method_charges': 1,  # mmff94
                'doConformers': True,
                "method_conf": 1,
                "number_conf": 10,
                "rmsd_cutoff": 0.5,
                }

        protocol = self.newProtocol(ProtChemOBabelPrepareLigands, **args)
        self.launchProtocol(protocol)
        small_1 = protocol.outputSmallMolecules


        self.assertIsNotNone(small_1,
                             "There was a problem with the import")

        for mol in small_1:
            self.assertTrue((mol.getFileName()).endswith(".mol2"),
                            "The format of first molecule is wrong. It must be in mol2 format")
            try:
                self.assertTrue((mol.getConformersFileName()).endswith("_conformers.mol2"),
                                "The format of conformers molecules is wrong. It must be in mol2 format")

            except:
                self.assertTrue((mol.getConformersFileName()).endswith("Not available"),
                                "Something was wrong in column of _ConformersFile")


class TestRDKitLigandPreparation(TestImportBase):
  @classmethod
  def _runPrep(cls, inProt, mode=0):
    protShape = cls.newProtocol(
      ProtChemRDKitPrepareLigands,
    )
    protShape.inputSmallMolecules.set(inProt)
    protShape.inputSmallMolecules.setExtended('outputSmallMolecules')
    if mode == 1:
        protShape.doConformers.set(True)

    cls.proj.launchProtocol(protShape, wait=False)
    return protShape

  def test(self):
    protsShape = []
    for i in range(2):
        protsShape.append(self._runPrep(inProt=self.protImportSmallMols, mode=i))

    for p in protsShape:
        self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
        self.assertIsNotNone(getattr(p, 'outputSmallMolecules', None))