# ***************************************************************************
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

# Scipion em imports
from pyworkflow.tests import *
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.protocols import *
from pwchem.tests.tests_imports import TestImportBase
from pwchem.utils import assertHandle

class TestOBLigandPreparation(TestImportBase):
	def test_1(self):
		""" Prepare a set of 4 ligands with conformer generation (genetic algotithm). """
		print("\n Prepare a set of 4 ligands with conformer generation (genetic algotithm) \n")

		inputMols = self.protImportSmallMols.outputSmallMolecules
		args = {
			'inputType': 0, # SmallMolecules
			'inputSmallMolecules': inputMols,
			'method_charges': 0, # gasteiger
			'doConformers': True,
			"method_conf": 0,
			"number_conf": 10,
			"rmsd_cutoff": 0.375,
		}

		protocol = self.newProtocol(ProtChemOBabelPrepareLigands, **args)
		self.launchProtocol(protocol)
		small1 = getattr(protocol, 'outputSmallMolecules', None)

		assertHandle(self.assertIsNotNone, small1, message="There was a problem with the import", cwd=protocol.getWorkingDir())

		for mol in small1:
			assertHandle(self.assertTrue, (mol.getFileName()).endswith(".mol2"),
										message="The format of first molecule is wrong. It must be in mol2 format", cwd=protocol.getWorkingDir())
			try:
				assertHandle(self.assertTrue, (mol.getConformersFileName()).endswith("_conformers.mol2"),
										message="The format of conformers molecules is wrong. It must be in mol2 format", cwd=protocol.getWorkingDir())
			except:
				assertHandle(self.assertTrue, (mol.getConformersFileName()).endswith("Not available"),
										message="Something was wrong in column of _ConformessFile", cwd=protocol.getWorkingDir())

	def test_2(self):
		""" Prepare a set of 4 ligands with conformer generation (confab algotithm). """
		print("\n Prepare a set of 4 ligands with conformer generation (confab algotithm) \n")

		inputMols = self.protImportSmallMols.outputSmallMolecules
		args = {
			'inputType': 0, # SmallMolecules
			'inputSmallMolecules': inputMols,
			'method_charges': 1, # mmff94
			'doConformers': True,
			"method_conf": 1,
			"number_conf": 10,
			"rmsd_cutoff": 0.5,
		}

		protocol = self.newProtocol(ProtChemOBabelPrepareLigands, **args)
		self.launchProtocol(protocol)
		small1 = getattr(protocol, 'outputSmallMolecules', None)

		assertHandle(self.assertIsNotNone, small1, message="There was a problem with the import", cwd=protocol.getWorkingDir())

		for mol in small1:
			assertHandle(self.assertTrue, (mol.getFileName()).endswith(".mol2"),
										message="The format of first molecule is wrong. It must be in mol2 format", cwd=protocol.getWorkingDir())
			try:
				assertHandle(self.assertTrue, (mol.getConformersFileName()).endswith("_conformers.mol2"),
										message="The format of conformers molecules is wrong. It must be in mol2 format", cwd=protocol.getWorkingDir())
			except:
				assertHandle(self.assertTrue, (mol.getConformersFileName()).endswith("Not available"),
										message="Something was wrong in column of _ConformessFile", cwd=protocol.getWorkingDir())

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
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())

class TestPrepareReceptor(BaseTest):
	@classmethod
	def setUpClass(cls):
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		setupTestProject(cls)
		cls._runImportPDB()
		cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

	@classmethod
	def _runImportPDB(cls):
		cls.protImportPDB = cls.newProtocol(
			ProtImportPdb,
			inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
		cls.proj.launchProtocol(cls.protImportPDB, wait=False)

	@classmethod
	def _runPrepareReceptor(cls):
		cls.protPrepareReceptor = cls.newProtocol(
			ProtChemPrepareReceptor,
			inputAtomStruct=cls.protImportPDB.outputPdb,
			HETATM=True, rchains=True,
			chain_name='{"model": 0, "chain": "C", "residues": 141}')

		cls.launchProtocol(cls.protPrepareReceptor)

	def test(self):
		self._runPrepareReceptor()

		self._waitOutput(self.protPrepareReceptor, 'outputStructure', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(self.protPrepareReceptor, 'outputStructure', None), cwd=self.protPrepareReceptor.getWorkingDir())
