# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *
# **************************************************************************

# General imports
import os, glob

# Scipion em imports
from pyworkflow.tests import BaseTest, setupTestProject, DataSet
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.tests.tests_imports import TestImportBase
from pwchem.protocols import ProtChemImportSmallMolecules, ConvertStructures, ProtChemExportCSV, ProtChemOperateSet
from pwchem.utils import assertHandle

class TestImportBoth(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		cls.dsLig = DataSet.getDataSet("smallMolecules")
		cls.lig_data = cls.dsLig.getFile('mix')
		cls.pdbProt, cls.smallProt = cls._importPDB(), cls._importSmallM(cls.lig_data)

	@classmethod
	def _importPDB(cls):
		inputPdbData = 1 # file
		args = {
			 'inputPdbData': inputPdbData,
				'pdbFile': cls.ds.getFile('PDBx_mmCIF/5ni1.pdb')
			}

		protocol = cls.newProtocol(ProtImportPdb, **args)
		cls.launchProtocol(protocol)
		return protocol

	@classmethod
	def _importSmallM(cls, path):
		args = {
			'multiple': True,
			'filesPath': path,
			'filesPattern': '*'
		}

		protocol = cls.newProtocol(ProtChemImportSmallMolecules, **args)
		cls.launchProtocol(protocol)
		return protocol

class TestConverter(TestImportBoth):
	def test_1(self):
		""" Convert a mix of small molecules file into pdb format. """
		print("\nConvert a mix of small molecules file into pdb format \n")

		# Import SetOfSmallMolecules
		smallM = self.smallProt.outputSmallMolecules
		args = {
			'inputObject': smallM, # SmallMolecules
			"useManager": 1,
			'outputFormatSmall': 0, # PDB
		}

		protocol = self.newProtocol(ConvertStructures, **args)
		self.launchProtocol(protocol)

		small1 = getattr(protocol, 'outputSmallMolecules', None)
		convertFile = glob.glob(protocol._getExtraPath("*"))

		assertHandle(self.assertIsNotNone, small1, message="There was a problem with the import", cwd=protocol.getWorkingDir())
		assertHandle(self.assertTrue, small1.getSize()==4,
								 message="There was a problem with the import or conversion and the SetOfSmallMolecules is empty", cwd=protocol.getWorkingDir())

		files = ""
		for file in convertFile:
			if not file.endswith(".pdb"):
				files += "%s; "

		assertHandle(self.assertTrue, files == "",
								 message="The conversion was incorrect and those files have a wrong format : %s" %files, cwd=protocol.getWorkingDir())

	def test_2(self):
		""" Convert a mix of small molecules file into smi format. """
		print("\nConvert a mix of small molecules file into smi format \n")

		# Import SetOfSmallMolecules
		smallM = self.smallProt.outputSmallMolecules

		args = {
			'inputObject': smallM, # SmallMolecules
			"useManager": 1,
			'outputFormatSmall': 3, # smiles or smi
		}

		protocol = self.newProtocol(ConvertStructures, **args)
		self.launchProtocol(protocol)
		small1 = getattr(protocol, 'outputSmallMolecules', None)
		convertFile = glob.glob(protocol._getExtraPath("*"))

		assertHandle(self.assertIsNotNone, small1, message="There was a problem with the import", cwd=protocol.getWorkingDir())
		assertHandle(self.assertTrue, small1.getSize()==4,
								 message="There was a problem with the import or conversion and the SetOfSmallMolecules is empty", cwd=protocol.getWorkingDir())

		files = ""
		for file in convertFile:
			if not file.endswith(".smi"):
				files += "%s; "

		assertHandle(self.assertTrue, files == "",
								 message="The conversion was incorrect and those files have a wrong format : %s" %files, cwd=protocol.getWorkingDir())

	def test_3(self):
		""" Convert a cif protein file into pdb format. """
		print("\nConvert a pdb protein file into mol2 format \n")

		# Import PDB as Scipion object
		target = self.pdbProt.outputPdb

		args = {
			'inputObject': target, # AtomStruct
			'outputFormatTarget': 0, # pdb
		}

		protocol = self.newProtocol(ConvertStructures, **args)
		self.launchProtocol(protocol)
		prot = getattr(protocol, 'outputStructure', None)

		assertHandle(self.assertIsNotNone, prot,
								 message="There was a problem with the conversion and the new file not exist", cwd=protocol.getWorkingDir())

		protEnd = prot.getFileName()
		assertHandle(self.assertTrue, protEnd.endswith(".pdb"),
								 message="The conversion was incorrect and this file have a wrong format : %s" %os.path.basename(protEnd), cwd=protocol.getWorkingDir())

class TestExportcsv(TestImportBase):
	@classmethod
	def _runExport(cls, inProt):
		protExport = cls.newProtocol(
				ProtChemExportCSV,
		)
		protExport.inputSet.set(inProt)
		protExport.inputSet.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protExport, wait=True)
		return protExport

	def test(self):
		pExp = self._runExport(inProt=self.protImportSmallMols)

		csvFile = os.path.abspath(pExp._getPath('output.csv'))
		assertHandle(self.assertTrue, os.path.exists(csvFile),
								 message="CSV file was not create. Check if its location changed", cwd=pExp.getWorkingDir())

class TestOperateSet(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.dsLig = DataSet.getDataSet("smallMolecules")

		cls.smallProtAll, cls.smallProtPart = cls._importSmallM(cls.dsLig.getFile('mol2')), \
																					cls._importSmallM(cls.dsLig.getFile('mol2'), pattern='*9.mol2')

	@classmethod
	def _importSmallM(cls, path, pattern='*'):
		args = {
			'multiple': True,
			'filesPath': path,
			'filesPattern': pattern
		}

		protocol = cls.newProtocol(ProtChemImportSmallMolecules, **args)
		cls.launchProtocol(protocol)
		return protocol

	@classmethod
	def _runSingleInputOperation(cls, inProt, op, refCol, extended='outputSmallMolecules', remDup=False, **kwargs):
		protOp = cls.newProtocol(
			ProtChemOperateSet,
			operation=op, refColumn=refCol, removeDuplicates=remDup, **kwargs
		)
		protOp.inputSet.set(inProt)
		protOp.inputSet.setExtended(extended)

		cls.proj.launchProtocol(protOp, wait=False)
		return protOp

	@classmethod
	def _runMultiInputOperation(cls, inProts, op, refCol, extended='outputSmallMolecules', remDup=False, **kwargs):
		protOp = cls.newProtocol(
			ProtChemOperateSet,
			operation=op, refColumn=refCol, removeDuplicates=remDup, **kwargs
		)
		for p in inProts:
			protOp.inputMultiSet.append(p)
			protOp.inputMultiSet[-1].setExtended(extended)

		cls.proj.launchProtocol(protOp, wait=False)
		return protOp

	@classmethod
	def _runDoubleInputOperation(cls, inProts, op, refCol, extended='outputSmallMolecules', remDup=False, **kwargs):
		protOp = cls.newProtocol(
			ProtChemOperateSet,
			operation=op, refColumn=refCol, removeDuplicates=remDup, **kwargs
		)
		protOp.inputSet.set(inProts[0])
		protOp.inputSet.setExtended(extended)

		protOp.secondSet.set(inProts[1])
		protOp.secondSet.setExtended(extended)

		cls.proj.launchProtocol(protOp, wait=False)
		return protOp

	def test(self):
		"""
		1. Union of several sets + 0. Unique operation
		2. Intersection of several sets
		3. Difference between a first and a second set
		4. Filter a set by the values of the items for one of the columns
		5. Remove specified columns
		6. Ranking according to reference column: higher / lower values
		"""
		inProts = [self.smallProtAll, self.smallProtPart]
		outProts = []

		protUnion = self._runMultiInputOperation(inProts, op=1, refCol='molName')
		outProts.append(protUnion)

		self._waitOutput(protUnion, 'outputSet', sleepTime=5)
		assertHandle(self.assertIsNotNone, getattr(protUnion, 'outputSet', None), cwd=protUnion.getWorkingDir())
		protUniq = self._runSingleInputOperation(protUnion, op=0, refCol='molName', extended='outputSet')
		outProts.append(protUniq)

		protInter = self._runMultiInputOperation(inProts, op=2, refCol='molName')
		outProts.append(protInter)

		protDiff = self._runDoubleInputOperation(inProts, op=3, refCol='molName')
		outProts.append(protDiff)

		kw = {'filterColumn': 'smallMoleculeFile', 'filterOp': 9, 'filterValue': '9.mol2'}
		protFilt = self._runSingleInputOperation(inProts[0], op=4, refCol='molName', **kw)
		outProts.append(protFilt)

		kw = {'remColumns': 'confId;gridId;poseId;dockId'}
		protRem = self._runSingleInputOperation(inProts[0], op=5, refCol='molName', **kw)
		outProts.append(protRem)

		kw = {'filterColumn': 'molName', 'threshold': '2'}
		protRank = self._runSingleInputOperation(inProts[0], op=6, refCol='molName', **kw)
		outProts.append(protRank)

		for p in outProts:
			self._waitOutput(p, 'outputSet', sleepTime=5)
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSet', None), cwd=p.getWorkingDir())
