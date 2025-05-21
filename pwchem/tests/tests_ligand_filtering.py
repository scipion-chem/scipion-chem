# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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

# Scipion chem imports
from pwchem.protocols import ProtocolADMEFiltering, ProtocolPainsRdkitFiltering, ProtocolGeneralLigandFiltering, \
	ProtocolShapeDistancesFiltering, ProtocolFingerprintFiltering, ProtocolOperateLibrary
from pwchem.tests.tests_imports import TestImportBase, TestImportSmallMoleculesLibrary
from pwchem.utils import assertHandle

generalFilter = '''Remove molecule if contains at least 1 atom type B
Keep molecule if contains at least 2 atoms 
Remove molecule if contains at least 3 cycles'''

class TestGeneralFiltering(TestImportBase):
	@classmethod
	def _runFilter(cls, inProt):
		protFilter = cls.newProtocol(
			ProtocolGeneralLigandFiltering,
			filterList=generalFilter
		)
		protFilter.inputSmallMolecules.set(inProt)
		protFilter.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protFilter, wait=False)
		return protFilter

	def test(self):
		protFilter = self._runFilter(inProt=self.protImportSmallMols)
		self._waitOutput(protFilter, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(protFilter, 'outputSmallMolecules', None),
								 cwd=protFilter.getWorkingDir())

class TestADMEFiltering(TestImportBase):
	@classmethod
	def _runADME(cls, inProt, mode=0):
		protADME = cls.newProtocol(
			ProtocolADMEFiltering,
			ruleChoice=mode
		)
		protADME.inputSmallMolecules.set(inProt)
		protADME.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protADME, wait=False)
		return protADME

	def test(self):
		protsADME = []
		for i in range(2):
			protsADME.append(self._runADME(inProt=self.protImportSmallMols, mode=i))

		for p in protsADME:
			self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())

class TestPAINSFiltering(TestImportBase):
	@classmethod
	def _runPAINS(cls, inProt):
		protPAINS = cls.newProtocol(
			ProtocolPainsRdkitFiltering,
		)
		protPAINS.inputSmallMolecules.set(inProt)
		protPAINS.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protPAINS, wait=False)
		return protPAINS

	def test(self):
		p = self._runPAINS(inProt=self.protImportSmallMols)

		self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())

class TestShapeFiltering(TestImportBase):
	@classmethod
	def _runShapeFilter(cls, inProt, mode=0):
		protShape = cls.newProtocol(
			ProtocolShapeDistancesFiltering,
			distanceType=mode, inputReferenceMolecule='SmallMolecule (ZINC00000480 molecule)'
		)
		protShape.inputSmallMolecules.set(inProt)
		protShape.inputSmallMolecules.setExtended('outputSmallMolecules')

		protShape.inputRefSmallMolecules.set(inProt)
		protShape.inputRefSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protShape, wait=False)
		return protShape

	def test(self):
		protsShape = []
		for i in range(3):
			protsShape.append(self._runShapeFilter(inProt=self.protImportSmallMols, mode=i))

		for p in protsShape:
			self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())

class TestFingerprintFiltering(TestImportBase):
	@classmethod
	def _runFingerprintFilter(cls, inProt, mode=0):
		protFinger = cls.newProtocol(
			ProtocolFingerprintFiltering,
			fpChoice=mode, inputReferenceMolecule='SmallMolecule (ZINC00000480 molecule)'
		)
		protFinger.inputSmallMolecules.set(inProt)
		protFinger.inputSmallMolecules.setExtended('outputSmallMolecules')

		protFinger.inputRefSmallMolecules.set(inProt)
		protFinger.inputRefSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protFinger, wait=False)
		return protFinger

	def test(self):
		protsShape = []
		for i in range(2):
			protsShape.append(self._runFingerprintFilter(inProt=self.protImportSmallMols, mode=i))

		for p in protsShape:
			self._waitOutput(p, 'outputSmallMolecules', sleepTime=10)
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())


class TestOperateSmallMoleculesLibrary(TestImportSmallMoleculesLibrary):

		@classmethod
		def _runOperateSets(cls, protLibs, op=0):
			protOperateLibrary = cls.newProtocol(ProtocolOperateLibrary, operation=op, refAttribute=1)
			for protLib in protLibs:
				protOperateLibrary.inputLibraries.append(protLib)
				protOperateLibrary.inputLibraries[-1].setExtended('outputLibrary')
	
			cls.proj.launchProtocol(protOperateLibrary, wait=False)
			return protOperateLibrary

		@classmethod
		def _runOperateLibrary(cls, protLib, op, kwargs):
			protOperateLibrary = cls.newProtocol(ProtocolOperateLibrary, operation=op, **kwargs)
			protOperateLibrary.inputLibrary.set(protLib)
			protOperateLibrary.inputLibrary.setExtended('outputLibrary')

			cls.proj.launchProtocol(protOperateLibrary, wait=False)
			return protOperateLibrary


		def test(self):
			protImport1 = self._runImportLibrary(defLib=True, nMols=1000, rSeed=44)
			protImport2 = self._runImportLibrary(defLib=True, nMols=1000, rSeed=4)
			self._waitOutput(protImport1, 'outputLibrary')
			self._waitOutput(protImport2, 'outputLibrary')
			
			opProts = []
			for op in range(3):
				opProts.append(self._runOperateSets([protImport1, protImport2], op))

			filtList = 'Keep molecule if LogP is below threshold 2.0\n'
			kwDic = {3: {"filterList": filtList}, 4: {"filterAttr": 'LogP', "filterValue": '-0.4'},
							 5: {"removeList": "Reactivity\nPurchasability"}}
			for op in range(3, 6):
				opProts.append(self._runOperateLibrary(protImport1, op, kwDic[op]))

			for opProt in opProts:
				self._waitOutput(opProt, 'outputLibrary')
				assertHandle(self.assertIsNotNone, getattr(opProt, 'outputLibrary', None),
										 cwd=opProt.getWorkingDir())
