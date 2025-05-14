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
	ProtocolShapeDistancesFiltering, ProtocolFingerprintFiltering
from pwchem.tests.tests_imports import TestImportBase
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