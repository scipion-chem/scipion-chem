# ***************************************************************************
# *
# * Authors:		 Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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

# Scipion em imports
from pyworkflow.tests import DataSet
import pyworkflow.tests as tests
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.protocols import ProtocolPharmacophoreGeneration, ProtocolPharmacophoreModification, ProtocolPharmacophoreFiltering
from pwchem.tests import TestExtractLigand, TestRDKitLigandPreparation
from pwchem.utils import assertHandle

modStr = '''REMOVE | Set 0 | PharmFeature 4 (Type: Aromatic. Coords: (-31.98, 10.50, -54.42). Radius: 1.00)
MODIFY | Set 0 | PharmFeature 2 TO: {"Type": "Acceptor", "Coords": "(-37.79, 9.96, -57.07)", "Radius": 2.0}'''

class TestPharmGeneration(TestExtractLigand):
	@classmethod
	def _runImportPDB(cls):
		protImportPDB = cls.newProtocol(
			ProtImportPdb,
			inputPdbData=0, pdbId='4erf')
		cls.launchProtocol(protImportPDB)
		cls.protImportPDB = protImportPDB

	@classmethod
	def _runGenPharm(cls, inProt):
		protGenPharm = cls.newProtocol(
			ProtocolPharmacophoreGeneration,
			eps=4.0, topClusters=1
		)
		protGenPharm.inputSmallMolecules.set(inProt)
		protGenPharm.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protGenPharm, wait=False)
		return protGenPharm

	def test(self):
		protExtract = self._runExtractLigand(self.protImportPDB)
		self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)

		protPharm = self._runGenPharm(inProt=protExtract)
		self._waitOutput(protPharm, 'outputPharmacophore', sleepTime=10)

		assertHandle(self.assertIsNotNone, getattr(protPharm, 'outputPharmacophore', None), cwd=protPharm.getWorkingDir())

class TestPharmModification(TestPharmGeneration):
	@classmethod
	def _runModPharm(cls, inProt):
		protModPharm = cls.newProtocol(
			ProtocolPharmacophoreModification,
			operationList=modStr
		)
		protModPharm.inputPharmacophores.append(inProt)
		protModPharm.inputPharmacophores[-1].setExtended('outputPharmacophore')

		cls.proj.launchProtocol(protModPharm, wait=False)
		return protModPharm

	def test(self):
		protExtract = self._runExtractLigand(self.protImportPDB)
		self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)

		protPharm = self._runGenPharm(inProt=protExtract)
		self._waitOutput(protPharm, 'outputPharmacophore', sleepTime=10)

		protModPharm = self._runModPharm(inProt=protPharm)
		self._waitOutput(protModPharm, 'outputPharmacophore', sleepTime=10)

		assertHandle(self.assertIsNotNone, getattr(protModPharm, 'outputPharmacophore', None), cwd=protModPharm.getWorkingDir())

class TestPharmFiltering(TestPharmModification, TestRDKitLigandPreparation):
	@classmethod
	def setUpClass(cls):
			tests.setupTestProject(cls)
			cls.ds = DataSet.getDataSet('model_building_tutorial')
			cls.dsLig = DataSet.getDataSet("smallMolecules")
			cls._runImportPDB()
			cls._runImportSmallMols()
			cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

	@classmethod
	def _runFilterPharm(cls, inPharmProt, inMolsProt):
		protFilterPharm = cls.newProtocol(
			ProtocolPharmacophoreFiltering,
		)
		protFilterPharm.inputPharmacophore.set(inPharmProt)
		protFilterPharm.inputPharmacophore.setExtended('outputPharmacophore')

		protFilterPharm.inputSmallMolecules.set(inMolsProt)
		protFilterPharm.inputSmallMolecules.setExtended('outputSmallMolecules')

		cls.proj.launchProtocol(protFilterPharm, wait=False)
		return protFilterPharm

	def test(self):
		pPrep = self._runPrep(inProt=self.protImportSmallMols, mode=0)
		protExtract = self._runExtractLigand(self.protImportPDB)

		self._waitOutput(pPrep, 'outputSmallMolecules', sleepTime=10)
		self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)

		protPharm = self._runGenPharm(inProt=protExtract)
		self._waitOutput(protPharm, 'outputPharmacophore', sleepTime=10)

		protModPharm = self._runModPharm(inProt=protPharm)
		self._waitOutput(protModPharm, 'outputPharmacophore', sleepTime=10)

		protFilPharm = self._runFilterPharm(inPharmProt=protModPharm, inMolsProt=pPrep)
		self._waitOutput(protFilPharm, 'outputSmallMolecules', sleepTime=10)

		assertHandle(self.assertIsNotNone, getattr(protFilPharm, 'outputSmallMolecules', None), cwd=protFilPharm.getWorkingDir())
