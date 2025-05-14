# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_EXPORT_CSV.PY
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
import os

# Scipion em imports
from pyworkflow.tests import *

# Scipion chem imports
from pwchem.protocols import *
from pwchem.tests import TestImportBase
from pwchem.constants import *
from pwchem.utils import assertHandle

ZINCFiltStr = '''Remove if in not-for-sale
Keep if in investigational-only'''

class TestImportDBIDs(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)

	@classmethod
	def getIdsFilePath(cls, inType=0, dbName='UniProt'):
		return os.path.abspath(cls.proj.getPath('inputIDs_{}_{}.txt'.format(inType, dbName)))

	@classmethod
	def _runImportDBIds(cls, inType=0, dbName='UniProt'):
		idsFile = cls.getIdsFilePath(inType, dbName)
		with open(idsFile, 'w') as f:
			f.write('\n'.join(DB_IDS[inType][dbName]))

		protImport = cls.newProtocol(
			ProtChemImportSetOfDatabaseIDs,
			filePath=idsFile, databaseName=dbName
		)
		cls.proj.launchProtocol(protImport)
		return protImport

	def test(self):
		pImp = self._runImportDBIds()

		self._waitOutput(pImp, 'outputDatabaseIDs')
		assertHandle(self.assertIsNotNone, getattr(pImp, 'outputDatabaseIDs', None), cwd=pImp.getWorkingDir())

class TestIdentifyLigands(TestImportBase):
	@classmethod
	def _runIdentify(cls, inProt):
		protIdentify = cls.newProtocol(
			ProtChemSmallMolIdentify,
			useManager=1, nameDatabase=3,
		)

		protIdentify.inputSet.set(inProt)
		protIdentify.inputSet.setExtended('outputSmallMolecules')
		cls.proj.launchProtocol(protIdentify)
		return protIdentify

	def test(self):
		pIdent = self._runIdentify(self.protImportSmallMols)
		self._waitOutput(pIdent, 'outputSmallMolecules')
		assertHandle(self.assertIsNotNone, getattr(pIdent, 'outputSmallMolecules', None), cwd=pIdent.getWorkingDir())

class TestZINCFilter(TestIdentifyLigands):
	@classmethod
	def _runZINCFilter(cls, inProt):
		protFilter = cls.newProtocol(
			ProtChemZINCFilter,
			filterList=ZINCFiltStr
		)

		protFilter.inputSet.set(inProt)
		protFilter.inputSet.setExtended('outputSmallMolecules')
		cls.proj.launchProtocol(protFilter)
		return protFilter

	def test(self):
		pIdent = self._runIdentify(self.protImportSmallMols)
		self._waitOutput(pIdent, 'outputSmallMolecules')
		pFilt = self._runZINCFilter(pIdent)
		self._waitOutput(pFilt, 'outputSmallMolecules')
		assertHandle(self.assertIsNotNone, getattr(pFilt, 'outputSmallMolecules', None), cwd=pFilt.getWorkingDir())

class TestUniProtCrossRef(TestImportDBIDs):
	@classmethod
	def _runCrossRef(cls, inProt):
		protFilter = cls.newProtocol(
			ProtChemUniprotCrossRef,
			filterList=ZINCFiltStr
		)

		protFilter.inputListID.set(inProt)
		protFilter.inputListID.setExtended('outputDatabaseIDs')
		cls.proj.launchProtocol(protFilter)
		return protFilter

	def test(self):
		pImp = self._runImportDBIds()
		self._waitOutput(pImp, 'outputDatabaseIDs')
		pCrossRef = self._runCrossRef(pImp)
		self._waitOutput(pCrossRef, 'outputDatabaseIDs')
		assertHandle(self.assertIsNotNone, getattr(pCrossRef, 'outputDatabaseIDs', None), cwd=pCrossRef.getWorkingDir())

class TestFetchLigands(TestImportDBIDs):
	@classmethod
	def _runFetchLigands(cls, inProt, inType=0, inDataBase=0, structDataBase=0):
		protFetch = cls.newProtocol(
			ProtocolLigandsFetching,
			inputType=inType,
		)

		protFetch.inputIDs.set(inProt)
		protFetch.inputIDs.setExtended('outputDatabaseIDs')

		if inType == 0:
			protFetch.inputDatabaseUniprot.set(inDataBase)
			protFetch.structDatabase.set(structDataBase)
		elif inType == 1:
			protFetch.inputDatabaseTarget.set(inDataBase)
		elif inType == 2:
			protFetch.inputDatabaseLigand.set(inDataBase)
			protFetch.structDatabase.set(structDataBase)

		cls.proj.launchProtocol(protFetch)
		return protFetch

	def test(self):
		impProts, fetchProts = {}, []
		for inType, inDic in enumerate(DB_IDS):
			for iBase, dbName in enumerate(inDic):
				impProts[(inType, iBase)] = self._runImportDBIds(inType=inType, dbName=dbName)

		for p in impProts.values():
			self._waitOutput(p, 'outputDatabaseIDs')
			assertHandle(self.assertIsNotNone, getattr(p, 'outputDatabaseIDs', None), cwd=p.getWorkingDir())

		for inType, inDic in enumerate(DB_IDS):
			for iBase, dbName in enumerate(inDic):
				pImp = impProts[(inType, iBase)]
				fetchProts.append(self._runFetchLigands(pImp, inType=inType, inDataBase=iBase))

		for p in fetchProts:
			self._waitOutput(p, 'outputSmallMolecules')
			assertHandle(self.assertIsNotNone, getattr(p, 'outputSmallMolecules', None), cwd=p.getWorkingDir())
