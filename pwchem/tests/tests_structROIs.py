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

# Scipion em imports
from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.protocols import ProtDefineStructROIs, ProtocolConsensusStructROIs
from pwchem.utils import assertHandle

defROIsStr = '''1) Coordinate: {"X": 45, "Y": 65, "Z": 60}
2) Residues: {"model": 0, "chain": "A", "index": "1-4", "residues": "VLSP"}
3) Ligand: {"molName": "HEM", "remove": "True"}
4) PPI: {"chain1": "0-A", "chain2": "0-B", "interDist": "5.0"}
5) Near_Residues: {"residues": "cys, cys", "distance": "5.0", "linkage": "Single"}'''

defROIsStr2 = '''1) Residues: {"model": 0, "chain": "A", "index": "61-63", "residues": "KVA"}
2) Residues: {"model": 0, "chain": "A", "index": "80-84", "residues": "LSALS"}
3) Ligand: {"molName": "HEM", "remove": "True"}'''

class TestDefineStructROIs(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		cls._runImportPDB()

	@classmethod
	def _runImportPDB(cls):
		protImportPDB = cls.newProtocol(
			ProtImportPdb,
			inputPdbData=1,
			pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
		cls.launchProtocol(protImportPDB)
		cls.protImportPDB = protImportPDB

	@classmethod
	def _runDefStructROIs(cls, defStr):
		protDef = cls.newProtocol(
			ProtDefineStructROIs,
			inROIs=defStr)

		protDef.inputAtomStruct.set(cls.protImportPDB)
		protDef.inputAtomStruct.setExtended('outputPdb')

		cls.proj.launchProtocol(protDef)
		return protDef

	def test(self):
		pDef = self._runDefStructROIs(defROIsStr)
		self._waitOutput(pDef, 'outputStructROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(pDef, 'outputStructROIs', None), cwd=pDef.getWorkingDir())

class TestConsensusStructROIs(TestDefineStructROIs):
	def _runConsensusStructROIs(self, inputSetsOfPockets):
		prot = self.newProtocol(ProtocolConsensusStructROIs)

		prot.inputStructROIsSets.set(inputSetsOfPockets)
		for i in range(len(prot.inputStructROIsSets)):
			prot.inputStructROIsSets[i].setExtended('outputStructROIs')

		self.launchProtocol(prot)
		return prot

	def test(self):
		pDef1 = self._runDefStructROIs(defROIsStr)
		pDef2 = self._runDefStructROIs(defROIsStr2)
		self._waitOutput(pDef1, 'outputStructROIs', sleepTime=10)
		self._waitOutput(pDef2, 'outputStructROIs', sleepTime=10)

		pCons = self._runConsensusStructROIs([pDef1, pDef2])
		self._waitOutput(pCons, 'outputStructROIs', sleepTime=10)
		assertHandle(self.assertIsNotNone, getattr(pCons, 'outputStructROIs', None),
	      					message="There was a problem with the consensus", cwd=pCons.getWorkingDir())
