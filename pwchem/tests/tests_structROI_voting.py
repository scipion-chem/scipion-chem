# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
# *              Irene Sánchez Martín (100495638@alumnos.uc3m.es)
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

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportPdb

from pwchem.protocols import ProtDefineStructROIs, ProtROIVoting
from pwchem.utils import assertHandle

defROIsStr = '''1) Coordinate: {"X": 45, "Y": 65, "Z": 60}
2) Residues: {"model": 0, "chain": "A", "index": "1-4", "residues": "VLSP"}'''

defROIsStr2 = '''1) Residues: {"model": 0, "chain": "A", "index": "61-63", "residues": "KVA"}
2) Residues: {"model": 0, "chain": "A", "index": "80-84", "residues": "LSALS"}'''


class TestROIVoting(BaseTest):
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

    @classmethod
    def _runROIVoting(cls, protDefs):
        prot = cls.newProtocol(ProtROIVoting)
        for protDef in protDefs:
            prot.roisList.append(protDef)
            prot.roisList[-1].setExtended('outputStructROIs')
        cls.proj.launchProtocol(prot)
        return prot

    def test(self):
        pDef1 = self._runDefStructROIs(defROIsStr)
        pDef2 = self._runDefStructROIs(defROIsStr2)
        self._waitOutput(pDef1, 'outputStructROIs', sleepTime=10)
        self._waitOutput(pDef2, 'outputStructROIs', sleepTime=10)

        pVoting = self._runROIVoting([pDef1, pDef2])
        self._waitOutput(pVoting, 'outputStructROIs', sleepTime=10)
        assertHandle(self.assertIsNotNone, getattr(pVoting, 'outputStructROIs', None),
                     cwd=pVoting.getWorkingDir())