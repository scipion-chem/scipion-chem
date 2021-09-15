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
# ***************************************************************************

import time

from pyworkflow.tests import BaseTest, DataSet
import pyworkflow.tests as tests
from pwem.protocols import ProtImportPdb
from pyworkflow.object import Pointer

from ..protocols import ProtocolConsensusPockets
try:
  doP2Rank = True
  from p2rank.protocols import P2RankFindPockets
except:
  doP2Rank = False
  print('P2Rank plugin not found. Skipping its part')

try:
  doFPocket = True
  from fpocket.protocols import FpocketFindPockets
except:
  doFPocket = False
  print('FPocket plugin not found. Skipping its part')

try:
  doAutoLigand = True
  from autodock.protocols import ProtChemAutoLigand, Autodock_GridGeneration
except:
  doAutoLigand = False
  print('Autodock plugin not found. Skipping AutoLigand part')


class TestConsensusPockets(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls._runImportPDB()

    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=1,
        pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
      cls.launchProtocol(protImportPDB)
      cls.protImportPDB = protImportPDB

    def _runFPocketFind(self):
        protFPocket = self.newProtocol(
            FpocketFindPockets,
            inputAtomStruct=self.protImportPDB.outputPdb)

        self.launchProtocol(protFPocket)
        pocketsOut = getattr(protFPocket, 'outputPockets', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut
      
    def _runP2RankFind(self):
        protP2Rank = self.newProtocol(
            P2RankFindPockets,
            inputAtomStruct=self.protImportPDB.outputPdb)

        self.launchProtocol(protP2Rank)
        pocketsOut = getattr(protP2Rank, 'outputPockets', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut

    def _runCreateGrid(self):
        protGrid = self.newProtocol(
            Autodock_GridGeneration,
            inputAtomStruct=self.protImportPDB.outputPdb,
            radius=37.0,
            spacing=1.0)

        self.launchProtocol(protGrid)
        gridOut = getattr(protGrid, 'outputGrid', None)
        self.assertIsNotNone(gridOut)
        return protGrid

    def _runAutoLigandFind(self, protGrid):
        protAutoLigand = self.newProtocol(
            ProtChemAutoLigand,
            inputGrid=protGrid.outputGrid,
            nFillPoints=10)

        self.launchProtocol(protAutoLigand)
        pocketsOut = getattr(protAutoLigand, 'outputPockets', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut

    def _runConsensusPockets(self, inputSetsOfPockets):
      prot = self.newProtocol(ProtocolConsensusPockets,
                              inputPocketSets= inputSetsOfPockets)

      self.launchProtocol(prot, wait=True)
      self.assertIsNotNone(prot.outputPocketsAll,
                           "There was a problem with the consensus")

    def testConsensusPockets(self):
      inputPockets = []
      if doAutoLigand:
        protGrid = self._runCreateGrid()
        autoligandPockets = self._runAutoLigandFind(protGrid)
        inputPockets.append(autoligandPockets)

      if doFPocket:
        fpocketPockets = self._runFPocketFind()
        inputPockets.append(fpocketPockets)

      if doP2Rank:
        p2rankPockets = self._runP2RankFind()
        inputPockets.append(p2rankPockets)

      self._runConsensusPockets(inputPockets)


