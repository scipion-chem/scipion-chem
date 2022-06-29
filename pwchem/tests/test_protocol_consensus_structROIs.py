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

from ..protocols import ProtocolConsensusStructROIs

class TestConsensusStructROIs(BaseTest):
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

    def _runFPocketFind(self, prot):
        protFPocket = self.newProtocol(
            prot,
            inputAtomStruct=self.protImportPDB.outputPdb)

        self.launchProtocol(protFPocket)
        pocketsOut = getattr(protFPocket, 'outputStructROIs', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut
      
    def _runP2RankFind(self, prot):
        protP2Rank = self.newProtocol(
            prot,
            inputAtomStruct=self.protImportPDB.outputPdb)

        self.launchProtocol(protP2Rank)
        pocketsOut = getattr(protP2Rank, 'outputStructROIs', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut

    def _runCreateGrid(self, prot):
        protGrid = self.newProtocol(
            prot,
            inputAtomStruct=self.protImportPDB.outputPdb,
            radius=37.0,
            spacing=1.0)

        self.launchProtocol(protGrid)
        gridOut = getattr(protGrid, 'outputGrid', None)
        self.assertIsNotNone(gridOut)
        return protGrid

    def _runAutoLigandFind(self, protGrid, prot):
        protAutoLigand = self.newProtocol(
            prot,
            prevGrid=True,
            inputGrid=protGrid.outputGrid,
            nFillPoints=10)

        self.launchProtocol(protAutoLigand)
        pocketsOut = getattr(protAutoLigand, 'outputStructROIs', None)
        self.assertIsNotNone(pocketsOut)
        return pocketsOut

    def _runConsensusStructROIs(self, inputSetsOfPockets):
      prot = self.newProtocol(ProtocolConsensusStructROIs,
                              inputPocketSets= inputSetsOfPockets)

      self.launchProtocol(prot, wait=True)
      self.assertIsNotNone(prot.outputStructROIs,
                           "There was a problem with the consensus")

    def testConsensusStructROIs(self):
      inputPockets = []
      doFPocket, doAutoLigand, doP2Rank = False, False, False

      try:
        from p2rank.protocols import P2RankFindPockets
        p2rankPockets = self._runP2RankFind(prot=P2RankFindPockets)
        inputPockets.append(p2rankPockets)
        doP2Rank = True
      except:
        print('P2Rank plugin not found. Skipping its part')

      try:
        from fpocket.protocols import FpocketFindPockets
        fpocketPockets = self._runFPocketFind(prot=FpocketFindPockets)
        inputPockets.append(fpocketPockets)
        doFPocket = True
      except:
        print('FPocket plugin not found. Skipping its part')

      try:
        from autodock.protocols import ProtChemAutoLigand, Autodock_GridGeneration
        protGrid = self._runCreateGrid(prot=Autodock_GridGeneration)
        autoligandPockets = self._runAutoLigandFind(protGrid, prot=ProtChemAutoLigand)
        inputPockets.append(autoligandPockets)
        doAutoLigand = True
      except:
        print('Autodock plugin not found. Skipping AutoLigand part')

      if doFPocket or doAutoLigand or doP2Rank:
          self._runConsensusStructROIs(inputPockets)
      else:
          print('No pocket plugins found installed. Try installing P2Rank or Fpocket')


