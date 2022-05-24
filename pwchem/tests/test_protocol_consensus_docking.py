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

from ..protocols import ProtocolConsensusDocking, ProtChemOBabelPrepareLigands, \
  ProtChemImportSmallMolecules, ProtDefinePockets

class TestConsensusDocking(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        cls._runImportPDB()
        cls._runImportSmallMols()

        cls._runPrepareLigandsOBabel()
        cls._runPrepareReceptorADT()
        cls._waitOutput(cls.protOBabel, 'outputSmallMolecules', sleepTime=5)
        cls._waitOutput(cls.protPrepareReceptor, 'outputStructure', sleepTime=5)

        cls._runPocketsSearch()
        cls._waitOutput(cls.protDefPockets, 'outputPockets', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=1,
        pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1_noHETATM.pdb'))
      cls.launchProtocol(protImportPDB)
      cls.protImportPDB = protImportPDB

    @classmethod
    def _runImportSmallMols(cls):
      cls.protImportSmallMols = cls.newProtocol(
        ProtChemImportSmallMolecules,
        filesPath=cls.dsLig.getFile('mol2'))
      cls.launchProtocol(cls.protImportSmallMols)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
      cls.protOBabel = cls.newProtocol(
        ProtChemOBabelPrepareLigands,
        inputType=0, method_charges=0,
        inputSmallMols=cls.protImportSmallMols.outputSmallMolecules,
        doConformers=True, method_conf=0, number_conf=2,
        rmsd_cutoff=0.375)

      cls.proj.launchProtocol(cls.protOBabel, wait=False)

    @classmethod
    def _runPrepareReceptorADT(cls):
      try:
          from autodock.protocols import ProtChemADTPrepareReceptor
          cls.protPrepareReceptor = cls.newProtocol(
            ProtChemADTPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            HETATM=True, rchains=True,
            chain_name='{"model": 0, "chain": "C", "residues": 93}',
            repair=3)

          cls.launchProtocol(cls.protPrepareReceptor)
      except:
          print('Autodock plugin is necesssary to run this test')

    @classmethod
    def _runPocketsSearch(cls):
      cls.protDefPockets = cls.newProtocol(
        ProtDefinePockets,
        inputAtomStruct=cls.protPrepareReceptor.outputStructure,
        inResidues='{"model": 0, "chain": "C", "index": "58-58", "residues": "H"}\n'
                   '{"model": 0, "chain": "C", "index": "101-101", "residues": "L"}')

      cls.proj.launchProtocol(cls.protDefPockets, wait=False)

    def _runVina(self, prot, pockets=None):
      if pockets == None:
        protVina = self.newProtocol(
          prot,
          wholeProt=True,
          inputAtomStruct=self.protPrepareReceptor.outputStructure,
          inputLibrary=self.protOBabel.outputSmallMolecules,
          radius=24, nPos=2,
          numberOfThreads=8)
        self.proj.launchProtocol(protVina, wait=False)

      else:
        protVina = self.newProtocol(
          prot,
          wholeProt=False,
          inputPockets=pockets,
          inputLibrary=self.protOBabel.outputSmallMolecules,
          pocketRadiusN=1.4, nPos=2,
          mergeOutput=True,
          numberOfThreads=4)
        self.proj.launchProtocol(protVina, wait=False)
      return protVina

    def _runAutoDock(self, prot, pockets=None):
        if pockets == None:
            protAutoDock = self.newProtocol(
                prot,
                wholeProt=True,
                inputAtomStruct=self.protPrepareReceptor.outputStructure,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                radius=24, gaRun=2,
                numberOfThreads=8)
            self.proj.launchProtocol(protAutoDock, wait=False)

        else:
            protAutoDock = self.newProtocol(
                prot,
                wholeProt=False,
                inputPockets=pockets,
                inputLibrary=self.protOBabel.outputSmallMolecules,
                pocketRadiusN=1.4, gaRun=2,
                mergeOutput=True,
                numberOfThreads=4)
            self.proj.launchProtocol(protAutoDock, wait=False)

        return protAutoDock

    def _runConsensusDocking(self, inputDockProts):
      protConsDocks = self.newProtocol(ProtocolConsensusDocking,
                                       action='_energy', maxmin=False)

      for i in range(len(inputDockProts)):
          protConsDocks.inputMoleculesSets.append(inputDockProts[i])
          protConsDocks.inputMoleculesSets[i].setExtended('outputSmallMolecules')

      self.launchProtocol(protConsDocks, wait=True)
      self.assertIsNotNone(protConsDocks.outputSmallMolecules,
                           "There was a problem with the consensus")

    def testConsensusDocking(self):
      inputDockProts = []
      doAD4, doVina = False, False

      from autodock.protocols import ProtChemAutodock, ProtChemVina
      protVina = self._runVina(prot=ProtChemVina, pockets=self.protDefPockets.outputPockets)
      inputDockProts.append(protVina)
      doVina = True

      protAD4 = self._runAutoDock(prot=ProtChemAutodock, pockets=self.protDefPockets.outputPockets)
      inputDockProts.append(protAD4)
      doAD4 = True

      self._waitOutput(protVina, 'outputSmallMolecules', sleepTime=5)
      self._waitOutput(protAD4, 'outputSmallMolecules', sleepTime=5)


      #TODO: alternative docking protocols (leDock, rDock)

      if doVina or doAD4:
        self._runConsensusDocking(inputDockProts)
      else:
        print('No pocket plugins found installed. Try installing P2Rank or Fpocket')

