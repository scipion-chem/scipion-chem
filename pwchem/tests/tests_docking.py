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

from pwchem.protocols import ProtocolConsensusDocking, ProtChemOBabelPrepareLigands, \
  ProtChemImportSmallMolecules, ProtDefineStructROIs, ProtocolScoreDocking, ProtocolRMSDDocking
from pwchem.tests import TestDefineStructROIs

defRoiStr = '''1) Coordinate: {"X": 51, "Y": 76, "Z": 61}'''

wSteps = "{'rfSpr': 0, 'depthProt': 5, 'depthLig': 1, 'fingerSize': 65536, 'isReference': False, 'scoreChoice': 'Vina', 'scoreVersionRF': '1', 'scoreVersionPLEC': 'linear', 'trainData': '2016'}\n" \
         "{'rfSpr': 0, 'depthProt': 5, 'depthLig': 1, 'fingerSize': 65536, 'isReference': False, 'scoreChoice': 'RFScore', 'scoreVersionRF': '1', 'scoreVersionPLEC': 'linear', 'trainData': '2016'}"


class TestScoreDocking(TestDefineStructROIs):
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
      doConformers=False)

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
  def _runDefStructROIs(cls, defStr):
    protDef = cls.newProtocol(
      ProtDefineStructROIs,
      inROIs=defStr)

    protDef.inputAtomStruct.set(cls.protPrepareReceptor)
    protDef.inputAtomStruct.setExtended('outputStructure')

    cls.proj.launchProtocol(protDef, wait=False)
    return protDef

  def _runVina(self, prot, pocketsProt=None):
    if pocketsProt == None:
      protVina = self.newProtocol(
        prot,
        fromReceptor=0,
        radius=24, nRuns=10,
        numberOfThreads=4)

      protVina.inputAtomStruct.set(self.protPrepareReceptor)
      protVina.inputAtomStruct.setExtended('outputStructure')

    else:
      protVina = self.newProtocol(
        prot,
        fromReceptor=1,
        pocketRadiusN=3, nRuns=3,
        numberOfThreads=4)

      protVina.inputStructROIs.set(pocketsProt)
      protVina.inputStructROIs.setExtended('outputStructROIs')

    protVina.inputSmallMolecules.set(self.protOBabel)
    protVina.inputSmallMolecules.setExtended('outputSmallMolecules')

    self.proj.launchProtocol(protVina, wait=False)
    return protVina

  def _runLeDock(self, prot, pocketsProt=None):
    if pocketsProt == None:
      protLeDock = self.newProtocol(
        prot,
        wholeProt=True,
        radius=24, nRuns=10,
        numberOfThreads=4)

      protLeDock.inputAtomStruct.set(self.protPrepareReceptor)
      protLeDock.inputAtomStruct.setExtended('outputStructure')

    else:
      protLeDock = self.newProtocol(
        prot,
        wholeProt=False,
        pocketRadiusN=3, nRuns=3,
        numberOfThreads=4)

      protLeDock.inputStructROIs.set(pocketsProt)
      protLeDock.inputStructROIs.setExtended('outputStructROIs')

    protLeDock.inputSmallMolecules.set(self.protOBabel)
    protLeDock.inputSmallMolecules.setExtended('outputSmallMolecules')

    self.proj.launchProtocol(protLeDock, wait=False)
    return protLeDock

  def _runScoreDocking(self, inputDockProts):
      protScoreDocks = self.newProtocol(ProtocolScoreDocking,
                                        numberOfThreads=4,
                                        summarySteps='1) Score: Vina\n'
                                                     '2) Score: RFScore, version 1. PDBbind 2016',
                                        workFlowSteps=wSteps)

      for i in range(len(inputDockProts)):
        protScoreDocks.inputMoleculesSets.append(inputDockProts[i])
        protScoreDocks.inputMoleculesSets[i].setExtended('outputSmallMolecules')

      self.proj.launchProtocol(protScoreDocks, wait=True)
      return protScoreDocks

  def test(self):
    protPockets = self._runDefStructROIs(defRoiStr)
    self._waitOutput(protPockets, 'outputStructROIs', sleepTime=5)

    inputDockProts = []
    doVina, doLe = False, False

    try:
      from autodock.protocols import ProtChemVinaDocking
      doVina = True
    except:
      pass
    try:
      from lephar.protocols import ProtChemLeDock
      doLe = True
    except:
      pass

    if doVina:
      protVina1 = self._runVina(prot=ProtChemVinaDocking)
      inputDockProts.append(protVina1)
      if not doLe:
        protVina2 = self._runVina(prot=ProtChemVinaDocking, pocketsProt=protPockets)
        inputDockProts.append(protVina2)

    if doLe:
      protLe1 = self._runLeDock(prot=ProtChemLeDock, pocketsProt=protPockets)
      inputDockProts.append(protLe1)
      if not doVina:
        protLe2 = self._runLeDock(prot=ProtChemLeDock)
        inputDockProts.append(protLe2)

    for p in inputDockProts:
      self._waitOutput(p, 'outputSmallMolecules', sleepTime=5)

    if len(inputDockProts) >= 1:
        pScore = self._runScoreDocking(inputDockProts)
        self.assertIsNotNone(pScore.outputSmallMolecules,
                             "There was a problem with the consensus")
    else:
        print('No docking plugins found installed. Try installing AutoDock or LePhar')


class TestConsensusDocking(TestScoreDocking):

    def _runConsensusDocking(self, inputDockProts):
      protConsDocks = self.newProtocol(ProtocolConsensusDocking, numberOfThreads=1,
                                       repAttr='_energy', maxmin=False)

      for i in range(len(inputDockProts)):
          protConsDocks.inputMoleculesSets.append(inputDockProts[i])
          protConsDocks.inputMoleculesSets[i].setExtended('outputSmallMolecules')

      self.launchProtocol(protConsDocks, wait=True)
      return protConsDocks

    def test(self):
      protPockets = self._runDefStructROIs(defRoiStr)
      self._waitOutput(protPockets, 'outputStructROIs', sleepTime=5)

      inputDockProts = []
      doVina, doLe = False, False

      try:
          from autodock.protocols import ProtChemVinaDocking
          doVina = True
      except:
          pass
      try:
          from lephar.protocols import ProtChemLeDock
          doLe = True
      except:
          pass

      if doVina:
          protVina1 = self._runVina(prot=ProtChemVinaDocking)
          inputDockProts.append(protVina1)
          if not doLe:
              protVina2 = self._runVina(prot=ProtChemVinaDocking, pocketsProt=protPockets)
              inputDockProts.append(protVina2)

      if doLe:
        protLe1 = self._runLeDock(prot=ProtChemLeDock, pocketsProt=protPockets)
        inputDockProts.append(protLe1)
        if not doVina:
            protLe2 = self._runLeDock(prot=ProtChemLeDock)
            inputDockProts.append(protLe2)

      for p in inputDockProts:
          self._waitOutput(p, 'outputSmallMolecules', sleepTime=5)

      if len(inputDockProts) >= 2:
          pCons = self._runConsensusDocking(inputDockProts)
          self.assertIsNotNone(pCons.outputSmallMolecules,
                               "There was a problem with the consensus")
      else:
          print('No docking plugins found installed. Try installing AutoDock')


class TestRMSDDocking(TestScoreDocking):

  def _runRMSDDocking(self, inputMolsProt, inputAS=None):
    pRMSD = self.newProtocol(ProtocolRMSDDocking,
                             refOrigin=1)

    molName = inputMolsProt.outputSmallMolecules.getFirstItem().__str__()
    pRMSD.refMolName.set(molName)

    pRMSD.inputSmallMolecules.set(inputMolsProt)
    pRMSD.inputSmallMolecules.setExtended('outputSmallMolecules')
    pRMSD.refSmallMolecules.set(inputMolsProt)
    pRMSD.refSmallMolecules.setExtended('outputSmallMolecules')

    self.launchProtocol(pRMSD, wait=True)
    return pRMSD

  def test(self):
    protPockets = self._runDefStructROIs(defRoiStr)
    self._waitOutput(protPockets, 'outputStructROIs', sleepTime=5)

    inputDockProts = []
    doVina, doLe = False, False

    try:
      from autodock.protocols import ProtChemVinaDocking
      doVina = True
    except:
      pass
    try:
      from lephar.protocols import ProtChemLeDock
      doLe = True
    except:
      pass

    if doVina:
      protVina1 = self._runVina(prot=ProtChemVinaDocking, pocketsProt=protPockets)
      inputDockProts.append(protVina1)

    if doLe and not doVina:
      protLe1 = self._runLeDock(prot=ProtChemLeDock, pocketsProt=protPockets)
      inputDockProts.append(protLe1)

    for p in inputDockProts:
      self._waitOutput(p, 'outputSmallMolecules', sleepTime=5)

    if len(inputDockProts) >= 1:
      pRMSD = self._runRMSDDocking(inputDockProts[0])
      self.assertIsNotNone(pRMSD.outputSmallMolecules,
                           "There was a problem with the consensus")
    else:
      print('No docking plugins found installed. Try installing AutoDock')

