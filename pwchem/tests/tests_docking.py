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

from pwchem.protocols import *
from pwchem.tests import TestDefineStructROIs

defRoiStr = '''1) Coordinate: {"X": 51, "Y": 76, "Z": 61}'''

defRoiStrLig = '''1) Ext-Ligand: {"pointerIdx": "0", "ligName": "SmallMolecule (g1_4erf_0R3-1_1 molecule)"}'''

wSteps = "{'rfSpr': 0, 'depthProt': 5, 'depthLig': 1, 'fingerSize': 65536, 'isReference': False, 'scoreChoice': 'Vina', 'scoreVersionRF': '1', 'scoreVersionPLEC': 'linear', 'trainData': '2016'}\n" \
         "{'rfSpr': 0, 'depthProt': 5, 'depthLig': 1, 'fingerSize': 65536, 'isReference': False, 'scoreChoice': 'RFScore', 'scoreVersionRF': '1', 'scoreVersionPLEC': 'linear', 'trainData': '2016'}"

prepRec = '{"model": 0, "chain": "C", "residues": 93}'

class TestExtractLigand(BaseTest):
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)
    cls.ds = DataSet.getDataSet('model_building_tutorial')
    cls._runImportPDB()

  @classmethod
  def _runImportPDB(cls):
    protImportPDB = cls.newProtocol(
      ProtImportPdb,
      inputPdbData=0, pdbId='4erf')
    cls.launchProtocol(protImportPDB)
    cls.protImportPDB = protImportPDB

  @classmethod
  def _runExtractLigand(cls, inputProt):
      protExtLig = cls.newProtocol(
        ProtExtractLigands,
        cleanPDB=True, rchains=True, chain_name='{"model": 0, "chain": "C", "residues": 93}')

      protExtLig.inputStructure.set(inputProt)
      protExtLig.inputStructure.setExtended('outputPdb')

      cls.proj.launchProtocol(protExtLig, wait=True)
      cls.protExtLig = protExtLig
      return protExtLig

  def test(self):
    protExtract = self._runExtractLigand(self.protImportPDB)
    self._waitOutput(protExtract, 'outputSmallMolecules', sleepTime=5)
    self.assertIsNotNone(protExtract.outputSmallMolecules,
                         "There was a problem with the ligand extraction")


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
        chain_name=prepRec,
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


class TestRMSDDocking(TestScoreDocking, TestExtractLigand):
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)
    cls.ds = DataSet.getDataSet('model_building_tutorial')
    cls._runImportPDB()

  @classmethod
  def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=0, pdbId='4erf')
      cls.launchProtocol(protImportPDB)
      cls.protImportPDB = protImportPDB

  @classmethod
  def _runPrepareTarget(cls, inProt):
    protPrepRec = cls.newProtocol(
      ProtChemPrepareReceptor,
      HETATM=False, rchains=True, chain_name=prepRec)

    protPrepRec.inputAtomStruct.set(inProt)
    protPrepRec.inputAtomStruct.setExtended('outputPdb')

    cls.proj.launchProtocol(protPrepRec, wait=False)
    return protPrepRec

  @classmethod
  def _runPrepareLigandsOBabel(cls, inProt):
    protOBabel = cls.newProtocol(
      ProtChemOBabelPrepareLigands,
      inputType=0, method_charges=0, doConformers=False)

    protOBabel.inputSmallMols.set(inProt)
    protOBabel.inputSmallMols.setExtended('outputSmallMolecules')

    cls.proj.launchProtocol(protOBabel, wait=False)
    return protOBabel

  @classmethod
  def _runDefStructROIs(cls, inProtAS, inProtLig, defStr):
    protDef = cls.newProtocol(
      ProtDefineStructROIs,
      inROIs=defStr, surfaceCoords=False)

    protDef.inputAtomStruct.set(inProtAS)
    protDef.inputAtomStruct.setExtended('outputStructure')

    protDef.inputPointers.append(inProtLig)
    protDef.inputPointers[0].setExtended('outputSmallMolecules')

    cls.proj.launchProtocol(protDef, wait=False)
    return protDef

  def _runRMSDDocking(self, inputMolsProt, inputProt, mode=0):
    pRMSD = self.newProtocol(ProtocolRMSDDocking,
                             refOrigin=mode)

    pRMSD.inputSmallMolecules.set(inputMolsProt)
    pRMSD.inputSmallMolecules.setExtended('outputSmallMolecules')

    if mode == 0:
        pRMSD.refAtomStruct.set(inputProt)
        pRMSD.refAtomStruct.setExtended('outputStructure')

        pRMSD.refLigName.set('0R3')
    else:
        pRMSD.refSmallMolecules.set(inputProt)
        pRMSD.refSmallMolecules.setExtended('outputSmallMolecules')

        molName = inputProt.outputSmallMolecules.getFirstItem().__str__()
        pRMSD.refMolName.set(molName)

    self.proj.launchProtocol(pRMSD, wait=False)
    return pRMSD

  def test(self):
    protPrep = self._runPrepareTarget(self.protImportPDB)
    protExtLig = self._runExtractLigand(self.protImportPDB)

    self._waitOutput(protPrep, 'outputStructure', sleepTime=5)
    self._waitOutput(protExtLig, 'outputSmallMolecules', sleepTime=5)

    protPockets = self._runDefStructROIs(protPrep, protExtLig, defRoiStrLig)
    self.protOBabel = self._runPrepareLigandsOBabel(protExtLig)

    self._waitOutput(protPockets, 'outputStructROIs', sleepTime=5)
    self._waitOutput(self.protOBabel, 'outputSmallMolecules', sleepTime=5)

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
      pRMSDs = []
      pInputs = [protPrep, protExtLig]
      for i in range(2):
          pRMSDs.append(self._runRMSDDocking(inputDockProts[0], pInputs[i], mode=i))

      for p in pRMSDs:
          self._waitOutput(p, 'outputSmallMolecules', sleepTime=5)
          self.assertIsNotNone(p.outputSmallMolecules,
                           "There was a problem with the consensus")
    else:
      print('No docking plugins found installed. Try installing AutoDock')
