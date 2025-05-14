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

# Scipion em imports
from pyworkflow.tests import BaseTest, DataSet
import pyworkflow.tests as tests
from pwem.protocols import ProtImportPdb, ProtImportSequence

# Scipion chem imports
from pwchem.protocols import *
from pwchem.tests import TestImportBase
from pwchem.utils import assertHandle

from pwchem.tests.tests_sequences import defSeqROIsSeq

solStr = '''1, 0.05
2, 0.24
3, 0.5
4, 0.15'''

solStr2 = '''1, 0.12
2, 0.2
3, 0.54
4, 0.7'''

INATTR = '''{"Set Idx": "0", "ID": "molName", "Values": "Solubility", "Higher_is_better": "True"}
{"Set Idx": "1", "ID": "molName", "Values": "Solubility", "Higher_is_better": "True"}'''

class TestAddAttribute(TestImportBase):
  @classmethod
  def _runAddAttribute(cls, inputProt, mode=0, inputFile=None):
    protAdd = cls.newProtocol(
      ProtAddAttribute,
      fromInput=mode, atKey='Solubility')

    if mode == 0:
      protAdd.atValue.set("0.5")
      protAdd.inputObject.set(inputProt)
      protAdd.inputObject.setExtended('outputSmallMolecules')
    else:
      protAdd.mapKey.set('_objId')
      protAdd.inputFile.set(inputFile)
      protAdd.inputSet.set(inputProt)
      protAdd.inputSet.setExtended('outputSmallMolecules')

    cls.proj.launchProtocol(protAdd)
    return protAdd

  def test(self):
    protAdd1 = self._runAddAttribute(self.protImportSmallMols, mode=0)

    inFile = self.proj.getTmpPath('attributeMapping.csv')
    with open(inFile, 'w') as f:
      f.write(solStr)
    protAdd2 = self._runAddAttribute(self.protImportSmallMols, mode=1, inputFile=inFile)

    self._waitOutput(protAdd1, 'outputObject', sleepTime=10)
    self._waitOutput(protAdd2, 'outputSet', sleepTime=10)
    assertHandle(self.assertIsNotNone, getattr(protAdd1, 'outputObject', None), cwd=protAdd1.getWorkingDir())
    assertHandle(self.assertIsNotNone, getattr(protAdd2, 'outputSet', None), cwd=protAdd2.getWorkingDir())


class TestCalculateSASA(BaseTest):
  @classmethod
  def setUpClass(cls):
    tests.setupTestProject(cls)
    cls.ds = DataSet.getDataSet('model_building_tutorial')
    cls._runImportPDB()

  @classmethod
  def _runImportPDB(cls):
    protImportPDB = cls.newProtocol(
      ProtImportPdb,
      inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/1aoi.cif'))
    cls.launchProtocol(protImportPDB)
    cls.protImportPDB = protImportPDB

  @classmethod
  def _runCalculateSASA(cls, inputProt):
      protSASA = cls.newProtocol(
        ProtCalculateSASA,
        extractSequence=True, chain_name='{"model": 0, "chain": "A", "residues": 98}')

      protSASA.inputAtomStruct.set(inputProt)
      protSASA.inputAtomStruct.setExtended('outputPdb')

      cls.proj.launchProtocol(protSASA, wait=True)
      return protSASA

  def test(self):
    protSASA = self._runCalculateSASA(self.protImportPDB)
    self._waitOutput(protSASA, 'outputAtomStruct', sleepTime=10)
    assertHandle(self.assertIsNotNone, getattr(protSASA, 'outputAtomStruct', None), cwd=protSASA.getWorkingDir())
    assertHandle(self.assertIsNotNone, getattr(protSASA, 'outputSequence', None), cwd=protSASA.getWorkingDir())

class TestMapAttributes(TestCalculateSASA):
  @classmethod
  def _runImportSequence(cls):
    args = {
      'inputSequenceName': 'User_seq',
      'inputProteinSequence': ProtImportSequence.IMPORT_FROM_FILES,
      'fileSequence': cls.ds.getFile('Sequences/1aoi_A_mutated.fasta')
    }
    cls.protImportSeq = cls.newProtocol(ProtImportSequence, **args)
    cls.launchProtocol(cls.protImportSeq)

  @classmethod
  def _runDefSeqROIs(cls, inProt):

    protDefSeqROIs = cls.newProtocol(
      ProtDefineSeqROI,
      chooseInput=0, inROIs=defSeqROIsSeq
    )
    protDefSeqROIs.inputSequence.set(inProt)
    protDefSeqROIs.inputSequence.setExtended('outputSequence')

    cls.proj.launchProtocol(protDefSeqROIs, wait=True)
    return protDefSeqROIs

  def _runMapROIsAttr(self, protROI, protAttr, mode=0):
      protMap = self.newProtocol(
        ProtMapAttributeToSeqROIs,
        inputFrom=mode, chain_name='{"model": 0, "chain": "A", "residues": 98}', attrName='SASA')

      protMap.inputSequenceROIs.set(protROI)
      protMap.inputSequenceROIs.setExtended('outputROIs')

      if mode == 0:
        protMap.inputAtomStruct.set(protAttr)
        protMap.inputAtomStruct.setExtended('outputAtomStruct')
      else:
        protMap.inputSequence.set(protAttr)
        protMap.inputSequence.setExtended('outputSequence')

      self.proj.launchProtocol(protMap, wait=False)
      return protMap


  def test(self):
    self._runImportSequence()
    protSASA = self._runCalculateSASA(self.protImportPDB)
    protROIs = self._runDefSeqROIs(self.protImportSeq)

    protMap0 = self._runMapROIsAttr(protROIs, protSASA, mode=0)
    protMap1 = self._runMapROIsAttr(protROIs, protSASA, mode=1)

    self._waitOutput(protMap0, 'outputROIs', sleepTime=10)
    self._waitOutput(protMap1, 'outputROIs', sleepTime=10)
    assertHandle(self.assertIsNotNone, getattr(protMap0, 'outputROIs', None), cwd=protMap0.getWorkingDir())
    assertHandle(self.assertIsNotNone, getattr(protMap1, 'outputROIs', None), cwd=protMap1.getWorkingDir())
    
class TestRanxFusion(TestAddAttribute):
  @classmethod
  def _runRanxFusion(cls, inputProts):
    protRanx = cls.newProtocol(
      ProtocolRANXFuse, inAttrs=INATTR)

    for i in range(len(inputProts)):
      protRanx.inputSets.append(inputProts[i])
      protRanx.inputSets[i].setExtended('outputSet')


    cls.proj.launchProtocol(protRanx)
    return protRanx

  def test(self):
    inFile1, inFile2 = self.proj.getTmpPath('attributeMapping.csv'), self.proj.getTmpPath('attributeMapping2.csv')
    with open(inFile1, 'w') as f:
      f.write(solStr)
    protAdd1 = self._runAddAttribute(self.protImportSmallMols, mode=1, inputFile=inFile1)

    with open(inFile2, 'w') as f:
      f.write(solStr2)
    protAdd2 = self._runAddAttribute(self.protImportSmallMols, mode=1, inputFile=inFile2)

    self._waitOutput(protAdd1, 'outputSet', sleepTime=10)
    self._waitOutput(protAdd2, 'outputSet', sleepTime=10)
    assertHandle(self.assertIsNotNone, getattr(protAdd1, 'outputSet', None), cwd=protAdd1.getWorkingDir())
    assertHandle(self.assertIsNotNone, getattr(protAdd2, 'outputSet', None), cwd=protAdd2.getWorkingDir())

    protRanx = self._runRanxFusion([protAdd1, protAdd2])
    self._waitOutput(protRanx, 'outputSet', sleepTime=10)
    assertHandle(self.assertIsNotNone, getattr(protRanx, 'outputSet', None), cwd=protRanx.getWorkingDir())

