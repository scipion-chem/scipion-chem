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
from pwem.protocols import ProtImportPdb

# Scipion chem imports
from pwchem.protocols import *
from pwchem.tests import TestImportBase
from pwchem.utils import assertHandle

solStr = '''1, 0.05
2, 0.24
3, 0.5
4, 0.15'''

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

    self._waitOutput(protAdd1, 'outputObject', timeOut=10)
    self._waitOutput(protAdd2, 'outputSet', timeOut=10)
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
    self._waitOutput(protSASA, 'outputAtomStruct', timeOut=10)
    assertHandle(self.assertIsNotNone, getattr(protSASA, 'outputAtomStruct', None), cwd=protSASA.getWorkingDir())
    assertHandle(self.assertIsNotNone, getattr(protSASA, 'outputSequence', None), cwd=protSASA.getWorkingDir())