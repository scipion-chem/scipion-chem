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

from pwem.protocols import ProtImportSequence, ProtImportPdb
from ..protocols import ProtChemImportSetOfSequences, ProtChemPairWiseAlignment

CLUSTALO, MUSCLE, MAFFT = 0, 1, 2

class TestPairAlignSequences(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')

    @classmethod
    def _runImportSeq(cls):
      protImportSeq = cls.newProtocol(
        ProtImportSequence,
        inputSequence=0, inputProteinSequence=3,
        uniProtSequence='P69905', inputSequenceName='hemo')
      cls.launchProtocol(protImportSeq)
      return protImportSeq

    @classmethod
    def _runImportPDB(cls):
      protImportPDB = cls.newProtocol(
        ProtImportPdb,
        inputPdbData=0, pdbId='5ni1')
      cls.launchProtocol(protImportPDB)
      return protImportPDB

    def _runAlignment(self, protSeq, protAS):
        protAlignSeqs = self.newProtocol(
          ProtChemPairWiseAlignment,
          condAtomStruct1=False, inputSequence1=protSeq.outputSequence,
          condAtomStruct2=True, inputAtomStruct2=protAS.outputPdb,
          chain_name2='{"model": 0, "chain": "A", "residues": 141}')
        self.launchProtocol(protAlignSeqs)
        return protAlignSeqs

    def testPairwiseAlignment(self):
        protSeq = self._runImportSeq()
        protAS = self._runImportPDB()

        protAlign = self._runAlignment(protSeq, protAS)
        self.assertIsNotNone(protAlign.outputAlignment)
