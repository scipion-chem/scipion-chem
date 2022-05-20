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

from ..protocols import ProtChemImportSetOfSequences, ProtChemMultipleSequenceAlignment

CLUSTALO, MUSCLE, MAFFT = 0, 1, 2

class TestAlignSequences(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls._runImportSeqs()

    @classmethod
    def _runImportSeqs(cls):
      protImportSeqs = cls.newProtocol(
        ProtChemImportSetOfSequences,
        multiple=False,
        filePath=cls.ds.getFile('Sequences/Several_sequences.fasta'))
      cls.launchProtocol(protImportSeqs)
      cls.protImportSeqs = protImportSeqs

    def _runAlignment(self, program=0):
        protAlignSeqs = self.newProtocol(
          ProtChemMultipleSequenceAlignment,
          setOfSequences=self.protImportSeqs.outputSequences,
          programList=program)
        self.proj.launchProtocol(protAlignSeqs, wait=False)
        return protAlignSeqs

    def testMultipleAlignment(self):
        protCLUSTALO = self._runAlignment(CLUSTALO)
        protMUSCLE = self._runAlignment(MUSCLE)
        protMAFFT = self._runAlignment(MAFFT)

        self._waitOutput(protCLUSTALO, 'outputSequences', sleepTime=5)
        self._waitOutput(protMUSCLE, 'outputSequences', sleepTime=5)
        self._waitOutput(protMAFFT, 'outputSequences', sleepTime=5)

        self.assertIsNotNone(protCLUSTALO.outputSequences)
        self.assertIsNotNone(protMUSCLE.outputSequences)
        self.assertIsNotNone(protMAFFT.outputSequences)

