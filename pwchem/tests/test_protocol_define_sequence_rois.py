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

from pyworkflow.tests import BaseTest
import pyworkflow.tests as tests
from pwem.objects import Sequence
from pwem.protocols import ProtImportSequence

from ..protocols import ProtDefineSeqROI, ProtChemImportVariants
from pwchem.objects import SequenceVariants

class TestDefineROIs(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)

    def _runImportVariants(self):
      protImportVariants = self.newProtocol(
        ProtChemImportVariants,
        fromID=True,
        inputUniProtKB='P69905')
      self.launchProtocol(protImportVariants)
      return protImportVariants

    def _runImportSeq(self):
      protImportSeq = self.newProtocol(
        ProtImportSequence,
        inputSequence=0, inputProteinSequence=3,
        uniProtSequence='P69905', inputSequenceName='hemo')
      self.launchProtocol(protImportSeq)
      return protImportSeq

    def _runDefineSeqROIs(self, inSeq):
        if type(inSeq) == Sequence:
            protAlignSeqs = self.newProtocol(
              ProtDefineSeqROI,
              inputSequence=inSeq,
              inROIs='1) Residues: {"index": "3-5", "residues": "LSP", "desc": "None"}\n'
                     '2) Residues: {"index": "48-51", "residues": "DLSH", "desc": "None"}\n')

        elif type(inSeq) == SequenceVariants:
          protAlignSeqs = self.newProtocol(
            ProtDefineSeqROI,
            inputSequenceVariants=inSeq, chooseInput=1,
            inROIs = '1) Residues: {"index": "52-55", "residues": "GSAQ", "desc": "None"}\n'
                     '2) Variant: Aichi\n'
                     '3) Mutations: H46Q\n')

        self.launchProtocol(protAlignSeqs)
        return protAlignSeqs

    def testDefineROIsFromSeq(self):
        protSeq = self._runImportSeq()
        protROIs1 = self._runDefineSeqROIs(protSeq.outputSequence)
        self.assertIsNotNone(protROIs1.outputROIs)

    def testDefineROIsFromVar(self):
        protVar = self._runImportVariants()
        protROIs2 = self._runDefineSeqROIs(protVar.outputVariants)
        self.assertIsNotNone(protROIs2.outputROIs)

