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

from ..protocols import ProtChemGenerateVariants, ProtChemImportVariants, \
  ProtExtractSeqsROI

class TestExtractROIs(BaseTest):
    @classmethod
    def setUpClass(cls):
        tests.setupTestProject(cls)
        cls._runImportVariants()

    @classmethod
    def _runImportVariants(cls):
      protImportVariants = cls.newProtocol(
        ProtChemImportVariants,
        fromID=True,
        inputUniProtKB='P69905')
      cls.launchProtocol(protImportVariants)
      cls.protImportVariants = protImportVariants

    def _runGenerateVariants(self):
        protGenerateVariants = self.newProtocol(
          ProtChemGenerateVariants,
          inputSequenceVariants=self.protImportVariants.outputVariants,
          toMutateList='1) Variant: Arbor\n'
                       '2) Variant: Kawachi\n'
                       '3) Mutations: A6D\n')
        self.launchProtocol(protGenerateVariants)
        return protGenerateVariants

    def _runExtractROIs(self, inputSeqs):
        protExtractROIs = self.newProtocol(
          ProtExtractSeqsROI,
          inputSequences=inputSeqs, direction=1, thres=0.2)
        self.launchProtocol(protExtractROIs)
        return protExtractROIs

    def testGenerateVariants(self):
        protVars = self._runGenerateVariants()
        protROIs = self._runExtractROIs(protVars.outputSequences)
        self.assertIsNotNone(protROIs.outputROIs)
