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

from ..protocols import ProtDefinePockets

class TestDefinePockets(BaseTest):
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

    @classmethod
    def _runPocketsSearch(cls):
      cls.protDefPockets = cls.newProtocol(
        ProtDefinePockets,
        inputAtomStruct=cls.protImportPDB.outputPdb,
        inResidues='{"model": 0, "chain": "A", "index": "58-58", "residues": "H"}\n'
                   '{"model": 0, "chain": "C", "index": "101-101", "residues": "L"}')

      cls.proj.launchProtocol(cls.protDefPockets, wait=False)

    def testConsensusDocking(self):
        self._runPocketsSearch()

