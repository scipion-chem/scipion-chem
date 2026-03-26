# **************************************************************************
# *
# *
# * Authors:    Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import *
from pwchem.constants import *
from pwchem.tests import MDSYSTEM, DataSetMDSystem
from pwchem.utils import assertHandle

class TestProLIFanalysis(BaseTest):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet(MDSYSTEM)
        setupTestProject(cls)
        cls._runImportSystem()

    @classmethod
    def _runImportSystem(cls):
        protImportMDSystem = cls.newProtocol(
            ProtocolImportMDSystem,
            inputCoords=cls.ds.getFile(DataSetMDSystem.amberSystemFile.value),
            inputCrd=cls.ds.getFile(DataSetMDSystem.amberCrdFile.value),
            inputTopology=cls.ds.getFile(DataSetMDSystem.amberTopoFile.value),
            addTraj=True,
            inputTrajectory=cls.ds.getFile(DataSetMDSystem.amberTrjFile.value),
            hasLig=True,
            ligID='LIG')
        cls.launchProtocol(protImportMDSystem)
        cls.protImportMDSystem = protImportMDSystem
        return protImportMDSystem

    @classmethod
    def _runProLIF(cls):
        protProLIF = cls.newProtocol(
            ProtocolProlif,
            inputMDSystem=cls.protImportMDSystem.outputSystem,
            frameNum=2
        )
        cls.launchProtocol(protProLIF)
        cls.prolifSystem = protProLIF
        return protProLIF

    def test_prolif(self):
        self._runProLIF()
        assertHandle(self.assertIsNotNone,getattr(self.prolifSystem, 'outputSystem', None))
        assertHandle(
            self.assertIsNotNone,
            hasattr(self.prolifSystem, '_prolifFile'),
            cwd=self.prolifSystem.getWorkingDir()
        )