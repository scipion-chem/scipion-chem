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

import glob
import os

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportPdb

from pwchem.protocols import (ProtocolImportMDSystem, ProtocolProlif, ProtocolTrajectoryClustering,
                              ProtocolInverseDistanceWeighting, ProtocolMDVideo)
from pwchem.protocols.MolecularDynamics.protocol_trajectory_clustering import MDSYSTEM, SETOFSTRUCTS
from pwchem.utils import assertHandle
from pwchem.tests import TestImportMDSystems

class TestProLIFanalysis(TestImportMDSystems):
    @classmethod
    def setUpClass(cls):
        cls.ds = DataSet.getDataSet('mdSystem')
        setupTestProject(cls)
        cls._runImportSystem()

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

class TestTrajClustering(TestImportMDSystems):

    @classmethod
    def _runClustering(cls, inputSetOfStructs=None, nGroup=10):
        inputArgs = {
            'inputFrom': MDSYSTEM,
            'inputMDSystem': cls.protImportMDSystem.outputSystem
        }
        if inputSetOfStructs is not None:
            inputArgs = {
                'inputFrom': SETOFSTRUCTS,
                'inputSetOfStructs': inputSetOfStructs,
                'refStructId': 1
            }

        protClust = cls.newProtocol(
            ProtocolTrajectoryClustering,
            clusterMode=1,
            nGroup=nGroup,
            **inputArgs
        )
        cls.launchProtocol(protClust)
        cls.outputSetOfAtomStructures = protClust
        return protClust

    def test_clustering(self):
        protClust = self._runClustering()
        assertHandle(self.assertIsNotNone, getattr(protClust, 'outputAtomStructs', None),
                     cwd=protClust.getWorkingDir())

        protClustFromSet = self._runClustering(inputSetOfStructs=protClust.outputAtomStructs, nGroup=3)
        assertHandle(self.assertIsNotNone, getattr(protClustFromSet, 'outputAtomStructs', None),
                     cwd=protClustFromSet.getWorkingDir())

class TestIDWreconstruct(TestTrajClustering):
    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb, inputPdbData=0,
            pdbId='3p6h')
        cls.launchProtocol(protImportPDB, wait=True)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runIDWreconstruct(cls, clust):
        protIDW = cls.newProtocol(
            ProtocolInverseDistanceWeighting,
            inputReference=cls.protImportPDB.outputPdb,
            inputEnsemble=clust
        )
        cls.launchProtocol(protIDW, wait=True)
        cls.protIDW = protIDW

    def test_idw(self):
        self._runImportSystem()
        self._runImportPDB()
        
        protClust = self._runClustering()
        self._runIDWreconstruct(protClust.outputAtomStructs)

        assertHandle(self.assertIsNotNone, getattr(self.protIDW, 'outputAtomStructs', None),
                     cwd=self.protIDW.getWorkingDir())


class TestMDVideo(TestImportMDSystems):

    @classmethod
    def _runVideo(cls, vidFormat=0, keepFrames=False):
        protVideo = cls.newProtocol(
            ProtocolMDVideo,
            inputMDSystem=cls.protImportMDSystem.outputSystem,
            vidResolution=0,        # 480p
            vidRay=False,
            maxFrames=10,
            vidFormat=vidFormat,
            keepFrames=keepFrames,
            numberOfThreads=2,
        )
        cls.launchProtocol(protVideo)
        return protVideo

    def _checkVideo(self, protVideo):
        assertHandle(self.assertIsNotNone, getattr(protVideo, 'outputVideo', None),
                     cwd=protVideo.getWorkingDir())
        videoFile = protVideo.outputVideo.getFileName()
        assertHandle(self.assertTrue, os.path.getsize(videoFile) > 0, cwd=protVideo.getWorkingDir())
        return videoFile

    def test_video_mp4(self):
        protVideo = self._runVideo(vidFormat=0, keepFrames=True)
        self._checkVideo(protVideo)

        plan = protVideo._readPlan()
        frameFiles = glob.glob(protVideo._getExtraPath('frames', 'chunk_*', 'frame_*.png'))
        assertHandle(self.assertEqual, len(frameFiles), plan['nTotal'], cwd=protVideo.getWorkingDir())

    def test_video_gif(self):
        protVideo = self._runVideo(vidFormat=1, keepFrames=False)
        self._checkVideo(protVideo)
