# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
import os.path
import time
from pathlib import Path

from fpocket.protocols import FpocketFindPockets
from pwem.protocols import ProtImportPdb
from pyworkflow.tests import BaseTest, DataSet

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules
from pwem.protocols.protocol_set_filter import ProtSetFilter
import pyworkflow.tests as tests

from ..objects import SetOfStructROIs
from ..protocols import ProtocolSCORCH2


class TestSCORCH2(BaseTest):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        tests.setupTestProject(cls)

        cls._runImportPDB()
        cls._runImportSmallMols()


        cls._runFPocketFind()
        cls._runFilterSet()
        cls._runDocking()


    @classmethod
    def _runFilterSet(cls):
        #cls.filteredSet = cls.newProtocol(
        #    ProtSetFilter, inputSet=cls.protFPocket,
        #    operation='ranking',
        #    threshold=0.2,
        #    rankingField='_score')
        cls.filteredSet = cls.newProtocol(ProtSetFilter)
        cls.filteredSet.inputSet.set(cls.protFPocket)
        cls.filteredSet.inputSet.setExtended('outputStructROIs')
        cls.filteredSet.operation.set(cls.filteredSet.CHOICE_RANKED)
        cls.filteredSet.threshold.set(2)
        cls.filteredSet.rankingField.set('_score')
        cls.launchProtocol(cls.filteredSet)

    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('pdb'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=True)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(cls.protImportPDB, wait=True)

    @classmethod
    def _runFPocketFind(cls):
        cls.protFPocket = cls.newProtocol(
            FpocketFindPockets,
            inputAtomStruct=cls.protImportPDB.outputPdb)

        cls.launchProtocol(cls.protFPocket)

    @classmethod
    def _runDocking(cls):
        try:
            from lephar.protocols import ProtChemLeDock
            doLe = True
        except:
            pass
        if doLe:
            cls.protLeDock = cls.newProtocol(
                ProtChemLeDock,
                wholeProt=False,
                pocketRadiusN=3, nRuns=3,
                numberOfThreads=4)

            cls.protLeDock.inputStructROIs.set(cls.filteredSet)
            cls.protLeDock.inputStructROIs.setExtended('outputStructROIs')

            cls.protLeDock.inputSmallMolecules.set(cls.protImportSmallMols)
            cls.protLeDock.inputSmallMolecules.setExtended('outputSmallMolecules')

            cls.launchProtocol(cls.protLeDock)
        try:
            from autodock.protocols import ProtChemVinaDocking
            doVina = True
        except:
            pass
        if doVina and not doLe:
            cls.protVina = cls.newProtocol(
                ProtChemVinaDocking,
                fromReceptor=0,
                radius=24, nRuns=10,
                numberOfThreads=4)

            cls.protVina.inputAtomStruct.set(cls.protFPocket)
            cls.protVina.inputAtomStruct.setExtended('outputStructROIs')

            cls.launchProtocol(cls.protVina)

        if not doLe and not doVina:
            exit()



    def _runSCORCH2(self):
        protSCORCH2 = self.newProtocol(ProtocolSCORCH2)

        protSCORCH2.useFeatures.set('False')
        #protSCORCH2.inputPDBproteinFile.set(self.targetProt.outputStructure)
        protSCORCH2.inputPDBligandFiles.set(self.protLeDock.outputSmallMolecules)

        self.proj.launchProtocol(protSCORCH2, wait=True)
        return protSCORCH2

    def test(self):
        protSCORCH2 = self._runSCORCH2()
        self._waitOutput(protSCORCH2, 'molecules', sleepTime=10)
        (self.assertIsNotNone, getattr(protSCORCH2, 'molecules', None))


