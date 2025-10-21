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
from pathlib import Path

from pwem.protocols import ProtImportPdb
from pyworkflow.tests import setupTestProject, DataSet

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules
from pwchem.tests import TestImportSequences, prepRec
from ..protocols import ProtocolSCORCH2


class TestSCORCH2(TestImportSequences):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        setupTestProject(cls)

        cls._runImportPDB()
        cls._runImportSmallMols()

        cls._runPrepareLigandsADT()
        cls._runPrepareReceptorADT()
        if not hasattr(cls, 'protPrepareLigands') or not hasattr(cls, 'protPrepareReceptor'):
            print("Autodock plugin not available ? skipping SCORCH2 test")
            return


    @classmethod
    def _runImportSmallMols(cls):
        cls.protImportSmallMols = cls.newProtocol(
            ProtChemImportSmallMolecules,
            filesPath=cls.dsLig.getFile('pdb'))
        cls.proj.launchProtocol(cls.protImportSmallMols, wait=True)

    @classmethod
    def _runImportPDB(cls):
        protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=1, pdbFile=cls.ds.getFile('PDBx_mmCIF/5ni1.pdb'))
        cls.launchProtocol(protImportPDB)
        cls.protImportPDB = protImportPDB

    @classmethod
    def _runPrepareReceptorADT(cls):
        try:
            from autodock.protocols import ProtChemADTPrepareReceptor
            cls.protPrepareReceptor = cls.newProtocol(
                ProtChemADTPrepareReceptor,
                inputAtomStruct=cls.protImportPDB.outputPdb,
                HETATM=True, rchains=True,
                chain_name=prepRec,
                repair=3)

            cls.launchProtocol(cls.protPrepareReceptor)
        except Exception as ex:
            print(f'Autodock plugin is necesssary to run this test: {ex}')

    @classmethod
    def _runPrepareLigandsADT(cls):
        try:
            from autodock.protocols import ProtChemADTPrepareLigands
            cls.protPrepareLigands = cls.newProtocol(
                ProtChemADTPrepareLigands,
                inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules)

            cls.launchProtocol(cls.protPrepareLigands)
        except Exception as ex:
            print(f'Autodock plugin is necesssary to run this test: {ex}')


    def _runSCORCH2(self):
        protSCORCH2 = self.newProtocol(ProtocolSCORCH2)

        protSCORCH2.useFeatures.set('False')
        protSCORCH2.inputPDBproteinFile.set(self.protPrepareReceptor.outputStructure)
        protSCORCH2.inputPDBligandFiles.set(self.protPrepareLigands.outputSmallMolecules)

        self.proj.launchProtocol(protSCORCH2, wait=True)
        return protSCORCH2

    def test(self):
        protSCORCH2 = self._runSCORCH2()
        self._waitOutput(protSCORCH2, '', sleepTime=10)
        extraPath = Path(protSCORCH2._getExtraPath())
        expectedCsv = extraPath / "scorch2_results.tsv"

        self.assertTrue(expectedCsv.exists(), f"Expected output TSV not found at: {expectedCsv}")
        self.assertGreater(expectedCsv.stat().st_size, 0, "Output TSV is empty")

