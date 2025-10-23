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

from pwem.protocols import ProtImportPdb
from pyworkflow.tests import setupTestProject, DataSet

# Scipion chem imports
from pwchem.protocols import ProtChemImportSmallMolecules, ProtChemOBabelPrepareLigands, ProtChemPrepareReceptor
from pwchem.tests import TestImportSequences
import pyworkflow.tests as tests
from ..protocols import ProtocolSCORCH2
from ..protocols.General.protocol_converter import ConvertStructures


class TestSCORCH2(TestImportSequences):
    @classmethod
    def setUpClass(cls):
        super().setUpClass()
        cls.ds = DataSet.getDataSet('model_building_tutorial')
        cls.dsLig = DataSet.getDataSet("smallMolecules")
        tests.setupTestProject(cls)

        cls._runImportPDB()
        cls._runImportSmallMols()

        #cls._runPrepareLigandsOBabel()
        cls._runConvertStructure()

        cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules')
        cls._waitOutput(cls.targetProt, 'outputStructure')


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
    def _runPrepareReceptorOBabel(cls):
        cls.receptorOBabel = cls.newProtocol(
            ProtChemPrepareReceptor,
            inputAtomStruct=cls.protImportPDB.outputPdb,
            waters=True, HETATM=True, rchains=False, PDBFixer=False)
        cls.proj.launchProtocol(cls.receptorOBabel)

    @classmethod
    def _runConvertStructure(cls):
        cls.targetProt = cls.newProtocol(
            ConvertStructures,
            inputObject=cls.protImportPDB.outputPdb,
            outputFormatTarget=0
        )
        cls.proj.launchProtocol(cls.targetProt)

    @classmethod
    def _runPrepareLigandsOBabel(cls):
        cls.protOBabel = cls.newProtocol(
            ProtChemOBabelPrepareLigands,
            inputType=0, method_charges=0,
            inputSmallMolecules=cls.protImportSmallMols.outputSmallMolecules,
            doConformers=False)

        cls.proj.launchProtocol(cls.protOBabel)


    def _runSCORCH2(self):
        protSCORCH2 = self.newProtocol(ProtocolSCORCH2)

        protSCORCH2.useFeatures.set('False')
        protSCORCH2.inputPDBproteinFile.set(self.targetProt.outputStructure)
        protSCORCH2.inputPDBligandFiles.set(self.protImportSmallMols.outputSmallMolecules)

        self.proj.launchProtocol(protSCORCH2, wait=True)
        return protSCORCH2

    def test(self):
        protSCORCH2 = self._runSCORCH2()
        self._waitOutput(protSCORCH2, '', sleepTime=10)
        extraPath = Path(protSCORCH2._getExtraPath())
        expectedCsv = extraPath / "scorch2_results.tsv"
        expected = os.path.abspath(expectedCsv)

        self.assertTrue(expectedCsv.exists(), f"Expected output TSV not found at: {expected}")
        self.assertGreater(expectedCsv.stat().st_size, 0, "Output TSV is empty")


