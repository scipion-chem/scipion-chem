# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to score docking positions obtained by several software using descriptors and score functions
available in the Open Drug Discovery Toolkit (ODDT, https://github.com/oddt/oddt)

"""
import os, glob

from pyworkflow.protocol import params
from pyworkflow.utils import Message, createLink
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *


class ProtocolPoseBusters(EMProtocol):
    """
    Performs plausibility checks for generated molecule poses. Info about each test: https://posebusters.readthedocs.io/en/latest/api.html
    """
    _label = 'PoseBusters docking tests'
    _tests = [
        'All (run all compatible tests)',
        'Distance geometry (ligand only)',
        'Energy ratio (ligand only)',
        'Flatness (ligand only)',
        'Identity (requires crystal ligand)',
        'Intermolecular distance (requires protein)',
        'RMSD (requires crystal ligand)',
        'Volume overlap (requires protein + crystal ligand)'
    ]


    def __init__(self, **args):
        super().__init__(**args)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('tests', params.EnumParam, choices=self._tests, default=0,
                      label="Select test: ",
                      help=(
                          "Select which PoseBusters test to run.\n\n"
                          "? All: runs all tests compatible with the provided inputs\n"
                          "? Ligand-only tests require only a predicted ligand\n"
                          "? Some tests require a crystal ligand (True molecule)\n"
                          "? Some tests require a protein structure (Protein)"
                      ))
        form.addParam('oneFile', params.BooleanParam, default=True,
                        label="Test on one docked molecule: ",
                        help='Choose whether to run the test on one docked molecule or a set.')

        form.addParam('inputMoleculesSets', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Docked Small Molecules: ",
                      help='Select the docked molecules to be tested.')
        form.addParam('molPred', params.StringParam,
                        label='Predicted molecule: ', condition='oneFile',
                        help='Choose the predicted molecule (docked ligand).')
        form.addParam('inputMoleculesRefSets', params.PointerParam, condition='oneFile',
                      pointerClass='SetOfSmallMolecules',
                      label="Input Reference Small Molecules: ",
                      help='Select the reference molecules to be tested against.')
        form.addParam('molTrue', params.StringParam,
                      label='True molecule: ', condition='tests in [0, 1, 4, 6, 7] and oneFile',
                      help='Choose the ground truth molecule (crystal ligand).')
        form.addParam('molCond', params.PointerParam, pointerClass="AtomStruct",
                      condition='tests in [0, 5, 7] and oneFile',
                       label='Protein:', help='Choose the conditioning molecule (protein).')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.tests.get() == 0:
            self._insertFunctionStep(self.allTestsStep)
        #self._insertFunctionStep(self.createOutputStep)

    def allTestsStep(self): #todo call pypi to run all tests and output in txt file
        # get args

        args = ['bust ']

        if self.oneFile.get():
            #docked molecule
            molPred = self.getSpecifiedMol('pred')
            args.append(os.path.abspath(molPred.getPoseFile()))

            if self.inputMoleculesRefSets.get() is not None:
                #true molecule
                molTrue = self.getSpecifiedMol('true')
                args.append(f'-l {os.path.abspath(molTrue.getPoseFile())}')

            if self.molCond.get() is not None:
                #protein
                args.append(f'-p {os.path.abspath(self.molCond.get().getFileName())}')

        else:
            for dockedMol in self.inputMoleculesSets.get():
                #todo mirar este caso
                args.append(f'{dockedMol.getPoseFile()} ')

        print(args)


    def createOutputStep(self):
        outDir = self.getOutDir()
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        mol = self.getSpecifiedMol()

        outFile = self.getOutMol2(outDir)
        mol.setFileName(outFile)
        newMols.append(mol)

        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["The results have been saved in the extra folder as scorch2.results.tsv.",
                   "Features and normalized features have been saved in the extra/results folder, along with the aggregated results and target results."]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        def _validate(self):
            errors = []

            if self.tests.get() in [4, 6, 7] and not self.molTrue.get():
                errors.append("Selected test requires a crystal ligand (True molecule).")

            if self.tests.get() in [5, 7] and not self.molCond.get():
                errors.append("Selected test requires a protein structure (Protein).")

            return errors

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getSpecifiedMol(self, string):
        myMol = None
        if string == 'pred':
            for mol in self.inputMoleculesSets.get():
                if mol.__str__() == self.molPred.get():
                    myMol = mol.clone()
                    break
        else :
            for mol in self.inputMoleculesRefSets.get():
                if mol.__str__() == self.molTrue.get():
                    myMol = mol.clone()
                    break

        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            return myMol

    def getOutDir(self):
        for d in os.listdir(self._getExtraPath()):
            d = self._getExtraPath(d)
            if os.path.isdir(d) and d.endswith('.acpype'):
                return d
        return None

    def getOutMol2(self, outDir):
        for d in os.listdir(outDir):
            d = os.path.join(outDir, d)
            if d.endswith('.mol2'):
                return d
        return None
