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
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC

from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
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

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        # todo do i give this option?
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
                      pointerClass='SetOfSmallMolecules', allowsNull=True,
                      label="Input Reference Small Molecules: ",
                      help='Select the reference molecules to be tested against.')
        form.addParam('molTrue', params.StringParam,
                      label='True molecule: ', condition='tests in [0, 1, 4, 6, 7] and oneFile',
                      help='Choose the ground truth molecule (crystal ligand).')
        form.addParam('molCond', params.PointerParam, pointerClass="AtomStruct",
                      condition='tests in [0, 5, 7] and oneFile', allowsNull=True,
                       label='Protein:', help='Choose the conditioning molecule (protein).')
        group = form.addGroup('Scoring function')
        group.addParam('fullReport', params.BooleanParam, default=True,
                      label="Output full report: ",
                      help='Choose whether to output full report or not.')
        group.addParam('outputFormat', params.EnumParam, choices=['short', 'long', 'csv'], default=2,
                       label="Output format: ",
                       help='Choose whether to output full report or not.')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.tests.get() == 0:
            self._insertFunctionStep(self.allTestsStep)
        else:
            pass
            #todo differentiate between tests and run them through the repo classes
        self._insertFunctionStep(self.createOutputStep)

    def allTestsStep(self):
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
                args.append(f'{dockedMol.getPoseFile()} ')

        outputFormat = self.getEnumText('outputFormat')
        args.append(f'--outfmt {outputFormat}')

        if (self.fullReport.get()):
            args.append('--full-report')

        if (self.outputFormat.get() == 2):
            resultsFile = self._getPath('results.csv')
        else:
            resultsFile = self._getPath('results.txt')
        args.append(f'--output {os.path.abspath(resultsFile)}')

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=SCORCH2_DIC,
            program="",
            cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
        )

    def createOutputStep(self):
        if (self.outputFormat.get() == 2):
            resultsFile = self._getPath('results.csv')
        else:
            resultsFile = self._getPath('results.txt')

        newMols = SetOfSmallMolecules.createCopy(self.inputMoleculesSets.get(), self._getPath(), copyInfo=True)

        predPose = self.getSpecifiedMol('pred').getPoseFile()

        for mol in self.inputMoleculesSets.get():
            mol.PoseBusters_file = String()
            if self.oneFile.get():
                if mol.getPoseFile() == predPose:
                    mol.setAttributeValue('PoseBusters_file', resultsFile)
            else:
                mol.setAttributeValue('PoseBusters_file', resultsFile)
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
