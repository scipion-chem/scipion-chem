# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This protocol is used to perform tests oer predicted molecule poses (https://github.com/maabuu/posebusters).

"""

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwem.convert import cifToPdb


class ProtocolPoseBusters(EMProtocol):
    """
    Performs plausibility checks for generated molecule poses.
    More info about the tests: https://posebusters.readthedocs.io/en/latest/cli.html
    """
    _label = 'PoseBusters docking tests'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
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
        form.addParam('useTrueMol', params.BooleanParam, default=False,
                      label="Use true molecule: ",
                      help='Choose whether to use the true molecule in the analysis.')
        form.addParam('inputMoleculesRefSets', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=True,
                      label="Input Reference Small Molecules: ", condition='useTrueMol',
                      help='Select the reference molecules to be tested against.')
        form.addParam('molTrue', params.StringParam,
                      label='True molecule: ', condition='useTrueMol',
                      help='Choose the ground truth molecule (crystal ligand).')
        form.addParam('molCond', params.BooleanParam, default=False,
                       label="Use protein to run tests: ",
                       help='Choose whether to output use protein to run tests. (If True, it will be used to run all default tests that require the protein as input.)')

        form.addParam('fullReport', params.BooleanParam, default=True,
                      label="Output full report: ",
                      help='Choose whether to output full report or not.')
        form.addParam('outputFormat', params.EnumParam, choices=['short', 'long', 'csv'], default=2,
                       label="Output format: ",
                       help='Choose whether to output full report or not.')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.poseBustersStep)
        self._insertFunctionStep(self.createOutputStep)

    def poseBustersStep(self):
        args = ['bust ']

        if self.oneFile.get():
            # docked molecule
            molPred = self.getSpecifiedMol('pred')
            inpFile = self.convertFormat(molPred)
            args.append(os.path.abspath(inpFile))

            if self.useTrueMol.get():
                # true molecule
                molTrue = self.getSpecifiedMol('true')
                inpFile = self.convertFormat(molTrue, type='crystal')
                args.append(f'-l {os.path.abspath(inpFile)}')

            if self.molCond.get():
                # protein
                inpFile = self.convertFormat(self.inputMoleculesSets.get().getProteinFile(), type='file')
                args.append(f'-p {os.path.abspath(inpFile)}')

        else:
            for dockedMol in self.inputMoleculesSets.get():
                inpFile = self.convertFormat(dockedMol)
                args.append(os.path.abspath(inpFile))

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

        predPose = self.getSpecifiedMol('pred')
        predPoseFile = predPose.getPoseFile() if predPose else None
        for mol in self.inputMoleculesSets.get():
            newMol = mol.clone()
            newMol.PoseBusters_file = String()
            if self.oneFile.get():
                if newMol.getPoseFile() == predPoseFile:
                    newMol.setAttributeValue('PoseBusters_file', resultsFile)
            else:
                newMol.setAttributeValue('PoseBusters_file', resultsFile)
            newMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputMoleculesSets.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def convertFormat(self, molPred, type=''):
        if type in ('AtomStruct', 'crystal'):
            basename = os.path.basename(molPred.getFileName()).split('.')[0]
            file = molPred.getFileName()
        elif type == 'file':
            file = molPred
            basename = os.path.splitext(os.path.basename(file))[0]
        else:
            basename = os.path.basename(molPred.getPoseFile()).split('.')[0]
            file = molPred.getPoseFile()

        if file.endswith('.cif'):
            inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            cifToPdb(file, inpFile)
        elif (file.endswith('.pdbqt')):
            inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            pdbqt2other(self, file, inpFile)
        else:
            inpFile = file
        return inpFile

    def getSpecifiedMol(self, string, one=False):
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
