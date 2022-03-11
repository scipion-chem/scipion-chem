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
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re

from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes

scriptName = 'scores_docking_oddt.py'
VINA, RFSCORE, NNSCORE, PLECSCORE = 0, 1, 2, 3

class ProtocolScoreDocking(EMProtocol):
    """
    Executes the scoring of a set of molecules which have been previously docked.
    """
    _label = 'ODDT score docking'
    actionChoices = ['MaxScore', 'MinEnergy']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMoleculesSets', params.MultiPointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Docked Small Molecules: ",
                       help='Select the docked molecules to be scored')
        group = form.addGroup('Scoring function')
        group.addParam('scoreChoice', params.EnumParam, default=VINA, label='Score to calculate: ',
                       choices=['Vina', 'RFScore', 'NNScore', 'PLECscore'],
                       help="Name of the score to calculate. \nIf the model has been trained previously, it is loaded "
                            "from {}".format(Plugin.getODDTModelsPath()))

        group.addParam('scoreVersionRF', params.EnumParam, default=0, label='Score version: ',
                       choices=['1', '2', '3'], condition='scoreChoice == {}'.format(RFSCORE),
                       help="Version of the score to calculate in RFScore")
        group.addParam('scoreVersionPLEC', params.EnumParam, default=0, label='Score version: ',
                       choices=['linear', 'nn', 'rf'], condition='scoreChoice == {}'.format(PLECSCORE),
                       help="Version of the score to calculate in PLECScore")

        group.addParam('trainData', params.EnumParam, default=5, label='PDBBind training dataset: ',
                      choices=[2007, 2012, 2013, 2014, 2015, 2016], condition='scoreChoice != {}'.format(VINA),
                      help="PDBBind dataset to train the model into", expertLevel=params.LEVEL_ADVANCED)
        group.addParam('rfSpr', params.IntParam, default=0, label='Minimum contacts: ',
                       condition='scoreChoice == {}'.format(RFSCORE), expertLevel=params.LEVEL_ADVANCED,
                       help="The minimum number of contacts in each pair of atom types in the training set for the "
                            "column to be included in training.\nThis is a way of removal of not frequent and "
                            "empty contacts.")
        group.addParam('depthProt', params.IntParam, default=5, label='ECFP env protein depth: ',
                       condition='scoreChoice == {}'.format(PLECSCORE), expertLevel=params.LEVEL_ADVANCED,
                       help="The depth of ECFP environments generated on the protein side of interaction.\n"
                            "By default 6 (0 to 5) environments are generated.")
        group.addParam('depthLig', params.IntParam, default=1, label='ECFP env ligand depth: ',
                       condition='scoreChoice == {}'.format(PLECSCORE), expertLevel=params.LEVEL_ADVANCED,
                       help="The depth of ECFP environments generated on the ligand side of interaction.\n"
                            "By default 2 (0 to 1) environments are generated.")
        group.addParam('fingerSize', params.IntParam, default=65536, label='Folded PLEC fingerprint size: ',
                       condition='scoreChoice == {}'.format(PLECSCORE), expertLevel=params.LEVEL_ADVANCED,
                       help="The final size of a folded PLEC fingerprint. This setting is not used to limit the data "
                            "encoded in PLEC fingerprint (for that tune the depths), but only the final lenght.\n"
                            "Setting it to too low value will lead to many collisions.")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('scoringStep')
        self._insertFunctionStep('createOutputStep')

    def scoringStep(self):
        mols = self.getInputMoleculesList()
        receptorFile = self.inputSmallMoleculesSets[0].get().getProteinFile()
        results = self.scoreDockings(mols, receptorFile)

    def createOutputStep(self):
        usedIds, i = [], 1
        resDic = self.parseResults()
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMoleculesSets[0].get(), self._getPath(), copyInfo=True)
        for mol in self.getInputMoleculesList():
            #Specific attribute name for each score?
            setattr(mol, "_oddtScore", Float(resDic[os.path.abspath(mol.getPoseFile())]))
            #Id management
            if mol.getObjId() in usedIds:
                newId = i
                i += 1
            else:
                newId = mol.getObjId()
                usedIds.append(newId)
                i = newId + 1

            mol.setObjId(newId)
            newMols.append(mol)

        newMols.updateMolClass()
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
        for molSet in self.inputSmallMoleculesSets:
            molSet = molSet.get()
            if not molSet.isDocked():
                validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        if len(self.inputSmallMoleculesSets) > 1:
            warnings.append("The ids of the molecules will be rewritten so they don't collide with those of the "
                            "other input sets")
        protFile = ''
        for molSet in self.inputSmallMoleculesSets:
            molSet = molSet.get()
            if protFile == '':
                protFile = os.path.abspath(molSet.getProteinFile())
            elif not protFile == os.path.abspath(molSet.getProteinFile()):
                validations += ['These sets of small molecules are docked to different proteins.\n'
                                'Be aware that their results might not be comparable']
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getInputMoleculesList(self):
        mols = []
        inputSets = fillEmptyAttributes(self.inputSmallMoleculesSets)
        for molSet in inputSets:
          for mol in molSet.get():
              mols.append(mol.clone())
        return mols

    def scoreDockings(self, molsScipion, receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion, receptorFile)
        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())

    def parseResults(self):
        resDic = {}
        with open(self._getPath('results.tsv')) as f:
            for line in f:
                resDic[line.split()[0]] = line.split()[1]
        return resDic

    def getScoreFunction(self):
        if self.scoreChoice.get() == RFSCORE:
            function = '{}_{}'.format(self.getEnumText('scoreChoice'), self.getEnumText('scoreVersionRF'))
        elif self.scoreChoice.get() == PLECSCORE:
            function = '{}_{}'.format(self.getEnumText('scoreChoice'), self.getEnumText('scoreVersionPLEC'))
        else:
            function = self.getEnumText('scoreChoice')
        return function

    def getModelFileName(self):
        pdbBind = self.getEnumText('trainData')
        if self.scoreChoice.get() == RFSCORE:
            fName = 'RFScore_v{}_pdbbind{}.pickle'.format(self.getEnumText('scoreVersionRF'), pdbBind)
        elif self.scoreChoice.get() == PLECSCORE:
            fName = 'PLEC{}_p{}_l{}_pdbbind{}_s{}.pickle'.format(self.getEnumText('scoreVersionPLEC'),
                                                                 self.depthProt.get(), self.depthLig.get(),
                                                                 pdbBind, self.fingerSize.get())
        elif self.scoreChoice.get() == NNSCORE:
            fName = 'NNScore_pdbbind{}.pickle'.format(pdbBind)
        else:
            fName = ''
        return fName

    def writeParamsFile(self, paramsFile, molsScipion, recFile):
        molFiles = []
        for mol in molsScipion:
          molFiles.append(os.path.abspath(mol.getPoseFile()))

        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            modelFile = os.path.abspath(Plugin.getODDTModelsPath(self.getModelFileName()))
            f.write('function: {}\n'.format(self.getScoreFunction()))
            if os.path.exists(modelFile) and self.scoreChoice.get() != VINA:
                f.write('model: {}\n'.format(modelFile))
            else:
                f.write('saveModel: {}\n'.format(modelFile))
                if self.scoreChoice.get() != VINA:
                    f.write('pdbbind: {}\n'.format(self.getEnumText('trainData')))
                if self.scoreChoice.get() == RFSCORE:
                    f.write('spr: {}\n'.format(self.rfSpr.get()))
                elif self.scoreChoice.get() == PLECSCORE:
                    f.write('depthProt: {}\n'.format(self.depthProt.get()))
                    f.write('depthLig: {}\n'.format(self.depthLig.get()))
                    f.write('fingerSize: {}\n'.format(self.fingerSize.get()))

            f.write('receptorFile: {}\n'.format(os.path.abspath(recFile)))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile









