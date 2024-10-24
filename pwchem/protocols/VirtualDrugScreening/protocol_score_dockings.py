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
import os, glob, time

from pyworkflow.protocol import params
from pyworkflow.utils import Message, createLink
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem.constants import RDKIT_DIC
from pwchem import Plugin

scriptName = 'scores_docking_oddt.py'
VINA, RFSCORE, NNSCORE, PLECSCORE = 0, 1, 2, 3

class ProtocolScoreDocking(EMProtocol):
    """
    Executes the scoring of a set of molecules which have been previously docked.
    """
    _label = 'ODDT score docking'
    actionChoices = ['MaxScore', 'MinEnergy']

    _paramNames = ['rfSpr', 'depthProt', 'depthLig', 'fingerSize']
    _enumParamNames = ['scoreChoice', 'scoreVersionRF', 'scoreVersionPLEC', 'trainData']
    _defParams = {'scoreChoice': 'Vina', 'scoreVersionRF': '1', 'scoreVersionPLEC': 'linear', 'trainData': '2016',
                  'rfSpr': 0, 'depthProt': 5, 'depthLig': 1, 'fingerSize': 65536, 'isReference': False}
    _omitParamNames = {}


    def __init__(self, **args):
        super().__init__(**args)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMoleculesSets', params.MultiPointerParam,
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
                       help="Version of the score to calculate in PLECscore")

        group.addParam('trainData', params.EnumParam, default=5, label='PDBBind training dataset: ',
                       choices=['2007', '2012', '2013', '2014', '2015', '2016'],
                       condition='scoreChoice != {}'.format(VINA),
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

        group = form.addGroup('Summary')
        group.addParam('insertStep', params.StringParam, default='',
                       label='Insert score number: ',
                       help='Insert the defined score into the workflow.\n'
                            'The default (when empty) is the last position')
        group.addParam('summarySteps', params.TextParam, width=120, readOnly=True,
                       label='Summary of score',
                       help='Summary of the defined score. \nManual modification will have no '
                            'effect, use the wizards to add / delete the scores')
        group.addParam('deleteStep', params.StringParam, default='',
                       label='Delete score number: ',
                       help='Delete the score of the specified index from the workflow.')
        group.addParam('workFlowSteps', params.TextParam, label='User transparent', condition='False')
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def getnThreads(self):
        '''Get the number of threads available for each scoring execution'''
        nPockets = len(self.workFlowSteps.get().strip().split('\n'))
        nThreads = self.numberOfThreads.get() // nPockets
        nThreads = 1 if nThreads == 0 else nThreads
        return nThreads

    def _insertAllSteps(self):
        sSteps, wSteps = [], []
        self.createGUISummary()

        cStep = self._insertFunctionStep('convertInputStep', prerequisites=[])
        #Performing every score listed in the form
        for i, wStep in enumerate(self.workFlowSteps.get().strip().split('\n')):
            if wStep.strip():
                if not wStep in wSteps:
                    sSteps.append(self._insertFunctionStep('scoringStep', wStep, i+1, prerequisites=[cStep]))
                    wSteps.append(wStep)

        self._insertFunctionStep('createOutputStep', prerequisites=sSteps)

    def convertInputStep(self):
        #Convert ligands
        allMols = self.getAllInputMols()
        outDir = self.getInputMolsDir()
        os.mkdir(outDir)

        maeMols, otherMols = getMAEMoleculeFiles(allMols)
        if len(maeMols) > 0:
            try:
                from pwchemSchrodinger.utils.utils import convertMAEMolSet
                convertMAEMolSet(maeMols, outDir, self.numberOfThreads.get(), updateSet=False, subset=False)
            except ImportError:
                print('Conversion of MAE input files could not be performed because schrodinger plugin is not installed')

        for mol in otherMols:
            molFile = mol.getPoseFile()
            molExt = os.path.splitext(molFile)[-1]
            molBaseFile = os.path.join(outDir, mol.getUniqueName() + molExt)
            createLink(molFile, molBaseFile)

        #Convert receptor
        recFile = self.getInputReceptorFile()
        if recFile.endswith('.mae') or recFile.endswith('.maegz'):
            outFn = self.getConvertReceptorFile()
            try:
                from pwchemSchrodinger.utils.utils import convertReceptor2PDB
                convertReceptor2PDB(recFile, os.path.abspath(outFn))
            except ImportError:
                print(
                    'Conversion of MAE input files could not be performed because schrodinger plugin is not installed')
        else:
            createLink(recFile, self.getConvertReceptorFile())


    def performScoring(self, molFiles, molLists, it, receptorFile, msjDic, iS):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams_{}_{}.txt'.format(iS, it)))
        self.writeParamsFile(paramsPath, receptorFile, molFiles, msjDic, iS, it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath(), popen=True)

    def scoringStep(self, wStep, i):
        # Perform the scoring using the ODDT package
        nt = self.getnThreads()
        receptorFile = self.getConvertReceptorFile()
        performBatchThreading(self.performScoring, self.getInputMolFiles(), nt, cloneItem=False,
                              receptorFile=receptorFile, msjDic=eval(wStep), iS=i)

    def createOutputStep(self):
        self.relabelDic = {}
        consensusMols = self.fillEmptyAttributes(self.getAllInputMols())
        consensusMols, idsDic = self.reorderIds(consensusMols)

        sDic = self.getScoresDic()
        newMols = self.defineOutputSet()
        for mol in consensusMols:
            molName = mol.getUniqueName()
            for resIdx, oddtScore in sDic[molName].items():
                setattr(mol, "_oddtScore_{}".format(resIdx), Float(oddtScore))
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
        for pSet in self.inputMoleculesSets:
            pSet = pSet.get()
            if not pSet.isDocked():
                validations.append('Sets of input molecules must be docked first.\n'
                                   'Set: {} has not been docked'.format(pSet))

        return validations

    def _warnings(self):
        warnings = []
        wSteps = []
        for i, wStep in enumerate(self.workFlowSteps.get().split('\n')):
            if wStep.strip():
                if wStep in wSteps:
                    warnings.append('{} Score line is repeated. It will not be performed twice'.format(i+1))
                else:
                    wSteps.append(wStep)

        return warnings

    def createGUISummary(self):
        with open(self._getExtraPath("summary.txt"), 'w') as f:
            f.write(self.createSummary())

    # --------------------------- UTILS functions -----------------------------------
    def getInputMolsDir(self):
        return os.path.abspath(self._getExtraPath('inputMolecules'))

    def getInputMolFiles(self):
        molDir = self.getInputMolsDir()
        molFiles = []
        for molFile in glob.glob(os.path.join(molDir, '*')):
            if '.mae' not in molFile:
                molFiles.append(molFile)
        return molFiles

    def getInputReceptorFile(self):
        return self.inputMoleculesSets[0].get().getProteinFile()

    def getConvertReceptorFile(self):
        recFile = self.getInputReceptorFile()
        recName = getBaseName(recFile)
        return self._getExtraPath(f'{recName}.pdb')

    def defineOutputSet(self):
        inputProteinFile = self.getInputReceptorFile()

        outDocked = SetOfSmallMolecules(filename=self._getPath('outputSmallMolecules.sqlite'))
        outDocked.setDocked(True)
        outDocked.setProteinFile(inputProteinFile)
        return outDocked

    def fillEmptyAttributes(self, inSet):
        '''Fill all items with empty attributes'''
        attributes = self.getAllAttributes(self.inputMoleculesSets)
        for item in inSet:
            for attr in attributes:
                if not hasattr(item, attr):
                    item.__setattr__(attr, attributes[attr])
        return inSet

    def reorderIds(self, inSet):
        '''Return the set with the reordered ids and a mapper dictionary {newId: oldId}'''
        idsDic = {}
        for i, item in enumerate(inSet):
            idsDic[i+1] = item.getObjId()
            item.setObjId(i+1)
        return inSet, idsDic

    def getAllAttributes(self, inSets):
        '''Return a dic with {attrName: ScipionObj=None}'''
        attributes = {}
        for inSet in inSets:
            item = inSet.get().getFirstItem()
            attrKeys = item.getObjDict().keys()
            for attrK in attrKeys:
                if not attrK in attributes:
                    value = item.__getattribute__(attrK)
                    attributes[attrK] = value.clone()
                    attributes[attrK].set(None)
        return attributes

    def getAllInputMols(self):
        mols = []
        for pSet in self.inputMoleculesSets:
            for mol in pSet.get():
                mols.append(mol.clone())
        return mols

    def yieldAllInputMols(self):
        for pSet in self.inputMoleculesSets:
            for mol in pSet.get():
                yield mol.clone()

    def parseResults(self, resFile):
        resDic = {}
        resIdx = resFile.split('/')[-1].split('_')[1]
        with open(resFile) as f:
            for line in f:
                file = line.split()[0]
                resDic[getBaseName(file)] = float(line.split()[1])
        return resIdx, resDic

    def getResultFiles(self):
        rFiles = []
        for file in os.listdir(self._getPath()):
            if file.startswith('results'):
                rFiles.append(self._getPath(file))
        rFiles = natural_sort(rFiles)
        return rFiles

    def getScoreTitles(self, idxs):
        steps = self.summarySteps.get().strip().split('\n')
        titles = []
        for i, st in enumerate(steps):
            if i in idxs:
                titles.append(st.split(':')[1].split('.')[0])
        return titles

    def getScoresDic(self):
        sDic = {}
        for resFile in self.getResultFiles():
            resIdx, resDic = self.parseResults(resFile)
            for fileBase in resDic:
                if fileBase in sDic:
                    sDic[fileBase][resIdx] = resDic[fileBase]
                else:
                    sDic[fileBase] = {resIdx: resDic[fileBase]}
        return sDic

    def getScoreFunction(self, msjDic):
        scoreChoice = msjDic['scoreChoice']
        if scoreChoice == 'RFScore':
            function = '{}_{}'.format(scoreChoice, msjDic['scoreVersionRF'])
        elif scoreChoice == 'PLECscore':
            function = '{}_{}'.format(scoreChoice, msjDic['scoreVersionPLEC'])
        else:
            function = scoreChoice
        return function

    def getModelFileName(self, msjDic):
        pdbBind = msjDic['trainData']
        scoreChoice = msjDic['scoreChoice']
        if scoreChoice == 'RFScore':
            fName = 'RFScore_v{}_pdbbind{}.pickle'.format(msjDic['scoreVersionRF'], pdbBind)
        elif scoreChoice == 'PLECscore':
            fName = 'PLEC{}_p{}_l{}_pdbbind{}_s{}.pickle'.format(msjDic['scoreVersionPLEC'],
                                                                 msjDic['depthProt'], msjDic['depthLig'],
                                                                 pdbBind, msjDic['fingerSize'])
        elif msjDic['scoreChoice'] == 'NNScore':
            fName = 'NNScore_pdbbind{}.pickle'.format(pdbBind)
        else:
            fName = ''
        return fName

    def getPoseFilesStr(self, molFiles, it):
        poseFilesFile = self._getTmpPath(f'poseFiles_{it}.txt')
        if not os.path.exists(poseFilesFile):
            poseFilesStr = ' '.join([molFile for molFile in molFiles])
            with open(poseFilesFile, 'w') as f:
                f.write(poseFilesStr)
        else:
            with open(poseFilesFile) as f:
                poseFilesStr = f.read()

        return poseFilesStr

    def writeParamsFile(self, paramsFile, recFile, molFiles, msjDic, iS, it):
        poseFilesStr = self.getPoseFilesStr(molFiles, it)
        with open(paramsFile, 'w') as f:
            scoreChoice = msjDic['scoreChoice']

            f.write('outputPath: results_{}_{}.tsv\n'.format(iS, it))
            modelFile = os.path.abspath(Plugin.getODDTModelsPath(self.getModelFileName(msjDic)))
            f.write('function: {}\n'.format(self.getScoreFunction(msjDic)))
            if os.path.exists(modelFile) and scoreChoice != 'Vina':
                f.write('model: {}\n'.format(modelFile))
            else:
                f.write('saveModel: {}\n'.format(modelFile))
                if scoreChoice != 'Vina':
                    f.write('pdbbind: {}\n'.format(msjDic['trainData']))
                if scoreChoice == 'RFScore':
                    f.write('spr: {}\n'.format(msjDic['rfSpr']))
                elif scoreChoice == 'PLECscore':
                    f.write('depthProt: {}\n'.format(msjDic['depthProt']))
                    f.write('depthLig: {}\n'.format(msjDic['depthLig']))
                    f.write('fingerSize: {}\n'.format(msjDic['fingerSize']))

            f.write('receptorFile: {}\n'.format(os.path.abspath(recFile)))
            f.write('ligandFiles: {}\n'.format(poseFilesStr))

        return paramsFile

####################### Listing functins ######################

    def countSteps(self):
        stepsStr = self.summarySteps.get() if self.summarySteps.get() is not None else ''
        steps = stepsStr.split('\n')
        return len(steps) - 1

    def createSummary(self):
        '''Creates the displayed summary from the internal state of the steps'''
        sumStr = ''
        for i, dicLine in enumerate(self.workFlowSteps.get().split('\n')):
            if dicLine.strip():
                msjDic = eval(dicLine.strip())
                msjDic = self.addDefaultForMissing(msjDic)
                scoreChoice = msjDic['scoreChoice']
                sumStr += '{}) Score: {}'.format(i+1, scoreChoice)

                if scoreChoice == 'RFScore':
                    sumStr += ', version {}'.format(msjDic['scoreVersionRF'])
                elif scoreChoice == 'PLECscore':
                    sumStr += ', version {}'.format(msjDic['scoreVersionPLEC'])

                if scoreChoice != 'Vina':
                    sumStr += '. PDBbind {}'.format(msjDic['trainData'])
                sumStr += '\n'
        return sumStr

    def addDefaultForMissing(self, msjDic):
        '''Add default values for missing parameters in the msjDic'''
        for pName in [*self._paramNames, *self._enumParamNames]:
            if not pName in msjDic:
                msjDic[pName] = self._defParams[pName]
        return msjDic

    def createMSJDic(self):
        msjDic = {}
        for pName in self._paramNames:
            if hasattr(self, pName):
                msjDic[pName] = getattr(self, pName).get()
            else:
                print('Something is wrong with parameter ', pName)

        for pName in self._enumParamNames:
            if hasattr(self, pName):
                msjDic[pName] = self.getEnumText(pName)
            else:
                print('Something is wrong with parameter ', pName)
        return msjDic

    def getStageParamsDic(self, type='All'):
        '''Return a dictionary as {paramName: param} of the stage parameters of the formulary.
        Type'''
        paramsDic = {}
        for paramName, param in self._definition.iterAllParams():
            if not paramName in self._omitParamNames and not isinstance(param, params.Group) and not isinstance(param,
                                                                                                                params.Line):
                if type == 'All':
                    paramsDic[paramName] = param
                elif type == 'Enum' and isinstance(param, params.EnumParam):
                    paramsDic[paramName] = param
                elif type == 'Normal' and not isinstance(param, params.EnumParam):
                    paramsDic[paramName] = param
        return paramsDic








