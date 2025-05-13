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

# General imports
import os, glob

# Scipion em imports
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem import Plugin
from pwchem.constants import RDKIT_DIC, WARNLIBBIG
from pwchem.utils import getBaseName, makeSubsets

ATYPE, SIZE, HASCYCLES = 'Contains at least x atom type', 'Contains at least x atoms', \
                                'Contains at least x cycles'
ATKEY, ANKEY, CKEY = 'typeAtom',  'numAtoms', 'numCycles'

scriptName = 'ligand_filter_script.py'

class ProtocolBaseLibraryToSetOfMols(EMProtocol):

    def addInputParams(self, form):
      form.addParam('useLibrary', params.BooleanParam, label='Use library as input : ', default=False,
                    help='Whether to use a SMI library SmallMoleculesLibrary object as input')

      form.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                    label='Input library: ', condition='useLibrary',
                    help="Input Small molecules library to predict")
      form.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules', allowsNull=False,
                    condition='not useLibrary', label="Input  Small Molecules: ",
                    help='Select the molecules to be filtered')
      return form

    def createInputStep(self, nt):
      if self.useLibrary.get():
        inDir = os.path.abspath(self._getTmpPath())
        ligFiles = self.inputLibrary.get().splitInFiles(inDir)
      else:
        ligFiles = [os.path.abspath(mol.getFileName()) for mol in self.inputSmallMolecules.get()]

      inputSubsets = makeSubsets(ligFiles, nt, cloneItem=False)
      for it, fileSet in enumerate(inputSubsets):
        with open(self.getInputFile(it), 'w') as f:
          f.write(' '.join(fileSet))

    def getInputFile(self, it):
      return self._getExtraPath(f'inputLigandFiles_{it}.txt')

    def getInputMolFiles(self, it):
      inFile = self.getInputFile(it)
      with open(inFile) as f:
        molFiles = f.read().strip().split()
      return molFiles

    def _validate(self):
        errors = []
        if self.useLibrary.get() and not self.inputLibrary.get().validateSplit():
            errors.append(WARNLIBBIG)
        return errors

class ProtocolGeneralLigandFiltering(ProtocolBaseLibraryToSetOfMols):
    """
    Filters a set of ligands by some user defined attributes: forbidden / necessary atom types, max/min size...
    """
    _label = 'ligand filtering'
    stepsExecutionMode = params.STEPS_PARALLEL

    def __init__(self, **args):
        EMProtocol.__init__(self, **args)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form = self.addInputParams(form)

        group = form.addGroup('Define filter')
        group.addParam('mode', params.EnumParam, choices=["Remove", "Keep"], label='Mode: ', default=0,
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to remove or keep the entry if it meets the attribute')

        group.addParam('filter', params.EnumParam, choices=[ATYPE, SIZE, HASCYCLES],
                       label='Filter by: ', default=0,
                       help='Define a filter by different attributes')
        group.addParam('filterValue', params.IntParam, label='With value (x): ', default=1,
                       help='Value x for the defined filter')
        group.addParam('atomTypeFilter', params.StringParam, label='Atom type: ', default='B',
                       condition='filter==0',
                       help='Atom type to keep / remove filter')
        group.addParam('scoreFilter', params.StringParam, label='Select attribute: ', default='',
                       condition='filter==3',
                       help='Score stored in the SetOfSmallMolecules or SmalloleculesLibrary to use as filter')

        group.addParam('addFilter', params.LabelParam, label='Add filter expression: ',
                       help='Add filter expression to the list')

        group = form.addGroup('Filter list')
        group.addParam('filterList', params.TextParam, width=70, default='', label='List of filtering expressions: ',
                       help='List of filtering expressions the molecules have to pass\n')

        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
      aSteps = []
      nt = self.numberOfThreads.get()
      if nt <= 1: nt = 2

      iStep = self._insertFunctionStep(self.createInputStep, nt-1, prerequisites=[])
      for it in range(nt-1):
        aSteps += [self._insertFunctionStep(self.filterStep, it, prerequisites=[iStep])]
      self._insertFunctionStep(self.createOutputStep, prerequisites=aSteps)

    def filterStep(self, it):
        '''Filter the Set of Small molecules with the defined filters
        '''
        filDic = self.parseFilter()

        paramsPath = os.path.abspath(self._getExtraPath('inputParams_{}.txt'.format(it)))
        self.writeParamsFile(paramsPath, filDic, it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())
        os.remove(self.getInputFile(it))

    def createOutputStep(self):
      passMolFiles = self.parseOutputFiles()
      if self.useLibrary.get():
          inLib = self.inputLibrary.get()
          outLib = inLib.clone()
          mapDic = outLib.getLibraryMap(inverted=True, fullLine=True)

          oLibFile = self._getPath('outputLibrary.smi')
          with open(oLibFile, 'w') as f:
            for molFile in passMolFiles:
              smiName = getBaseName(molFile)
              f.write(f'{mapDic[smiName]}\n')

          outputLib = inLib.clone()
          outputLib.setFileName(oLibFile)
          self._defineOutputs(outputLibrary=outputLib)
      else:
          inMols = self.inputSmallMolecules.get()
          outputSet = inMols.create(self._getPath())
          for mol in inMols:
            if os.path.abspath(mol.getFileName()) in passMolFiles:
              outputSet.append(mol)

          if len(outputSet) > 0:
              self._defineOutputs(outputSmallMolecules=outputSet)
              self._defineSourceRelation(self.inputSmallMolecules, outputSet)


    # --------------------------- UTILS functions -----------------------------------
    def createElementLine(self):
      keepStr = self.getEnumText('mode')
      filterStr, fValue = self.getEnumText('filter'), self.filterValue.get()

      towrite = f"{keepStr} molecule if {filterStr.replace('x', str(fValue)).lower()} "
      if filterStr == ATYPE:
        towrite += self.atomTypeFilter.get()
      towrite += '\n'
      return towrite

    def parseFilter(self):
        filDic = {ATKEY: [], ANKEY: [], CKEY: []}
        filterStr = self.filterList.get().strip()
        for fil in filterStr.split('\n'):
            fAction, fVal = fil.split()[0], float(fil.split('at least')[1].split()[0])
            if 'atom type' in fil:
              fType, atomType = ATKEY, fil.split()[-1]
            elif 'atoms' in fil:
              fType, atomType = ANKEY, None
            elif 'cycles' in fil:
              fType, atomType = CKEY, None

            filDic[fType].append([fAction, fVal, atomType])
        return filDic

    def writeParamsFile(self, paramsFile, filDic, it):
      molFiles = self.getInputMolFiles(it)
      with open(paramsFile, 'w') as f:
        f.write(f'ligandFiles:: {" ".join(molFiles)}\n')
        f.write('filters:: {}\n'.format(filDic))

        f.write('outputPath:: {}\n'.format(os.path.abspath(self._getExtraPath('passMolecules_{}.txt'.format(it)))))
      return paramsFile

    def parseOutputFiles(self):
      passFileNames = []
      for file in glob.glob(self._getExtraPath('passMolecules_*')):
        with open(file) as f:
          f.readline()
          for line in f:
            passFileNames.append(line.strip())
      return passFileNames

    def performScoreFilter(self, molSet, filtList):
      '''Filters the molSet based on score filters.
      '''
      if self.useLibrary.get():
        headers = self.inputLibrary.get().getHeaders()

      for action, value, scoreName in filtList:
        if self.useLibrary.get():
          scIdx = headers.index(scoreName)

        filtSet = molSet.copy()
        for mol in molSet:
          if not self.useLibrary.get():
            molScore = getattr(mol, scoreName).get()
          else:
            molScore = self.parseFileScore(mol, scIdx)

          if (molScore >= value and action == 'Remove') or (molScore < value and action == 'Keep'):
            filtSet.remove(mol)
        molSet = filtSet.copy()

      return molSet

    def parseFileScore(self, file, column):
      with open(file) as f:
        line = f.readline()
        molScore = float(line.split()[column])
      return molScore
