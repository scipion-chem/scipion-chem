#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# # -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
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
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import RDKIT_DIC
from pwchem.utils import getBaseName, runInParallel, makeSubsets

ATYPE, SIZE, HASCYCLES = 'Contains at least x atom type', 'Contains at least x atoms', 'Contains at least x cycles'
ATKEY, ANKEY, CKEY = 'typeAtom',  'numAtoms', 'numCycles'

scriptName = 'ligand_filter_script.py'

class ProtocolLigandFiltering(EMProtocol):
    """
    Filters a set of ligands by some user defined attributes: forbidden / necessary atom types, max/min size...
    """
    _label = 'ligand filtering'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input  Small Molecules: ",
                      help='Select the molecules to be filtered')

        group = form.addGroup('Define filter')
        group.addParam('mode', params.EnumParam, choices=["Remove", "Keep"], label='Mode: ', default=0,
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to remove or keep the entry if it meets the attribute')

        group.addParam('filter', params.EnumParam, choices=[ATYPE, SIZE, HASCYCLES], label='Filter by: ', default=0,
                       help='Define a filter by different attributes')
        group.addParam('filterValue', params.IntParam, label='With value (x): ', default=1,
                       help='Value x for the defined filter')
        group.addParam('atomTypeFilter', params.StringParam, label='Atom type: ', default='B',
                       condition='filter==0',
                       help='Atom type to keep / remove filter')

        group.addParam('addFilter',params.LabelParam, label='Add filter expression: ',
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
      inputSubsets = makeSubsets(self.inputSmallMolecules.get(), nt - 1)
      for it, subset in enumerate(inputSubsets):
        aSteps += [self._insertFunctionStep('filterStep', subset, it, prerequisites=[])]
      self._insertFunctionStep('createOutputStep')

    def filterStep(self, molSet, it):
        '''Filter the Set of Small molecules with the defined filters'''
        filDic = self.parseFilter()

        paramsPath = os.path.abspath(self._getExtraPath('inputParams_{}.txt'.format(it)))
        self.writeParamsFile(paramsPath, molSet, filDic, it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())


    def createOutputStep(self):
        outputSet = self.inputSmallMolecules.get().create(self._getPath())
        passMolFiles = self.parseOutputFiles()
        for mol in self.inputSmallMolecules.get():
          if os.path.abspath(mol.getFileName()) in passMolFiles:
            outputSet.append(mol)

        if len(outputSet) > 0:
            self._defineOutputs(outputSmallMolecules=outputSet)
            self._defineSourceRelation(self.inputSmallMolecules, outputSet)

    # --------------- INFO functions -------------------------
    def _citations(self):
        return []

    # --------------------------- UTILS functions -----------------------------------
    def parseFilter(self):
        filDic = {ATKEY: [], ANKEY: [], CKEY: []}
        filterStr = self.filterList.get().strip()
        for fil in filterStr.split('\n'):
            fAction, fVal = fil.split()[0], int(fil.split('at least')[1].split()[0])
            if 'atom type' in fil:
              fType, atomType = 'typeAtom', fil.split()[-1]
            elif 'atoms' in fil:
              fType, atomType = 'numAtoms', None
            elif 'cycles' in fil:
              fType, atomType = 'numCycles', None

            filDic[fType].append([fAction, fVal, atomType])
        return filDic


    def writeParamsFile(self, paramsFile, molsScipion, filDic, it):
      molFiles = []
      with open(paramsFile, 'w') as f:
        for mol in molsScipion:
          molFiles.append(os.path.abspath(mol.getFileName()))

        f.write('ligandFiles:: {}\n'.format(' '.join(molFiles)))
        f.write('filters:: {}\n'.format(filDic))

        f.write('outputPath:: {}\n'.format(os.path.abspath(self._getExtraPath('passMolecules_{}.txt'.format(it)))))

      return paramsFile

    def parseOutputFiles(self):
      passFileNames = []
      for file in glob.glob(self._getExtraPath('passMolecules_*')):
        print(file)
        with open(file) as f:
          f.readline()
          for line in f:
            passFileNames.append(line.strip())
      return passFileNames
