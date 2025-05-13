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
import os

# Scipion em imports
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem.utils import concatThreadFiles, splitFile, getBaseFileName, getBaseName

class ProtocolLibraryFiltering(EMProtocol):
    """
    Filters a small molecules library by some user defined attributes.
    """
    _label = 'library filtering'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                      label='Input library: ', help="Input Small molecules library to filter")

        group = form.addGroup('Define filter')
        group.addParam('scoreFilter', params.StringParam, label='Select attribute: ', default='',
                       help='Score stored in the SmallMoleculesLibrary to use as filter')
        group.addParam('mode', params.EnumParam, choices=["Over", "Below"], label='Mode: ', default=0,
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to keep the molecule if the values of the selected attribute are below/over '
                            'the determined threshold')
        group.addParam('filterValue', params.FloatParam, label='Threshold: ', default=1,
                       help='Threshold for the defined filter')

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

    def createInputStep(self, nt):
      inDir = self._getTmpPath()
      libFile = os.path.abspath(self.inputLibrary.get().getFileName())
      splitFile(libFile, n=nt, oDir=inDir, remove=False)

    def filterStep(self, it):
        '''Filter the smiFile with the defined filters
        '''
        filList = self.parseFilter()
        smiFile = self.getSmiFile(it)
        self.performScoreFilter(smiFile, filList)

    def createOutputStep(self):
      inLib = self.inputLibrary.get()
      oLibFile = self._getPath(getBaseFileName(inLib.getFileName()))
      inDir = os.path.abspath(self._getExtraPath())
      concatThreadFiles(oLibFile, inDir=inDir)

      outputLib = inLib.clone()
      outputLib.setFileName(oLibFile)
      self._defineOutputs(outputLibrary=outputLib)

    # --------------- INFO functions -------------------------
    def _citations(self):
        return []

    # --------------------------- UTILS functions -----------------------------------
    def createElementLine(self):
      overlow = self.getEnumText('mode').lower()
      fAttribute, fValue = self.scoreFilter.get(), self.filterValue.get()

      towrite = f"Keep molecule if {fAttribute} is {overlow} threshold {fValue}\n"
      return towrite

    def getSmiFile(self, it):
      inFile = self.inputLibrary.get().getFileName()
      return self._getTmpPath(f'{getBaseName(inFile)}_{it+1}.smi')

    def getStringBetween(self, line, before, after):
      return line.split(before)[1].split(after)[0].strip()

    def parseFilter(self):
        filList = []
        filterStr = self.filterList.get().strip()
        for fil in filterStr.split('\n'):
            fAction = self.getStringBetween(fil, ' is ', ' threshold ')
            fVal = float(self.getStringBetween(fil, ' threshold ', '\n'))
            fAttr = self.getStringBetween(fil, ' if ', ' is ')

            filList.append([fAction, fVal, fAttr])
        return filList

    def getInputAttributes(self):
      attrs = self.inputLibrary.get().getHeaders()
      return attrs

    def performScoreFilter(self, smiFile, filtList):
      '''Filters the molSet based on score filters.
      '''
      headers = self.inputLibrary.get().getHeaders()
      oSmiFile = self._getExtraPath(getBaseFileName(smiFile))
      with open(oSmiFile, 'w') as f:
        with open(smiFile) as fIn:
          for line in fIn:
            linePass = True
            sline = line.split()
            for action, value, scoreName in filtList:
              scIdx = headers.index(scoreName)
              if (float(sline[scIdx]) > value and action == 'below') or (float(sline[scIdx]) < value and action == 'over'):
                linePass = False

            if linePass:
              f.write(line)

      return oSmiFile
