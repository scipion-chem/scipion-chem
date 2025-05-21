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
import pandas as pd

# Scipion em imports
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem.utils import splitFile, getBaseFileName, getBaseName, concatFiles, findThreadFiles
from pwchem.objects import SmallMoleculesLibrary

def thresholdFunc(threshold, action='over'):
  """Returns a function that checks if values are above threshold"""

  def filterFunc(value):
    try:
      return float(value) < threshold if action == 'below' else float(value) > threshold
    except (ValueError, TypeError):
      return False  # or handle invalid values differently

  return filterFunc

UNION, INTER, DIFF, FILT, RANK, REM = 0, 1, 2, 3, 4, 5
SMI, NAME = 0, 1

class ProtocolOperateLibrary(EMProtocol):
    """
    Operates a set of small molecules libraries.
    """
    _label = 'operate library'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        group = form.addGroup('Operation')
        group.addParam('operation', params.EnumParam, label='Operation: ', default=UNION,
                       choices=['Union', 'Intersection', 'Difference', 'Filter', 'Ranking', 'Remove columns'],
                       help='-Sets operations: Duplicates share the same reference column value\n'
                            '\tUnique: keep just one item with the same reference column value.'
                            '\n\tUnion: merges two or more sets.\n\tIntersection: keep only the items '
                            'repeated in all the input sets.\n\tDifference: keep only the items in the first library '
                            'that are not present in the rest.'
                            '\n\n-Modification operations:\n\tFilter: outputs only '
                            'those items passing the filter\n\tRemove columns: remove the specified columns\n\t'
                            'Ranking: outputs only the top/bottom elements for the specified column.')
        
        group.addParam('refAttribute', params.EnumParam, label='Reference attribute: ',
                       choices=['SMI', 'Name'], default=SMI, condition=f'operation in [{UNION}, {INTER}, {DIFF}]',
                       help='Attribute of the library to use as reference for the set operations')

        group = form.addGroup('Input')
        group.addParam('inputLibraries', params.MultiPointerParam, pointerClass="SmallMoleculesLibrary",
                       condition=f'operation in [{UNION}, {INTER}, {DIFF}]', allowsNull=True,
                       label='Input libraries: ', help="Input Small molecules libraries to operate as sets")
        group.addParam('inputLibrary', params.PointerParam, pointerClass="SmallMoleculesLibrary",
                       condition=f'operation not in [{UNION}, {INTER}, {DIFF}]', allowsNull=True,
                       label='Input library: ', help="Input Small molecules library to filter")
        
        group = form.addGroup('Define filter', condition=f'operation in [{FILT}, {RANK}, {REM}]')
        group.addParam('filterAttr', params.StringParam, label='Filter attribute: ', default='',
                       condition=f'operation in [{FILT}, {RANK}, {REM}]',
                       help='Attribute to use as filter or to remove')
        group.addParam('mode', params.EnumParam, choices=["Over", "Below"], label='Mode: ', default=0,
                       display=params.EnumParam.DISPLAY_HLIST, condition=f'operation == {FILT}',
                       help='Whether to keep the molecule if the values of the selected attribute are below/over '
                            'the determined threshold')

        group.addParam('filterValue', params.StringParam, label='Ranking / threshold value: ', default=0.5,
                       condition=f'operation in [{FILT}, {RANK}]',
                       help='Ranking mode: outputs the best items for the attribute. The value defined as number of '
                            'outputs (e.g "100" for 100 higher values outputs), proportion from input or percentage '
                            '(e.g: "-0.3" or "-30%" for lower values outputs). If the value is positive, best is '
                            'highest; if negative, best are lowest values.'
                            '\nThreshold mode: threshold for the defined filter')

        group = form.addGroup('Filter list', condition=f'operation == {FILT}')
        group.addParam('addFilter', params.LabelParam, label='Add filter expression: ',
                       help='Add filter expression to the list')
        group.addParam('filterList', params.TextParam, width=70, default='', label='List of filtering expressions: ',
                       help='List of filtering expressions the molecules have to pass\n')

        group = form.addGroup('Remove list', condition=f'operation == {REM}')
        group.addParam('addAttribute', params.LabelParam, label='Add attribute to remove: ',
                       help='Add filter expression to the list')
        group.addParam('removeList', params.TextParam, width=70, default='', label='List of attributes to remove: ',
                       help='List of attributes to remove from the library\n')

        form.addParallelSection(threads=4, mpi=1)

        # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
      aSteps = []
      if self.operation.get() in [FILT, REM]:
          nt = self.numberOfThreads.get()
          if nt <= 1: nt = 2

          iStep = self._insertFunctionStep(self.createInputStep, nt-1, prerequisites=[])
          for it in range(nt-1):
            aSteps += [self._insertFunctionStep(self.filterStep, it, prerequisites=[iStep])]
      
      elif self.operation.get() in [RANK]:
          aSteps = [self._insertFunctionStep(self.rankStep, prerequisites=[])]
      
      else:
        aSteps = [self._insertFunctionStep(self.operateStep, prerequisites=[])]
      self._insertFunctionStep(self.createOutputStep, prerequisites=aSteps)

    def createInputStep(self, nt):
      inDir = self._getTmpPath()
      libFile = os.path.abspath(self.inputLibrary.get().getFileName())
      splitFile(libFile, n=nt, oDir=inDir, remove=False)

    def filterStep(self, it):
        '''Filter the smiFile with the defined filters
        '''
        smiFile = self.getSmiFile(it)
        if self.operation.get() == FILT:
            filList = self.parseFilter()
            self.performScoreFilter(smiFile, filList)
        else:
            remCols = self.removeList.get().strip().split('\n')
            self.performColRemoval(smiFile, remCols)

    
    def rankStep(self):
        import heapq
        inLib = self.inputLibrary.get()
        nMax = self.getRankMax()
        colIdx = inLib.getHeaders().index(self.filterAttr.get())
        mult = -1 if self.filterValue.get().startswith('-') else 1
        
        heap = []
        libFile = os.path.abspath(inLib.getFileName())
        with open(libFile) as fIn:
            for line in fIn:
                parts = line.rstrip('\n').split()
                try:
                    mVal = float(parts[colIdx])
                except (ValueError, IndexError):
                    continue

                if len(heap) < nMax:
                    heapq.heappush(heap, (mult * mVal, line))
                else:
                    heapq.heappushpop(heap, (mult * mVal, line))

        oFile = self.getOutputSmiFile()
        with open(oFile, 'w') as f:
            for _, line in heap:
                f.write(line)
        

    def operateStep(self):
        mapDic = self.inputLibraries[0].get().getLibraryMap(fullLine=True, inverted=self.refAttribute.get() == NAME, lineDic=True)
        for i, inPointer in enumerate(self.inputLibraries):
            if i > 0:
                auxDic = {}
                inLib = inPointer.get()
                for k, lDic in inLib.yieldLibraryMapItems(fullLine=True, inverted=self.refAttribute.get() == NAME, lineDic=True):
                    mapDic, auxDic = self.performOperate(mapDic, auxDic, k, lDic)
                if self.operation.get() == INTER:
                    mapDic = auxDic

        headers = self.getAllHeaders()
        oFile = self.getOutputSmiFile()
        with open(oFile, 'w') as f:
            for lDic in mapDic.values():
                line = f'{lDic["SMI"]}\t{lDic["Name"]}'
                for k in headers:
                    val = lDic[k] if k in lDic else None
                    line += f'\t{val}'
                line += '\n'
                f.write(line)

    def createOutputStep(self):
      headers = self.getAllHeaders()
      oFile = self.getOutputSmiFile()
      if self.operation.get() in [FILT, REM]:
          thFiles = findThreadFiles(self.inputLibrary.get().getFileName(), self._getExtraPath())
          concatFiles(thFiles, oFile, remove=True)
          if self.operation.get() == REM:
            [headers.remove(remCol) for remCol in self.removeList.get().strip().split('\n')]

      outputLib = SmallMoleculesLibrary(libraryFilename=oFile, headers=['SMI', 'Name']+headers)
      outputLib.calculateLength()
      self._defineOutputs(outputLibrary=outputLib)

    # --------------- INFO functions -------------------------
    def _citations(self):
        return []

    # --------------------------- UTILS functions -----------------------------------
    def createElementLine(self):
      overlow = self.getEnumText('mode').lower()
      fAttribute, fValue = self.filterAttr.get(), self.filterValue.get()

      towrite = f"Keep molecule if {fAttribute} is {overlow} threshold {fValue}\n"
      return towrite

    def getRankMax(self):
        nVal = self.filterValue.get()
        if '%' in nVal:
            nVal = float(nVal.replace('%', '')) / 100
        nVal = abs(float(nVal))
        if nVal < 1:
            nInput = self.inputLibrary.get().getLength()
            nMax = int(nInput * nVal)
        else:
            nMax = nVal
        return nMax

    def performOperate(self, mapDic, auxDic, k, lDic):
        '''Update the map dictionary depending on the operation selected
        '''
        if k in mapDic:
            # if element found
            if self.operation.get() == DIFF:
                # Remove if diff
                del mapDic[k]
            else:
                # Combine columns if union or intersect
                mapDic[k].update(lDic)
                if self.operation.get() == INTER:
                    # Auxiliary dic gets only the intersecting keys
                    auxDic[k] = mapDic[k]
        else:
            if self.operation.get() == UNION:
                # Add element if union
                mapDic[k] = lDic

        return mapDic, auxDic

    def getAllHeaders(self):
        if self.operation.get() in [UNION, INTER, DIFF]:
            oHeadFile = self._getTmpPath('headers.txt')
            if os.path.exists(oHeadFile):
                with open(oHeadFile) as f:
                    headers = f.read().split(',')
            else:
                headers = list(set([h for inLib in self.inputLibraries for h in inLib.get().getHeaders()]))
                headers.remove('SMI'), headers.remove('Name')
                with open(oHeadFile, 'w') as f:
                    f.write(f'{",".join(headers)}')
        else:
            headers = self.inputLibrary.get().getHeaders()
            headers.remove('SMI'), headers.remove('Name')

        return headers

    def getOutputSmiFile(self):
        return self._getPath('outputLibrary.smi')

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

    def filterChunk(self, chunk, filtList, headers):
      for action, threshold, scoreName in filtList:
        scIdx = headers.index(scoreName)
        filterFunc = thresholdFunc(threshold, action)
        chunk = chunk[chunk.iloc[:, scIdx].apply(filterFunc)]
      return chunk

    def performScoreFilter(self, smiFile, filtList, chunksize=100000):
      """
      Process a large file in chunks, filtering based on column values

      Args:
          smiFile: Path to input smi file
          filtList: List of filters that the elements of the chunks must pass
          chunksize: Number of rows per chunk
      """
      headers = self.inputLibrary.get().getHeaders()
      oSmiFile = self._getExtraPath(getBaseFileName(smiFile))

      chunkReader = pd.read_csv(smiFile, sep='\s+', chunksize=chunksize, header=None)
      firstChunk = True
      for chunk in chunkReader:
        filteredChunk = self.filterChunk(chunk, filtList, headers)

        if firstChunk:
          filteredChunk.to_csv(oSmiFile, sep='\t', index=False, header=None)
          firstChunk = False
        else:
          filteredChunk.to_csv(oSmiFile, mode='a', sep='\t', header=False, index=False)

    def removeColumsChunk(self, chunk, remCols, headers):
        '''Remove the columns of the pandas dataframe corresponding to the headers in the remCols list
        '''
        remIdxs = [headers.index(scoreName) for scoreName in remCols]
        keepCols = [i for i in range(chunk.shape[1]) if i not in remIdxs]
        return chunk.iloc[:, keepCols]

    def performColRemoval(self, smiFile, remCols, chunksize=100000):
      """
      Process a large file in chunks, filtering based on column values

      Args:
          smiFile: Path to input smi file
          remCols: List of columns to be removed
          chunksize: Number of rows per chunk
      """
      headers = self.inputLibrary.get().getHeaders()
      oSmiFile = self._getExtraPath(getBaseFileName(smiFile))

      chunkReader = pd.read_csv(smiFile, sep='\s+', chunksize=chunksize)
      firstChunk = True
      for chunk in chunkReader:
        filteredChunk = self.removeColumsChunk(chunk, remCols, headers)

        if firstChunk:
          filteredChunk.to_csv(oSmiFile, sep='\t', index=False)
          firstChunk = False
        else:
          filteredChunk.to_csv(oSmiFile, mode='a', sep='\t', header=False, index=False)