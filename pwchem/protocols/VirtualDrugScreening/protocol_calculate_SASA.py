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
This protocol is used to calculate the SASA of an AtomStruct

"""
import os, json

from pyworkflow.utils import Message
from pyworkflow.protocol import params
from pwem.convert.atom_struct import toPdb, toCIF, AtomicStructHandler, addScipionAttribute
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem import Plugin as pwchemPlugin
from pwchem.objects import SetOfStructROIs, StructROI, SetOfSequenceROIs, Sequence, SequenceROI
from pwchem.utils import *
from pwchem.constants import MGL_DIC

class ProtCalculateSASA(EMProtocol):
    """
    Calculate SASA of an AtomStruct
    """
    _label = 'Calculate SASA'
    _ATTRNAME = 'SASA'
    _OUTNAME = 'outputAtomStruct'
    _possibleOutputs = {_OUTNAME: AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      allowsNull=False, label="Input AtomStruct: ",
                      help='Select the AtomStruct object where the structural ROIs will be defined')
        form.addParam('mapStructure', params.BooleanParam, default=False,
                      label='Map SASA to structure: ',
                      help='Whether to map the SASA into an output AtomStruct which can be visualized')

        form.addParam('extractRegions', params.BooleanParam, default=False,
                      label='Extract sequence conservation regions: ',
                      help='Whether to output sequence regions with the conservation values over/under a '
                           'defined threshold')

        group = form.addGroup('Sequence accessibility regions')
        group.addParam('chain_name', params.StringParam, condition='extractRegions',
                       label='Structure chain:', help="Chain of the structure to extract regions from")
        group.addParam('direction', params.EnumParam, default=0, condition='extractRegions',
                       choices=['Low accessibility', 'High accessibility'],
                       label='Look for regions with: ', display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to look for regions with high or low accessibility.')

        group.addParam('highLabel', params.LabelParam, condition='extractRegions and direction==0',
                       label='Keep below threshold',
                       help='Will look for high accessibility values')
        group.addParam('lowLabel', params.LabelParam, condition='extractRegions and direction==1',
                       label='Keep over threshold',
                       help='Will look for low accessibility values')

        group.addParam('thres', params.FloatParam, condition='extractRegions',
                       label='Main threshold for accessibility region: ',
                       help='Main threshold that checks the values of the accessibility and defines that a position is '
                            'or is not accessible depending on its accessibility values.')
        group.addParam('flexThres', params.FloatParam, condition='extractRegions',
                       label='Flexible threshold for accessibility region: ', allowsNull=True,
                       help='This threshold allows an already started accessible region to keep growing if the '
                            'accessibility values passes it even if it does not pass the main threshold.')
        group.addParam('numFlex', params.IntParam, default=1,
                       label='Number of consecutive flexible: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Number of consecutive times that the main threshold can be violated and the flexible '
                            'threshold stands.')
        group.addParam('minSize', params.IntParam, default=1,
                       label='Minimum region size: ', condition='extractRegions',
                       help='Minimum size for a region to be considered. Gaps do not count.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculateSASAStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateSASAStep(self):
        outSASA = self.getSASAFile()
        if not os.path.exists(outSASA):
            inputAS = self.inputAtomStruct.get().getFileName()
            # Removes HETATM (sometimes yiled error in Biopython)
            clean_PDB(inputAS, self._getExtraPath('inpdb.pdb'), waters=True, HETATM=True)

            args = '{} {}'.format(os.path.abspath(self._getExtraPath('inpdb.pdb')), outSASA)
            pwchemPlugin.runScript(self, 'calculate_SASA.py', args, env='plip', cwd=self._getExtraPath())

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        outStructFileName = self._getPath('outputStructureSASA.cif')
        # Write conservation in a section of the output cif file
        ASH = AtomicStructHandler()

        if self.mapStructure.get():
            sasaDic = self.getSASADic()
            inpAS = toCIF(self.inputAtomStruct.get().getFileName(), self._getTmpPath('inputStruct.cif'))
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, sasaDic, self._ATTRNAME)
            ASH._writeLowLevel(outStructFileName, cifDic)

            AS = AtomStruct(filename=outStructFileName)
            self._defineOutputs(outputAtomStruct=AS)

        if self.extractRegions.get():
            inFile = self.inputAtomStruct.get().getFileName()
            ASH.read(inFile)
            outStr = str(ASH.getSequenceFromChain(modelID=0, chainID=self.getChain()))
            outSeq = Sequence(name='{}_{}'.format(getBaseFileName(inFile), self.getChain()), sequence=outStr)

            roiIdxs = self.getROIs()
            outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
            for idxs in roiIdxs:
              roi = outStr[idxs[0] - 1:idxs[1] - 1]
              if len(roi.replace('-', '')) >= self.minSize.get():
                roiSeq = Sequence(sequence=roi, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs))
                seqROI = SequenceROI(sequence=outSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
                outROIs.append(seqROI)

            if len(outROIs) > 0:
              self._defineOutputs(outputROIs=outROIs)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getSASAFile(self):
        return os.path.abspath(self._getPath('sasa.txt'))

    def getSASADic(self):
        sasaDic = {}
        with open(self.getSASAFile()) as fIn:
            for line in fIn:
                sasaDic[line.split()[0]] = float(line.split()[1])
        return sasaDic

    def getChain(self):
      chainJson = self.chain_name.get()
      if chainJson:
          chain = json.loads(chainJson)['chain']
      else:
          chain = ''
      return chain

    def getAccesibilityValues(self):
      chain = self.getChain()
      if not os.path.exists(self.getSASAFile()):
          self.calculateSASAStep()
      with open(self.getSASAFile()) as f:
        consValues = []
        for line in f:
          spec, value = line.split()
          if spec.split(':')[0][-1] == chain:
            consValues.append(value)
      return list(map(float, consValues))

    def getROIs(self):
        consValues = self.getAccesibilityValues()

        flexThres = self.flexThres.get()
        if not flexThres:
            flexThres = self.thres.get()

        rois = []
        inRoi, fails = False, 0
        for i, v in enumerate(consValues):
            v = float(v)
            if (v > self.thres.get() and self.direction.get() == 1) or \
                    (v < self.thres.get() and self.direction.get() == 0):
                fails = 0
                if not inRoi:
                    inRoi = True
                    iniRoi = i + 1

            elif ((v > flexThres and self.direction.get() == 1) or \
                    (v < flexThres and self.direction.get() == 0)) and inRoi:
                if fails >= self.numFlex.get():
                    rois.append([iniRoi, i + 1])
                    inRoi = False
                else:
                    fails += 1

            elif inRoi:
                rois.append([iniRoi, i + 1])
                inRoi = False
        return rois