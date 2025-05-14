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
This protocol is used to operate sequence ROIs (union, intrsection, difference and symmetric difference)
The operations is applied on the sets iteratively, in the order they are inputted (1 op 2, 12 op 3, 123 op 4, ...)

"""
import sys

from itertools import groupby
from operator import itemgetter

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *

UNION, INTERSECTION, DIFF, SYM_DIFF, BEST = 0, 1, 2, 3, 4

class ProtOperateSeqROI(EMProtocol):
    """
    Implements operations for several sets of sequence ROIs: union, intersection, difference...
    """
    _label = 'Operate sequence ROIs'
    _operations = ['Union', 'Intersection', 'Difference', 'Symmetric_difference', 'Keep best']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputROIsSets', params.MultiPointerParam,
                       pointerClass='SetOfSequenceROIs', allowsNull=False,
                       label="Input Sequence ROIs Sets",
                       help='Select the sequence ROIs sets to operate. If the input is only one set, '
                            'the operation will be performed over its own elements')
        group.addParam('operateIntraSet', params.BooleanParam, label='Operate items intra set: ', default=True,
                       help='Whether to perform the operation "intra set": between each item for each set')
        group.addParam('operateInterSet', params.BooleanParam, label='Operate items inter set: ', default=True,
                       help='Whether to perform the operation "inter set": between items from each different set.'
                            'This operation is performed after the intra set, so if both true, the inter set '
                            'operation will be performed on the resulting ROIs of each intra operated set')

        group = form.addGroup('Operation')
        group.addParam('operation', params.EnumParam, choices=self._operations, display=params.EnumParam.DISPLAY_HLIST,
                       label='Operation: ', default=UNION,
                       help='Operation to perform, in the order of the inputed sets')

        group.addParam('minOverlap', params.FloatParam, label='Minimum overlap proportion: ', default=0,
                       help='Minimum overlap between two ROIs to consider them for the operation.'
                            'The proportion will be calculated as the number of overlapped residues divided by the '
                            'size of the smallest of the two ROIs.\nIf set to 0, '
                            'the operation will be performed as long as there is at least 1 residue overlap')
        group.addParam('bestAttribute', params.StringParam, label='Best based on: ', default='',
                       condition=f'operation=={BEST}',
                       help='Which attribute of the ROI to use to determine best for overlapping ROIs')
        group.addParam('bestDirection', params.BooleanParam, label='Higher is better: ', default=True,
                       condition=f'operation=={BEST}',
                       help='For the selected attribute, whether higher values mean better')

        group.addParam('keepNonOverlaping', params.BooleanParam, label='Keep non overlaping ROIs: ', default=True,
                       help='Whether to include in the final set those ROIs which do not overlap '
                            'with others in the other sets')
        group.addParam('keepAttributes', params.BooleanParam, expertLevel=params.LEVEL_ADVANCED,
                       label='Keep attributes of ROIs: ', default=True,
                       help='Whether to keep the attributes of the output ROIs. The attributes will come from the '
                            'first original ROI the output ROI comes from')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        inSets = self.getInputListOfSets()
        if self.operateIntraSet.get():
            inSets = self.getIntraOperatedROIs(inSets)
        if self.operateInterSet.get():
            inSets = self.getInterOperatedROIs(inSets)
        else:
            oSets = []
            for roiList in inSets:
                oSets += [roi.clone() for roi in roiList]
            inSets = oSets

        if len(inSets) > 0:
            outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
            for i, roi in enumerate(inSets):
                roi.setObjId(i)
                outROIs.append(roi)
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
        if self.operation.get() == BEST:
            for inSet in self.inputROIsSets:
                bAttr = self.bestAttribute.get()
                if not hasattr(inSet.get().getFirstItem(), bAttr):
                    errors.append(f'{inSet.get()} input set has no attribute {bAttr}')

        if not self.operateIntraSet.get() and not self.operateInterSet.get():
            errors.append('You have to choose at least one from intra/inter operation to be true')

        if self.operateInterSet.get() and len(self.inputROIsSets) < 2:
            errors.append('You must specify at least two input sets to operate inter set')
        return errors

    def _warnings(self):
        ws = []
        inpSeqs = []
        for roiSet in self.inputROIsSets:
            inpSeqs.append(roiSet.get().getSequence())
        if len(set(inpSeqs)) > 1:
            ws.append('Not all the sequences in the input sets are the same, this might lead to errors')
        return ws

    # --------------------------- UTILS functions -----------------------------------
    def getInputListOfSets(self):
      inSets = []
      for inPoint in self.inputROIsSets:
        inSets.append(inPoint.get())
      return inSets

    def getIntraOperatedROIs(self, inSets):
        newSets = []
        for curSet in inSets:
            curROIs = [roi.clone() for roi in curSet]
            newSet = [curROIs.pop(0)]
            finalSet = []
            while curROIs:
                for newROI in newSet:
                    reset = False
                    for curROI in curROIs:
                        operatedROIs = self.operateROIs(newROI, curROI, operation=self.getEnumText('operation'))
                        if operatedROIs:
                            newSet.remove(newROI), curROIs.remove(curROI)
                            newSet += operatedROIs
                            reset = True
                            break

                    if reset:
                        break
                    else:
                        finalSet += newSet
                        newSet = [curROIs.pop(0)]

            finalSet += newSet

            newSets.append(finalSet)
        return newSets

    def getInterOperatedROIs(self, inSets):
        prevSet = inSets[0]
        for curSet in inSets[1:]:
            newSet = []
            curOver = [False] * len(curSet)
            prevOver = [False] * len(prevSet)
            for j, prevROI in enumerate(prevSet):
                for i, curROI in enumerate(curSet):
                    newROIs = self.operateROIs(prevROI, curROI, operation=self.getEnumText('operation'))

                    if newROIs:
                        prevOver[j], curOver[i] = True, True
                    for newROI in newROIs:
                        newSet.append(newROI.clone())

            if self.keepNonOverlaping.get():
                for i, curROI in enumerate(curSet):
                    if not curOver[i] and not self.operation.get() == DIFF:
                        newSet.append(curROI.clone())

                for i, prevROI in enumerate(prevSet):
                    if not prevOver[i]:
                        newSet.append(prevROI.clone())
            prevSet = newSet

        return newSet


    def getOperationIdxs(self, roi1, roi2, operation):
        r1, r2 = range(roi1.getROIIdx(), roi1.getROIIdx2() + 1), range(roi2.getROIIdx(), roi2.getROIIdx2() + 1)
        if operation == self._operations[UNION]:
            return list(set(r1).union(r2))

        elif operation == self._operations[INTERSECTION]:
            return list(set(r1).intersection(r2))

        elif operation == self._operations[DIFF]:
            return list(set(r1).difference(r2))

        elif operation == self._operations[SYM_DIFF]:
            return list(set(r1).symmetric_difference(r2))

        elif operation == self._operations[BEST]:
            attr1, attr2 = getattr(roi1, self.bestAttribute.get()), getattr(roi2, self.bestAttribute.get())
            bestROI = roi1 if attr1 >= attr2 and self.bestDirection.get() else roi2
            return list(range(bestROI.getROIIdx(), bestROI.getROIIdx2()+1))

    def operateROIs(self, roi1, roi2, operation):
        newROIs = []
        overIdxs = self.getOperationIdxs(roi1, roi2, operation='Intersection')
        overIdxs.sort()
        overProp = len(overIdxs) / min(roi1.getROILength(), roi2.getROILength())

        if overProp > self.minOverlap.get():
            operIdxs = self.getOperationIdxs(roi1, roi2, operation)
            operIdxs.sort()
            if len(operIdxs) > 0:
                totalSequence = roi1._sequence
                consRanges = self.getConsecutiveRanges(operIdxs)
                for consIdxs in consRanges:
                    roiStr = totalSequence.getSequence()[consIdxs[0]-1:consIdxs[-1]]
                    roiId = 'ROI_{}-{}'.format(consIdxs[0], consIdxs[-1])
                    roiSeq = Sequence(sequence=roiStr, id=roiId)
                    if self.keepAttributes or operation == self._operations[BEST]:
                        roi2Idxs = [roi2.getROIIdx(), roi2.getROIIdx2()]
                        newROI = roi1
                        if consIdxs[0] == roi2Idxs[0] and consIdxs[-1] == roi2Idxs[1]:
                            newROI = roi2
                        setattr(newROI, '_ROISequence', roiSeq)
                        newROI.setROIIdx(consIdxs[0]), newROI.setROIIdx2(consIdxs[-1])
                    else:
                        # todo: option to merge attributes
                        newROI = SequenceROI(sequence=totalSequence, seqROI=roiSeq,
                                             roiIdx=consIdxs[0], roiIdx2=consIdxs[-1])
                    newROIs.append(newROI)

        return newROIs



    def getConsecutiveRanges(self, rang):
        ranges = []
        for k, g in groupby(enumerate(rang), lambda ix : ix[0] - ix[1]):
            ranges.append(list(map(itemgetter(1), g)))
        return ranges


