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
import os, json
from itertools import groupby
from operator import itemgetter
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, min_dist, residue_depth
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *
from pwchem import Plugin

UNION, INTERSECTION, DIFF, SYM_DIFF = 0, 1, 2, 3

class ProtOperateSeqROI(EMProtocol):
    """
    Implements operations for several sets of sequence ROIs: union, intersection, difference,...
    """
    _label = 'Operate sequence ROIs'
    _operations = ['Union', 'Intersection', 'Difference', 'Symmetric_difference']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        form.addParam('inputROIsSets', params.MultiPointerParam,
                      pointerClass='SetOfSequenceROIs', allowsNull=False,
                      label="Input Sequence ROIs Sets",
                      help='Select the sequence ROIs sets to operate. If the input is only one set, '
                           'the operation will be performed over its own elements')

        group = form.addGroup('Operation')
        group.addParam('operation', params.EnumParam, choices=self._operations,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Operation: ', default=UNION,
                       help='Operation to perform, in the order of the inputed sets')
        group.addParam('keepNonOverlaping', params.BooleanParam,
                       label='Keep non overlaping ROIs: ', default=True,
                       help='Whether to include in the final set those ROIs which do not overlap '
                            'with others in the other sets')
        group.addParam('keepName', params.BooleanParam, expertLevel=params.LEVEL_ADVANCED,
                       label='Keep name of ROIs: ', default=False,
                       help='Whether to keep the name of those output ROIs which totally correspond to '
                            'one in the input')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        operatedROIs = self.getOperatedROIs()

        if len(operatedROIs) > 0:
            outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
            for i, roi in enumerate(operatedROIs):
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
        if len(self.inputROIsSets) < 2:
            errors.append('You must specify at least two input sets to operate')

        return errors

    def _warnings(self):
        ws = []
        inpSeqs = []
        for roiSet in self.inputROIsSets:
            inpSeqs.append(roiSet.get().getSequence())
        if len(set(inpSeqs)) > 1:
            ws.append('Not all the sequences in the input sets are the same, this might lead to errors')
        return  ws

    # --------------------------- UTILS functions -----------------------------------
    def getOperatedROIs(self):
        inSets = self.inputROIsSets

        prevSet = inSets[0]
        for curSet in inSets[1:]:
            newSet = []
            curOver = [False] * len(curSet.get())
            prevOver = [False] * len(prevSet.get())
            for j, prevROI in enumerate(prevSet.get()):
                for i, curROI in enumerate(curSet.get()):
                    newROIs = self.operateROIs(prevROI, curROI)

                    if len(self.getOperationIdxs(prevROI, curROI, operation='Intersection')) > 0:
                        prevOver[j], curOver[i] = True, True

                    for newROI in newROIs:
                        newSet.append(newROI.clone())

            if self.keepNonOverlaping.get():
                for i, curROI in enumerate(curSet.get()):
                    if not curOver[i] and not self.operation.get() == DIFF:
                        newSet.append(curROI.clone())

                for i, prevROI in enumerate(prevSet.get()):
                    if not prevOver[i]:
                        newSet.append(prevROI.clone())

        return newSet


    def getOperationIdxs(self, roi1, roi2, operation=None):
        r1, r2 = range(roi1.getROIIdx(), roi1.getROIIdx2() + 1), range(roi2.getROIIdx(), roi2.getROIIdx2() + 1)
        if (operation is None and self.operation.get() == UNION) or operation==self._operations[UNION]:
            return list(set(r1).union(r2))

        elif (operation is None and self.operation.get() == INTERSECTION) or operation==self._operations[INTERSECTION]:
            return list(set(r1).intersection(r2))

        elif (operation is None and self.operation.get() == DIFF) or operation==self._operations[DIFF]:
            return list(set(r1).difference(r2))

        elif (operation is None and self.operation.get() == SYM_DIFF) or operation==self._operations[SYM_DIFF]:
            return list(set(r1).symmetric_difference(r2))

    def operateROIs(self, roi1, roi2, operation=None):
        overIdxs = self.getOperationIdxs(roi1, roi2, operation='Intersection')
        overIdxs.sort()
        operIdxs = self.getOperationIdxs(roi1, roi2, operation)
        operIdxs.sort()
        if len(overIdxs) > 0 and len(operIdxs) > 0:
            totalSequence = roi1._sequence
            newROIs = []
            for consIdxs in self.getConsecutiveRanges(operIdxs):
                roiStr = totalSequence.getSequence()[consIdxs[0]-1:consIdxs[-1]]
                roiId = 'ROI_{}-{}'.format(consIdxs[0], consIdxs[-1])
                if self.keepName:
                    roi1Idxs, roi2Idxs = [roi1.getROIIdx(), roi1.getROIIdx2()], [roi2.getROIIdx(), roi2.getROIIdx2()]
                    if consIdxs[0] == roi1Idxs[0] and consIdxs[-1] == roi1Idxs[1]:
                        roiId = roi1.getROIId()
                    elif consIdxs[0] == roi2Idxs[0] and consIdxs[-1] == roi2Idxs[1]:
                        roiId = roi2.getROIId()

                roiSeq = Sequence(sequence=roiStr, id=roiId)
                newROI = SequenceROI(sequence=totalSequence, seqROI=roiSeq,
                                     roiIdx=consIdxs[0], roiIdx2=consIdxs[-1])
                newROIs.append(newROI)
            return newROIs
        else:
            return []


    def getConsecutiveRanges(self, rang):
        ranges = []
        for k, g in groupby(enumerate(rang), lambda ix : ix[0] - ix[1]):
            ranges.append(list(map(itemgetter(1), g)))
        return ranges


