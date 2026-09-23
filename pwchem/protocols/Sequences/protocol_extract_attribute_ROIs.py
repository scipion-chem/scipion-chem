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
This protocol is used define sequence ROIs based on attribute values stored in a SequenceChem object

"""
import os, math, json
import numpy as np
from scipy.stats import entropy

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb, toCIF, AtomicStructHandler, addScipionAttribute

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *
from pwchem.utils.utilsFasta import pairwiseAlign
from pwchem import Plugin


class ProtExtractSeqsROI(EMProtocol):
    """
    AI Generated:

    This protocol is used to define SequenceROIs based on attribute values stored in a SequenceChem object.

    The protocol scans a selected numerical attribute along a sequence and identifies continuous regions
    that satisfy a threshold condition (above or below). These regions are converted into SequenceROI
    objects with their corresponding sequence coordinates.

    A Gaussian smoothing option can be applied to the attribute profile before ROI detection to reduce
    noise and improve region continuity.

    Overview
    --------
    1. Input a SequenceChem object containing annotated attributes
    2. Select the attribute to analyze
    3. Optionally apply Gaussian smoothing to the attribute signal
    4. Detect continuous regions based on a threshold
    5. Filter regions by minimum size
    6. Generate SequenceROI objects for each detected region

    Input
    -----
    inputSequence:
        SequenceChem object containing sequence and associated attributes.

    inputAttribute:
        Name of the attribute used to define regions (e.g. conservation, variability, score).

    direction:
        Defines selection mode:
        - Below threshold: select low values
        - Over threshold: select high values

    thres:
        Threshold value used to define ROI boundaries.

    minSize:
        Minimum number of residues required for a region to be considered valid.

    performSoftening:
        If True, applies Gaussian smoothing to attribute values before ROI detection.

    wSize:
        Window size for Gaussian smoothing (must be odd).

    gStd:
        Standard deviation for Gaussian kernel.

    Workflow
    --------
    1. Retrieve attribute values from SequenceChem object
    2. Optionally apply Gaussian smoothing
    3. Scan sequence to identify regions above/below threshold
    4. Merge consecutive positions into continuous ROIs
    5. Filter regions by minimum size
    6. Create SequenceROI objects with coordinates

    Output
    ------
    outputROIs:
        SetOfSequenceROIs containing:
        - Extracted sequence regions
        - Start and end coordinates
        - Sequence fragments corresponding to selected attribute regions

    Key Features
    ------------
    - Attribute-driven ROI detection
    - Supports both high and low threshold extraction
    - Gaussian smoothing for noise reduction
    - Minimum region size filtering
    - Compatible with downstream structural and sequence analysis workflows

    Notes
    -----
    - Coordinates are 1-based indexing
    - Gaps are not counted in region size filtering
    - Smoothing uses a Gaussian kernel applied in a sliding window manner
    """
    _label = 'Extract sequences ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSequence', params.PointerParam, pointerClass='SequenceChem',
                      label="Input sequence: ",
                      help='Select sequence you want to extract ROIs from based on its attributes')
        group.addParam('inputAttribute', params.StringParam, label="Input attribute: ",
                       help='Select the attribute you want to extract regions from.')

        group = form.addGroup('Cluster regions')
        group.addParam('direction', params.EnumParam, default=0,
                       choices=['Below threshold ', 'Over threshold'],
                       label='Extract regions: ', display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to extract regions with high or low attribute values '
                            '(over or below the threshold)')

        group.addParam('highLabel', params.LabelParam, condition='direction==0',
                       label='Keep below threshold',
                       help='High variability will look for high conservation values of the metrics')
        group.addParam('lowLabel', params.LabelParam, condition='direction==1',
                       label='Keep over threshold',
                       help='Low variability will look for low conservation values of the metrics')

        group.addParam('thres', params.FloatParam, default=0.5, label='Threshold for extracting regions: ',
                       help='Threshold that checks the values of the attribute and defines that a position is '
                            'or is not extracted depending on them.\n')
        group.addParam('minSize', params.IntParam, default=1, label='Minimum region size: ',
                       help='Minimum size for a region to be considered. Gaps do not count.')

        group.addParam('performSoftening', params.BooleanParam, default=True, label='Perform softening: ',
                       help='Use a Gaussian sliding window filter to soften the attribute values')
        line = group.addLine('Gaussian softening: ', condition='performSoftening',
                             help='Perform Gaussian softening with this parameters')
        line.addParam('wSize', params.IntParam, label="Window size (odd number): ", default=3)
        line.addParam('gStd', params.FloatParam, label="Gaussian deviation: ", default=1)



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        outSeq = self.getInputSequence()
        outStr = outSeq.getSequence()
        
        attrValues = self.getAttributeValues()
        roiIdxs = self.getROIIndexes(attrValues, thres=self.thres.get(), direction=self.direction.get())
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for idxs in roiIdxs:
            roi = outStr[idxs[0]-1: idxs[1]-1]
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

    def _warnings(self):
        warns = []
        if self.performSoftening.get() and self.wSize.get() % 2 != 1:
            warns.append('The size of the softenning window must be an odd number. Using {} instead of {}'.
                          format(self.wSize.get() - 1, self.wSize.get()))
        return warns

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getInputSequence(self):
        return self.inputSequence.get()
    
    def performSoft(self, values, wSize=3, gStd=1):
        from scipy import signal
        if wSize % 2 != 1:
            wSize = wSize - 1
        if wSize <= 0:
            return values

        values = list(map(float, values))
        window, values = signal.windows.gaussian(wSize, std=gStd), np.array(values)
        nValues, size = [], len(window) // 2
        for i in range(len(values)):
            cValues = values[i - size: i + size + 1]
            if i - size < 0:
                # Left border
                cValues = values[:i + size + 1]
                cWindow = window[-len(cValues):]
            elif len(cValues) < len(window):
                #  Right border
                cWindow = window[:len(cValues)]
            else:
                cWindow = window

            cWindow = cWindow / sum(window)
            nValues.append(np.dot(cValues, cWindow))
        return nValues

    def getROIIndexes(self, attrValues, thres, direction):
        roiIdxs = []
        inRoi = False
        for i, v in enumerate(attrValues):
            v = float(v)
            if (v > thres and direction == 1) or (v < thres and direction == 0):
                if not inRoi:
                    inRoi, iniRoi = True, i + 1

            elif inRoi:
                roiIdxs.append([iniRoi, i + 1])
                inRoi = False
        return roiIdxs
    
    def getAttributeValues(self):
        inSeq = self.getInputSequence()
        attrDic = inSeq.getAttributesDic()
        attrValues = attrDic[self.inputAttribute.get()]
        if self.performSoftening.get():
            attrValues = self.performSoft(attrValues, self.wSize.get(), self.gStd.get())

        return attrValues
