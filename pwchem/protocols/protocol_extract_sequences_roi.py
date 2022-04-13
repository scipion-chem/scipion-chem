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
This protocol is used to import a set of pockets (of fpocket, p2rank, autoligand) from some files

"""
import os, math
import numpy as np
from scipy.stats import entropy

from Bio import AlignIO
from Bio.Align import AlignInfo

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *
from pwchem import Plugin

SHANNON, SIMPSON, KABAT = 'Shannon Entropy', 'Simpson Diversity Index', 'Wu-kabat Variability coefficient'

def simpson(counts):
  sumi = 0
  N = sum(counts)
  for c in counts:
    sumi += c * (c - 1)
  return 1 - (sumi / (N * (N - 1)))


def kabat(counts):
  return (sum(counts) * len(counts)) / max(counts)

class ProtExtractSeqsROI(EMProtocol):
    """
    Extract a set of ROIs from a set of sequences based on the conservation in a alignment
    """
    _label = 'Extract sequences ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSequences', params.PointerParam, pointerClass='SetOfSequences',
                      allowsNull=False, label="Input sequences: ",
                      help='Select the set of sequences object where the ROI will be defined')

        group = form.addGroup('Conservation measure')
        group.addParam('method', params.EnumParam,
                       choices=[SHANNON, SIMPSON, KABAT],
                       label='Method to measure conservation: ',
                       help='Method to measure the conservation in each position of the sequences.\n'
                            'http://imed.med.ucm.es/PVS/pvs-help.html#vmth')

        group = form.addGroup('Cluster regions')
        group.addParam('direction', params.EnumParam, default=0,
                       choices=['High conservation', 'Low conservation'],
                       label='Look for regions with: ', display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to look for regions with high or low conservation.'
                            'high conservation will look for low conservation values of the metrics and viceversa')

        group.addParam('highLabel', params.LabelParam, condition='direction==0',
                       label='Keep below threshold',
                       help='High conservation will look for low conservation values of the metrics')
        group.addParam('lowLabel', params.LabelParam, condition='direction==1',
                       label='Keep over threshold',
                       help='Low conservation will look for high conservation values of the metrics')

        group.addParam('thres', params.FloatParam,
                       label='Main threshold for conserved region: ',
                       help='Main threshold that checks the values of the conservation and defines that a position is '
                            'or is not conserved depending on its conservation values.\n'
                            'Conservation values are normalized from 0 (high conservation) to 1 (low conservation)')
        group.addParam('flexThres', params.FloatParam,
                       label='Flexible threshold for conserved region: ',
                       help='This threshold allows an already started conserved region to keep growing if the '
                            'conservation values passes it even if it does not pass the main threshold.\n'
                            'Conservation values are normalized from 0 (high conservation) to 1 (low conservation)')
        group.addParam('numFlex', params.IntParam, default=1,
                       label='Number of consecutive flexible: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Number of consecutive times that the main threshold can be violated and the flexible '
                            'threshold stands.')

        group.addParam('minSize', params.IntParam, default=1,
                       label='Minimum region size: ',
                       help='Minimum size for a region to be considered. Gaps do not count.')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculateConservationStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateConservationStep(self):
        outFile = self.calcConservation()

    def defineOutputStep(self):
        consSeq = str(self.calcConsensus())
        outSeq = self.inputSequences.get().getFirstItem()
        outSeq.setSequence(consSeq)

        roiIdxs = self.getROIs()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for idxs in roiIdxs:
            roi = consSeq[idxs[0]-1:idxs[1]-1]
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
    #Main functions
    def calcConsensus(self):
        inFasta = self._getPath('alignedSequences.fasta')
        inSeqs = self.fillInputSequences()
        for seq in inSeqs:
            seq.appendToFile(inFasta)
        alignment = AlignIO.read(inFasta, 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        return summary_align.gap_consensus(0.5)

    def calcConservation(self):
        seqsArr = self.getSequencesArray()
        if self.getEnumText('method') == SHANNON:
            values = self.calcShannon(seqsArr)
        elif self.getEnumText('method') == SIMPSON:
            values = self.calcSimpson(seqsArr)
        elif self.getEnumText('method') == KABAT:
            values = self.calcKabat(seqsArr)

        print('Conservation: ', values)

        outFile = self.getConservationFile()
        with open(outFile, 'w') as f:
            for value in values:
                f.write('{}\t'.format(value))
        return outFile

    def getROIs(self):
        with open(self.getConservationFile()) as f:
            consValues = f.readline().split()

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

            elif ((v > self.flexThres.get() and self.direction.get() == 1) or \
                    (v < self.flexThres.get() and self.direction.get() == 0)) and inRoi:
                if fails >= self.numFlex.get():
                    rois.append([iniRoi, i + 1])
                    inRoi = False
                else:
                    fails += 1

            elif inRoi:
                rois.append([iniRoi, i + 1])
                inRoi = False
        return rois

    ##################

    def getMaxLenSeq(self, seqSet):
      maxi = 0
      for seq in seqSet:
          leni = len(seq.getSequence())
          if leni > maxi:
              maxi = leni
      return maxi

    def fillToLen(self, listi, leni, val='-'):
        if len(listi) < leni:
            if type(listi) == list:
                listi += [val] * (leni-len(listi))
            elif type(listi) == str:
                listi += val * (leni - len(listi))
        return listi

    def fillInputSequences(self):
        seqSet = []
        inSet = self.inputSequences.get()
        maxLen = self.getMaxLenSeq(inSet)
        for seq in inSet:
            seqSet.append(Sequence(sequence=self.fillToLen(seq.getSequence(), maxLen)))
        return seqSet

    def getElementNumber(self, listi):
        props = []
        for elem in set(listi):
          props.append(list(listi).count(elem))
        return props

    def getSequencesArray(self):
        seqMat = []
        inSeqs = self.fillInputSequences()
        for seq in inSeqs:
            seqMat.append(np.array(list(seq.getSequence())))
        return np.array(seqMat)

    def getConservationFile(self):
        return self._getPath('conservationValues.tsv')

    def calcShannon(self, seqsArr, normalized=True):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            ent = entropy(counts, base=2)
            if normalized:
                ent = ent / math.log2(21)
            entrs.append(ent)
        return entrs

    def calcSimpson(self, seqsArr):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            entrs.append(simpson(counts))
        return entrs

    def calcKabat(self, seqsArr, normalized=True):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            ent = kabat(counts)
            if normalized:
                ent = (ent-1) / 440
            entrs.append(ent)
        return entrs


