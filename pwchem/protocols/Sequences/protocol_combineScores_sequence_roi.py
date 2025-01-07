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
This protocol is used to filter a set of sequenceROIs based on the values of their selected attributes

"""
import json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSequenceROIs
from pwchem.utils import Float

class ProtCombineScoresSeqROI(EMProtocol):
    """
    Score a set of sequence ROIs based on the values of their selected attributes and combining with information from
    other set of ROIs of the same sequence
    """
    _label = 'Combine sequence ROIs scores'
    _overlapOptions = ['Input ROI contains conditional ROI', 'Conditional ROI contains input ROI', 'Any contains any',
                       'Partial overlap on input ROI', 'Partial overlap on conditional ROI', 'Partial overlap on any']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSequenceROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                      label="Input sequenceROIs: ",
                      help='Select the set of sequenceROIs to be filtered')
        group.addParam('scoreName', params.StringParam, label='Score name: ', default='',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Define a name for the calculated score the ROIs will be filtered by.')


        group = form.addGroup('Define filters')
        group.addParam('selectAttribute', params.StringParam, label='Select filter attribute: ', default='',
                       help='Select the filter attribute you want to filter by')
        group.addParam('attrWeight', params.FloatParam, label='Attribute weight: ', default=1,
                       help='Define a weight for the attribute. The bigger the absolute value of the weight, the more '
                            'important it will be. If negative, the filter will look for smaller values.')
        group.addParam('addFilter', params.LabelParam, label='Add filter: ',
                       help='Add defined filter')

        group = form.addGroup('Filters summary')
        group.addParam('filtSummary', params.TextParam, width=100, label='Filters summary:', default='',
                       help='Summary of the selected filters')

        form.addSection(label='Conditional filters')
        group = form.addGroup('Conditional filters')
        group.addParam('conditionalROIs', params.MultiPointerParam, pointerClass='SetOfSequenceROIs',
                       label="Input conditional sequenceROIs: ",
                       help='Select sets of ROIs containing extra information about the sequence that will '
                            'modify the score of the input ROIs if they overlap. For example, this ROIs '
                            'can contain secondary structural information and you may want to give more value to those'
                            ' input ROIs that overlap with helix structures.')
        group.addParam('scoreCondName', params.StringParam, label='Conditional score name: ', default='',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Define a name for the calculated score the ROIs will be filtered by.')

        group.addParam('inSet', params.StringParam, label='From set: ', default='',
                       help='Add input sequence ROI as epitope from this set')
        group.addParam('condAttribute', params.StringParam, label='Select conditional attribute: ', default='',
                       help='Select the conditional attribute you want to filter by')
        group.addParam('condWeight', params.FloatParam, label='Conditional attribute weight: ', default=1,
                       help='Define a weight for the conditional attribute. '
                            'The bigger the absolute value of the weight, the more important it will be. '
                            'If negative, the filter will look for smaller values.')
        group.addParam('overlapDef', params.EnumParam, label='Define overlap: ', default=0,
                       choices=self._overlapOptions,
                       help='How to define whether the input and the conditional ROIs overlap to count the conditional '
                            'score on the overlapping input ROI.\nContains: a ROI is totally contained into the other.'
                            '\nPartial overlap on X: the number of residues overlapping / len(X) is at least the '
                            'defined value')
        group.addParam('overlapSize', params.FloatParam, label='Minimum overlap: ', default=0.75,
                       condition='overlapDef>2',
                       help='Minimum proportion of overlap to be considered for partial overlap options.\n'
                            'e.g: Partial overlap on input ROI: 5 residues are overlapping for an input ROI of 8 '
                            'residues: overlap=5/8. If the value is over this threshold, we consider they overlap.')

        group.addParam('addCondFilter', params.LabelParam, label='Add conditional filter: ',
                       help='Add defined conditional filter')

        group = form.addGroup('Conditional filters summary')
        group.addParam('condSummary', params.TextParam, width=100, label='Conditional filters summary:', default='',
                       help='Summary of the selected conditional filters')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
      # 1) Scoring the input ROIs based on the filters defined
      scDic = self.getInputScores()

      # 2) Conditional scores for those input ROIs overlapping with the conditional input ROIs
      condScDic = self.getConditionalScores()

      scName = self.scoreName.get().strip() if self.scoreName.get().strip() else f'Score_{self.getObjId()}'
      condScName = self.scoreCondName.get().strip() if self.scoreCondName.get().strip() \
        else f'CondScore_{self.getObjId()}'

      outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
      for roi in self.inputSequenceROIs.get():
        curROI = roi.clone()
        roiId = curROI.getObjId()
        setattr(curROI, scName, Float(scDic[roiId]))
        setattr(curROI, condScName, Float(condScDic[roiId]))
        outROIs.append(curROI)

      if len(outROIs) > 0:
        self._defineOutputs(outputROIs=outROIs)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.filtSummary.get().strip():
          summary.append(f'Score summary: \n{self.filtSummary.get().strip()}')
          if self.condSummary.get().strip():
            summary.append(f'\nConditional score summary: \n{self.condSummary.get().strip()}')

        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getInputScores(self):
      '''Returns a dictionary of the form: {roiId: combinedScore}'''
      scDic = {}
      wDic = self.parseFilters()

      for roi in self.inputSequenceROIs.get():
        curROI = roi.clone()
        scDic[curROI.getObjId()] = 0
        for wKey, weight in wDic.items():
          scDic[curROI.getObjId()] += getattr(curROI, wKey).get() * weight
      return scDic

    def getConditionalScores(self):
      '''Returns a dictionary of the form: {roiId: combinedScore}'''
      condScDic = {}
      inROIsDic = {roi.clone().getObjId(): roi.clone() for roi in self.inputSequenceROIs.get()}
      cwDic = self.parseCondFilters()
      for roi in self.inputSequenceROIs.get():
        roiId = roi.clone().getObjId()
        condScDic[roiId] = 0
        roi = inROIsDic[roiId].clone()
        roiIdxs = roi.getROIIdxs()
        for (condAttr, condW) in cwDic:
          condROIs = cwDic[(condAttr, condW)]['ROIs']
          for cROI in condROIs:
            cROI = cROI.clone()
            condIdxs = cROI.getROIIdxs()
            if self.doOverlap(roiIdxs, condIdxs, cwDic[(condAttr, condW)]):
              condScDic[roi.getObjId()] += getattr(cROI, condAttr).get() * condW
      return condScDic

    def buildSumLine(self, params=None):
      return f'{{"Attribute": "{self.selectAttribute.get()}", "Weight": {self.attrWeight.get()}}}'

    def getCondSetId(self):
      return self.getCondSet().getObjId()

    def getCondSet(self, setId=None):
      setId = self.inSet.get() if setId is None else setId
      condSet = None
      for inPointer in self.conditionalROIs:
        inSet = inPointer.get()
        if inSet.__str__() == setId or inSet.getObjId() == setId:
          condSet = inSet
      return condSet

    def buildCondSumLine(self, params=None):
      sLine = f'{{"Set": {self.getCondSetId()}, "Attribute": "{self.condAttribute.get()}", ' \
             f'"Weight": {self.condWeight.get()}, "OverlapMode": {self.overlapDef.get()}'
      if self.overlapDef.get() > 2:
        sLine += f', "OverlapThreshold": {self.overlapSize.get()}'

      return sLine + '}'


    def parseFilters(self):
      '''Returs a dictionary of the form {attribute: weight} for the filter attributes
      '''
      wDic = {}
      for sumLine in self.filtSummary.get().split('\n'):
        sumLine = sumLine.strip()
        if sumLine:
          filterInfo = json.loads(sumLine[3:])
          attribute, weight = filterInfo["Attribute"], filterInfo["Weight"]
          wDic[attribute] = weight
      return wDic

    def parseCondFilters(self):
      '''Returns a dictionary of the form {(attribute, weight): [conditional rois]} for the conditional filter
      attributes
      '''
      wDic = {}
      for sumLine in self.condSummary.get().split('\n'):
        sumLine = sumLine.strip()
        if sumLine:
          filterInfo = json.loads(sumLine[3:])
          setId, attribute, weight = filterInfo["Set"], filterInfo["Attribute"], filterInfo["Weight"]
          ovMode = filterInfo["OverlapMode"]
          condSet = self.getCondSet(setId)
          wDic[(attribute, weight)] = {'ROIs': [roi.clone() for roi in condSet], 'ovMode': ovMode}
          if 'OverlapThreshold' in filterInfo:
            wDic[(attribute, weight)].update({'ovThres': filterInfo["OverlapThreshold"]})
          
      return wDic

    def doContain(self, idxs1, idxs2):
      '''Returns whether idxs1 contain idxs2'''
      return idxs1[0] <= idxs2[0] and idxs1[1] >= idxs2[1]

    def doPartialOverlap(self, idxs1, idxs2, overThres):
      '''Returns whether idxs1 and idxs2 overlap at least overThres proportion with respect to idxs1 length'''
      overNum = len(range(max(idxs1[0], idxs2[0]), min(idxs1[-1], idxs2[-1]) + 1))
      return overNum / len(range(idxs1[0], idxs1[1]+1)) >= overThres

    def doOverlap(self, idxs1, idxs2, overDic):
      overlap = False
      overMode = overDic['ovMode']
      if overMode == 0:
        overlap = self.doContain(idxs1, idxs2)
      elif overMode == 1:
        overlap = self.doContain(idxs2, idxs1)
      elif overMode == 2:
        overlap = self.doContain(idxs1, idxs2) or self.doContain(idxs2, idxs1)
      elif overMode == 3:
        overlap = self.doPartialOverlap(idxs1, idxs2, overDic['ovThres'])
      elif overMode == 4:
        overlap = self.doPartialOverlap(idxs2, idxs1, overDic['ovThres'])
      elif overMode == 5:
        overlap = self.doPartialOverlap(idxs2, idxs1, overDic['ovThres']) or \
                  self.doPartialOverlap(idxs1, idxs2, overDic['ovThres'])
      return overlap