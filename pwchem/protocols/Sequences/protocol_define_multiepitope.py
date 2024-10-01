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
This protocol is used to manually define a multiepitope from one or several SetOfROIs by choosing among these ROIs
and linkers
"""
import re

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import MultiEpitope, Sequence, SequenceROI
from pwchem.utils import unifyAttributes

SROI, MANU = 0, 1

class ProtDefineMultiEpitope(EMProtocol):
    """
    Allows to manually define a MultiEpitope object from several SetOfSequenceROIs
    """
    _label = 'Define multiepitope'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputROIsSets', params.MultiPointerParam,
                       pointerClass='SetOfSequenceROIs, MultiEpitope', allowsNull=False,
                       label="Input Sequence ROIs Sets",
                       help='Select the sequence ROIs sets to operate. If the input is only one set, '
                            'the operation will be performed over its own elements')

        group = form.addGroup('Add epitope/linker')
        group.addParam('origin', params.EnumParam, label='Define from: ', default=SROI,
                       choices=['SequenceROI', 'Manually'], display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to define the epitope/linker from a previously defined sequence ROI or manually')
        group.addParam('inSet', params.StringParam, label='From set: ', default='', condition='origin==0',
                       help='Add input sequence ROI as epitope from this set')
        group.addParam('inROI', params.StringParam, label='Sequence ROI as epitope: ', default='',
                       condition='origin==0', help='Choose a sequence ROI to add as epitope')

        group.addParam('manualROI', params.StringParam, label='Manual sequence: ', default='', condition='origin==1',
                       help='Define a epitope/linker sequence manually to add to the multiepitope')

        group.addParam('roiType', params.EnumParam, label='ROI type: ', default=0,
                       choices=['Epitope', 'Linker'], display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether the new ROI is defined as epitope or linker')
        group.addParam('roiName', params.StringParam, label='Name for the epitope: ', default='',
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Choose a name for the selected epitope')
        group.addParam('addROI', params.LabelParam, label='Add selected epitope: ',
                       help='Add the selected epitope to the multiepitope')

        group = form.addGroup('Multiepitope summary')
        group.addParam('multiSummary', params.TextParam, width=100, label='Multiepitope summary:', default='',
                       help='Summary of the added elements to build the multiepitope')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('defineOutputStep')
      
    def getNewROIName(self):
      if self.roiName.get().strip():
        roiName = self.roiName.get().strip()
      elif self.origin.get() == SROI:
        inROI = self.getInputROI()
        roiName = inROI.getROIId()
      else:
        roiName = ''
      return roiName

    def buildSumLine(self):
      '''Returns the information line to be written in the summary'''
      roiType = self.getEnumText('roiType')
      roiName = self.getNewROIName()

      if self.origin.get() == SROI:
        inSet, inROI = self.getInputSet(), self.getInputROI()
        sumLine = f'{roiType}: {inROI.__str__()} (Set {inSet.getObjId()}, Item {inROI.getObjId()})'
        if roiName:
          sumLine = sumLine[:-1] + f', Name {roiName})'

      else:
        sumLine = f'{roiType}: {self.manualROI.get().strip()}'
        if roiName:
          sumLine = sumLine + f' (Name {roiName})'

      return sumLine

    def defineOutputStep(self):
        outROIs = self.parseSummary()
        outROIs = unifyAttributes(outROIs)
        if len(outROIs) > 0:
            outMultiEp = MultiEpitope(filename=self._getPath('multiEpitope.sqlite'))
            for i, roi in enumerate(outROIs):
                roi.setObjId(i)
                outMultiEp.append(roi)
            self._defineOutputs(outputROIs=outMultiEp)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.multiSummary.get().strip():
          summary.append(self.multiSummary.get().strip())
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    def _warnings(self):
        ws = []
        return ws

    # --------------------------- UTILS functions -----------------------------------
    def getInputListOfSets(self):
      inSets = []
      for inPoint in self.inputROIsSets:
        inSets.append(inPoint.get())
      return inSets

    def getInputSet(self, setId=None):
        multiPointer = self.inputROIsSets
        inSet = None
        setId = self.inSet.get() if setId is None else setId
        for inPointer in multiPointer:
            curSet = inPointer.get()
            if curSet.getObjId() == setId or curSet.__str__() == setId:
                inSet = curSet
        return inSet

    def getInputROI(self, setId=None, roiId=None):
        inROI = None
        inSet = self.getInputSet(setId)
        roiId = self.inROI.get() if roiId is None else roiId
        for item in inSet:
            if item.getObjId() == roiId or item.__str__() == roiId:
                inROI = item.clone()
        return inROI

    def getElementName(self, sumLine):
      pat = re.compile(r"Name (.+)\)")
      search = re.search(pat, sumLine)
      name = search.group(1) if search else None
      return name

    def getROIIds(self, sumLine=None):
      '''Returns the ids of the set and the ROI in the summary line (if provided) or the selected in the form'''
      if sumLine:
        pat = re.compile(r"Set (\d*), Item (\d*)")
        search = re.search(pat, sumLine)
        setId, itemId = int(search.group(1)), int(search.group(2))
      else:
        inSet, inROI = self.getInputSet(), self.getInputROI()
        setId, itemId = inSet.getObjId(), inROI.getObjId()
      return int(setId), int(itemId)

    def parseSummary(self):
      epCount, lkCount = 0, 0
      outROIs, mSeq = [], ''
      for sumLine in self.multiSummary.get().split('\n'):
        sumLine = sumLine.strip()
        if sumLine:
          roiType = sumLine.split(') ')[1].split(':')[0]
          curCount = lkCount if roiType == 'Linker' else epCount
          curCount += 1

          if 'SequenceROI' in sumLine:
            # From SequenceROI
            setId, roiId = self.getROIIds(sumLine)
            roiName = self.getElementName(sumLine)
            seqROI = self.getInputROI(setId, roiId)
            roiSequence = seqROI.getROISequence()

            idxs = [len(mSeq) + 1, len(mSeq) + len(roiSequence) + 1]
            seqROI.setType(roiType)
            seqROI.setROIIdxs(idxs)
            if roiName:
              seqROI.setROIId(roiName), seqROI.setROIName(roiName)

          else:
            # Manual definition
            roiSequence = sumLine.split(f'{roiType}:')[-1].split('(')[0].strip()
            idxs = [len(mSeq)+1, len(mSeq) + len(roiSequence)+1]
            roiName = self.getElementName(sumLine) if self.getElementName(sumLine) else f'{roiType}_{curCount}'
            roiSeq = Sequence(sequence=roiSequence, name=roiName, id=roiName,
                              description='Manually added')
            seqROI = SequenceROI(seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1], type=roiType)

          mSeq += roiSequence
          outROIs.append(seqROI)

      multiSeq = Sequence(sequence=mSeq, name='MultiEpitope', id='MultiEpitope',
                          description='MultiEpitope defined manually')
      for roi in outROIs:
        roi.setSequence(multiSeq)

      return outROIs



