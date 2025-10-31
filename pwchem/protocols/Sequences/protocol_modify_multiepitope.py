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
This protocol is used to manually modify a multiepitope object
"""
import re

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import MultiEpitope, Sequence, SequenceROI
from pwchem.utils import unifyAttributes

REM, ADD, MOD = 0, 1, 2

class ProtModifyMultiEpitope(EMProtocol):
    """
    Allows to manually modify a MultiEpitope object
    """
    _label = 'Modify multiepitope'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputMultiEpitope', params.PointerParam,
                       pointerClass='MultiEpitope', allowsNull=False, label="Input multiepitope: ",
                       help='Select the MultiEpitope object to modify')

        group = form.addGroup('Add action')
        group.addParam('action', params.EnumParam, label='Action: ', default=0, choices=['Remove', 'Add', 'Modify'],
                       display=params.EnumParam.DISPLAY_HLIST,
                       help='Select an action to perform on the MultiEpitope over one of their epitopes/linkers')
        group.addParam('inROI', params.StringParam, label='Select a epitope/linker: ', default='',
                       condition=f'action in [{REM}, {MOD}]',
                       help='Choose a epitope/linker to remove/modify')
        group.addParam('addType', params.EnumParam, label='Add as: ', default=0, choices=['Epitope', 'Linker'],
                       condition=f'action in [{ADD}]', display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to add the defined sequence as epitope or linker')
        group.addParam('newROI', params.StringParam, label='Define a new epitope/linker sequence: ', default='',
                       condition=f'action in [{ADD}, {MOD}]',
                       help='Define a epitope/linker to add/modify')
        group.addParam('newName', params.StringParam, label='Define a new name: ', default='',
                       condition=f'action in [{ADD}, {MOD}]', expertLevel=params.LEVEL_ADVANCED,
                       help='Define a new name for the added/modified epitope/linker')
        group.addParam('addMod', params.LabelParam, label='Add modification: ',
                       help='Add the defined modification to the multiepitope')

        group = form.addGroup('Modifications summary')
        group.addParam('modSummary', params.TextParam, width=100, label='Modifications summary:', default='',
                       help='Summary of the modified elements in the multiepitope')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('defineOutputStep')

    def buildSumLine(self, extraPar=''):
      if self.action.get() == REM:
        sumLine = f'Remove {self.inROI.get()}'
      elif self.action.get() == ADD:
        newName = f'(Name {self.newName.get()})' if self.newName.get().strip() else ''
        sumLine = f'Add {self.getEnumText("addType")} {newName}: {self.newROI.get()}\n'
      elif self.action.get() == MOD:
        newName = f'(Name {self.newName.get()})' if self.newName.get().strip() else ''
        sumLine = f'Modify {self.inROI.get()} to {newName}: {self.newROI.get()}\n'
      return sumLine

    def defineOutputStep(self):
        actionsDic = self.parseSummary()

        outROIs, multiStr = [], ''
        for roi in self.inputMultiEpitope.get():
          if roi.getObjId() in actionsDic['Add']:
            roi = roi.clone()
            idxs = [len(multiStr) + 1, len(multiStr) + len(roi.getROISequence()) + 1]
            multiStr += roi.getROISequence()
            roi.setROIIdxs(idxs)
            outROIs.append(roi)
          elif roi.getObjId() in actionsDic['Modify']:
            roi = roi.clone()
            seqROIObj = actionsDic['Modify'][roi.getObjId()]
            if not seqROIObj.getId() and not seqROIObj.getSeqName():
              newName = f'{roi.getType()}_{len(outROIs)+1}'
              seqROIObj.setName(newName), seqROIObj.setId(newName)

            roiStr = seqROIObj.getSequence()
            idxs = [len(multiStr) + 1, len(multiStr) + len(roiStr) + 1]
            seqROI = SequenceROI(seqROI=seqROIObj, roiIdx=idxs[0], roiIdx2=idxs[1], type=roi.getType())
            multiStr += roiStr
            outROIs.append(seqROI)

        for seqROIObj in actionsDic['Add']:
          if not seqROIObj.getId() and not seqROIObj.getSeqName():
            newName = f'{seqROIObj.getDescription()}_{len(outROIs) + 1}'
            seqROIObj.setName(newName), seqROIObj.setId(newName)

          roiStr = seqROIObj.getSequence()
          idxs = [len(multiStr) + 1, len(multiStr) + len(roiStr) + 1]
          seqROI = SequenceROI(seqROI=seqROIObj, roiIdx=idxs[0], roiIdx2=idxs[1], type=seqROIObj.getDescription())
          multiStr += roiStr
          outROIs.append(seqROI)

        multiSeq = Sequence(sequence=multiStr, name='MultiEpitope', id='MultiEpitope',
                            description='MultiEpitope defined manually')

        outROIs = unifyAttributes(outROIs)
        if len(outROIs) > 0:
            outMultiEp = MultiEpitope(filename=self._getPath('multiEpitope.sqlite'))
            for i, roi in enumerate(outROIs):
                roi.setObjId(i)
                roi.setSequence(multiSeq)
                outMultiEp.append(roi)
            self._defineOutputs(outputROIs=outMultiEp)

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

    def _warnings(self):
        ws = []
        return ws

    # --------------------------- UTILS functions -----------------------------------

    def getInputROI(self, roiId=None):
        inROI = None
        inSet = self.inputMultiEpitope.get()
        roiId = self.inROI.get() if not roiId else roiId
        for item in inSet:
            if item.getObjId() == roiId or item.__str__() == roiId:
                inROI = item.clone()
        return inROI

    def parseSummary(self):
      actionsDic = {'Remove': [], 'Add': [], 'Modify': {}}
      for sumLine in self.modSummary.get().split('\n'):
        sumLine = sumLine.strip()[3:]
        if sumLine:
          if sumLine.startswith('Add'):
            addType = sumLine.split()[1].strip()
            nameSearch = re.search(r"\(Name (.*)\)", sumLine)
            newName = nameSearch.group(1) if nameSearch else ''
            newStr = sumLine.split(': ')[-1].strip()
            newSeq = Sequence(sequence=newStr, name=newName, id=newName, description=addType)
            actionsDic['Add'].append(newSeq)

          elif sumLine.startswith('Remove'):
            roiId = re.search(r"\(ID (\d+)\):", sumLine).group(1)
            actionsDic['Remove'].append(int(roiId))

          elif sumLine.startswith('Modify'):
            roiId = re.search(r"\(ID (\d+)\):", sumLine).group(1)

            nameSearch = re.search(r"\(Name (.*)\)", sumLine)
            newName = nameSearch.group(1) if nameSearch else ''
            newStr = sumLine.split(': ')[-1].strip()
            newSeq = Sequence(sequence=newStr, name=newName, id=newName)

            actionsDic['Modify'][int(roiId)] = newSeq
      return actionsDic



