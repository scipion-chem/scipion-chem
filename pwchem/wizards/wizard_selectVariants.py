# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *
# * Biocomputing Unit, CNB-CSIC
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""

"""

import json
import numpy as np
import matplotlib.pyplot as plt

from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog
from pyworkflow.object import String

from pwem.wizards.wizard import VariableWizard, EmWizard
from pwem.objects import SetOfSequences, Pointer

from pwchem.protocols import ProtChemGenerateVariants, ProtDefineSeqROI

################# Sequence variants ###############################
class SelectVariant(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    protocol = form.protocol
    inputParams, outputParam = self.getInputOutput(form)
    inpObj = getattr(protocol, inputParams[0]).get()

    auxList = list(inpObj.getMutationsInLineage().keys())
    auxList.sort()
    varsList = [String('Original')]
    for var in auxList:
        varsList.append(String(var.strip()))

    provider = ListTreeProviderString(varsList)
    dlg = dialog.ListDialog(form.root, "Sequence Variants", provider,
                            "Select one Variant")
    form.setVar(outputParam[0], dlg.values[0].get())


SelectVariant().addTarget(protocol=ProtChemGenerateVariants,
                          targets=['selectVariant'],
                          inputs=['inputSequenceVariants'],
                          outputs=['selectVariant'])

SelectVariant().addTarget(protocol=ProtDefineSeqROI,
                          targets=['selectVariant'],
                          inputs=['inputSequenceVariants'],
                          outputs=['selectVariant'])


class SelectMutation(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getMutations(self, protocol, inputName):
    if hasattr(protocol, inputName):
      vars = []
      inputVarObj = getattr(protocol, inputName).get()
      varFile = inputVarObj.getVariantsFileName()
      with open(varFile) as f:
        for line in f:
          vars.append(String(line.split(';')[0].strip()))
      return vars

  def show(self, form, *params):
    protocol = form.protocol
    inputParams, outputParam = self.getInputOutput(form)

    mutList = self.getMutations(protocol, inputParams[0])
    provider = ListTreeProviderString(mutList)
    dlg = dialog.ListDialog(form.root, "Sequence Mutations", provider,
                            "Select one Mutation (prevResidue-residueNumber-postResidue)")
    muts = ''
    for mut in dlg.values:
        muts += mut.get().split()[0] + ', '

    form.setVar(outputParam[0], muts[:-2])


SelectMutation().addTarget(protocol=ProtChemGenerateVariants,
                           targets=['selectMutation'],
                           inputs=['inputSequenceVariants'],
                           outputs=['selectMutation'])

SelectMutation().addTarget(protocol=ProtDefineSeqROI,
                           targets=['selectMutation'],
                           inputs=['inputSequenceVariants'],
                           outputs=['selectMutation'])


class AddSequenceROIWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getOriginLabel(self, inParamName):
        if inParamName == 'resPosition':
            return 'Residues'
        elif inParamName == 'selectVariant':
            return 'Variant'
        elif inParamName == 'selectMutation':
            return 'Mutations'
        elif inParamName == 'inputSubsequences':
            return 'SubSequences'

    def getPrevPointersIds(self, prevPointers):
      ids = []
      for p in prevPointers:
        ids.append(p.get().getObjId())
      return ids

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParam = self.getInputOutput(form)

        inList = getattr(protocol, inputParams[1]).get()
        num = len(inList.strip().split('\n'))
        if inList.strip() != '':
          num += 1

        inParamName = inputParams[0]
        label = self.getOriginLabel(inParamName)

        if label == 'SubSequences':
            prevPointers = getattr(protocol, outputParam[1])
            prevIds = self.getPrevPointersIds(prevPointers)
            newSet = getattr(protocol, inputParams[0]).get()
            newId = newSet.getObjId()

            if not newId in prevIds:
                newIndex = len(prevPointers)
                prevPointers.append(Pointer(newSet))
            else:
                newIndex = prevIds.index(newId)
            form.setVar(outputParam[1], prevPointers)

        roiInfo = getattr(protocol, inParamName).get()
        if label == 'SubSequences':
            roiInfo = '{"PointerIdx": "%s", "Name": "%s"}' % (newIndex, roiInfo.__str__())
        elif len(inputParams) > 2 and hasattr(protocol, inputParams[2]):
            try:
                roiInfo = roiInfo.replace('}', ', "desc": "%s"}' % (getattr(protocol, inputParams[2]).get()))
            except:
                return

        form.setVar(outputParam[0], getattr(protocol, inputParams[1]).get() + '{}) {}: {}\n'.format(num, label, roiInfo))


AddSequenceROIWizard().addTarget(protocol=ProtDefineSeqROI,
                                 targets=['addROI'],
                                 inputs=[{'whichToAdd': ['resPosition', 'inputSubsequences',
                                                         'selectVariant', 'selectMutation']},
                                         'inROIs', 'descrip'],
                                 outputs=['inROIs', 'inputPointers'])

AddSequenceROIWizard().addTarget(protocol=ProtChemGenerateVariants,
                                 targets=['addVariant'],
                                 inputs=[{'fromVariant': ['selectVariant', 'selectMutation']},
                                         'toMutateList'],
                                 outputs=['toMutateList'])

