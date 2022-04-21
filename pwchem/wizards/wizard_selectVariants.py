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

import numpy as np
import matplotlib.pyplot as plt

from pwem.wizards.wizard import EmWizard
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog
from pyworkflow.object import String
from pwchem.protocols import ProtChemGenerateVariants, ProtExtractSeqsROI, ProtDefineSeqROI
from pwchem.wizards import VariableWizard

################# Sequence variants ###############################
class SelectVariant(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)
    inpObj = getattr(protocol, inputParam[0]).get()

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
    inputParam, outputParam = self.getInputOutput(form)

    mutList = self.getMutations(protocol, inputParam[0])
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

    def show(self, form, *params):
      protocol = form.protocol
      inputParam, outputParam = self.getInputOutput(form)

      inList = getattr(protocol, inputParam[2]).get()
      num = len(inList.strip().split('\n'))
      if inList.strip() != '':
        num += 1

      addIdx = getattr(protocol, inputParam[0]).get()
      inParamName = inputParam[1][addIdx]
      label = protocol._originOptions[addIdx]

      roiInfo = getattr(protocol, inParamName).get()
      if len(inputParam) > 3 and hasattr(protocol, inputParam[3]):
          roiInfo = roiInfo.replace('}', ', "desc": "%s"}' % (getattr(protocol, inputParam[3]).get()))

      form.setVar(outputParam[0], getattr(protocol, inputParam[2]).get() + '{}) {}: {}\n'.format(num, label, roiInfo))


AddSequenceROIWizard().addTarget(protocol=ProtDefineSeqROI,
                                 targets=['addROI'],
                                 inputs=['whichToAdd', ['resPosition', 'selectVariant', 'selectMutation'],
                                         'inROIs', 'descrip'],
                                 outputs=['inROIs'])

AddSequenceROIWizard().addTarget(protocol=ProtChemGenerateVariants,
                                 targets=['addVariant'],
                                 inputs=['fromVariant', ['selectVariant', 'selectMutation'],
                                         'toMutateList'],
                                 outputs=['toMutateList'])




########################## Sequence conservation ####################################

class CheckSequencesConservation(EmWizard):
  _targets = [(ProtExtractSeqsROI, ['thres'])]

  def getConservationValues(self, protocol):
      protocol.calcConservation()
      with open(protocol.getConservationFile()) as f:
        consValues = f.readline().split()
      return list(map(float, consValues))

  def show(self, form, *params):
    protocol = form.protocol
    consValues = self.getConservationValues(protocol)

    xs = list(map(str, range(1, len(consValues) + 1)))
    y_pos = np.arange(len(xs))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(y_pos, consValues)
    yloc = plt.MaxNLocator(10)
    ax.yaxis.set_major_locator(yloc)
    ax.set_ylim(0, 1)

    plt.xlabel("Sequence position")
    plt.ylabel("Conservation value")
    plt.title('Conservation values along sequence')

    plt.show()

