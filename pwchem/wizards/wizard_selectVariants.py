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

