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
from pwchem.protocols import ProtChemGenerateVariants, ProtExtractSeqsROI

################# Sequence variants ###############################
def getMutations(protocol, inputName='inputSequenceVariants'):
  if hasattr(protocol, inputName):
    vars = []
    inputVarObj = getattr(protocol, inputName).get()
    varFile = inputVarObj.getVariantsFileName()
    with open(varFile) as f:
      for line in f:
        vars.append(String(line.split(';')[0].strip()))
    return vars

class SelectVariant(EmWizard):
  _targets = [(ProtChemGenerateVariants,['selectVariant'])]

  def getInputSeqVar(self, protocol):
      return protocol.inputSequenceVariants.get()

  def show(self, form, *params):
    protocol = form.protocol
    inpObj = self.getInputSeqVar(protocol)

    auxList = list(inpObj.getMutationsInLineage().keys())
    auxList.sort()
    varsList = [String('Original')]
    for var in auxList:
        varsList.append(String(var.strip()))

    provider = ListTreeProviderString(varsList)
    dlg = dialog.ListDialog(form.root, "Sequence Variants", provider,
                            "Select one Variant")
    form.setVar('selectVariant', dlg.values[0].get())

class SelectMutation(EmWizard):
  _targets = [(ProtChemGenerateVariants,['selectMutation'])]

  def show(self, form, *params):
    protocol = form.protocol
    mutList = getMutations(protocol)
    provider = ListTreeProviderString(mutList)
    dlg = dialog.ListDialog(form.root, "Sequence Mutations", provider,
                            "Select one Mutation (prevResidue-residueNumber-postResidue)")
    muts = ''
    for mut in dlg.values:
        muts += mut.get().split()[0] + ', '

    form.setVar('selectMutation', muts[:-2])

class AddVariantToListWizard(EmWizard):
  _targets = [(ProtChemGenerateVariants, ['addVariant'])]

  def show(self, form, *params):
    protocol = form.protocol
    if protocol.toMutateList.get() == None or protocol.toMutateList.get().strip() == '':
        written = ''
        num = 1
    else:
        written = protocol.toMutateList.get()
        num = len(written.strip().split('\n')) + 1

    if protocol.fromVariant:
        form.setVar('toMutateList', written + '{}) Var: {}\n'.format(num, protocol.selectVariant.get()))
    else:
        form.setVar('toMutateList', written + '{}) Muts: {}\n'.format(num, protocol.selectMutation.get()))


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

