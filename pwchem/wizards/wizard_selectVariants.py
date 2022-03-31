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

# Imports
from pwem.wizards.wizard import EmWizard
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog
from pyworkflow.object import String
from pwchem.protocols import ProtChemInsertVariants, ProtChemGenerateSequences


def getMutations(protocol):
  if hasattr(protocol, 'inputSequenceVariants'):
    vars = []
    inputVarObj = protocol.inputSequenceVariants.get()
    varFile = inputVarObj.getVariantsFileName()
    with open(varFile) as f:
      for line in f:
        vars.append(String(line.strip()))
    return vars

class SelectVariant(EmWizard):
  _targets = [(ProtChemInsertVariants, ['selectVariant']), (ProtChemGenerateSequences,['selectVariant'])]

  def getInputSeqVar(self, protocol):
      return protocol.inputSequenceVariants.get()

  def show(self, form, *params):
    protocol = form.protocol
    inpObj = self.getInputSeqVar(protocol)

    auxList = list(inpObj.getMutationsInLineage().keys())
    auxList.sort()
    varsList = []
    for var in auxList:
        varsList.append(String(var.strip()))

    provider = ListTreeProviderString(varsList)
    dlg = dialog.ListDialog(form.root, "Sequence Variants", provider,
                            "Select one Variant")
    form.setVar('selectVariant', dlg.values[0].get())

class SelectMutation(EmWizard):
  _targets = [(ProtChemInsertVariants, ['selectMutation'])]

  def show(self, form, *params):
    protocol = form.protocol
    mutList = getMutations(protocol)
    provider = ListTreeProviderString(mutList)
    dlg = dialog.ListDialog(form.root, "Sequence Mutations", provider,
                            "Select one Mutation (prevResidue-residueNumber-postResidue)")
    form.setVar('selectMutation', dlg.values[0].get())

class AddVarMutationWizard(EmWizard):
  _targets = [(ProtChemInsertVariants, ['addMutation'])]

  def show(self, form, *params):
    protocol = form.protocol
    if protocol.toMutateList.get() == None:
        written = ''
    else:
        written = protocol.toMutateList.get()
    form.setVar('toMutateList', written + '{}\n'.format(protocol.selectMutation.get()))

class AddVariantToListWizard(EmWizard):
  _targets = [(ProtChemGenerateSequences, ['addVariant'])]

  def show(self, form, *params):
    protocol = form.protocol
    if protocol.toMutateList.get() == None:
        written = ''
    else:
        written = protocol.toMutateList.get()
    form.setVar('toMutateList', written + '{}\n'.format(protocol.selectVariant.get()))
