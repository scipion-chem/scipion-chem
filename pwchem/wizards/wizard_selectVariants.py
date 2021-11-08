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
from pwchem.protocols import InsertVariants

class SelectVariant(EmWizard):
  _targets = [(InsertVariants, ['selectVariant'])]    #Apply...: protocol object / 'selectVariant': StringParam for variant choose

  def getVariants(self, protocol):
    if hasattr(protocol, 'inputSequence'):  #'inputVariants': PointerParam for natural variants object
      vars = []
      inputVarObj = protocol.inputSequence.get()
      varFile = inputVarObj.getVariantsFileName()
      with open(varFile) as f:
        for line in f:
          vars.append(line.strip())
      return vars

  def show(self, form, *params):
    protocol = form.protocol
    varsList = self.getVariants(protocol)
    provider = ListTreeProviderString(varsList)
    dlg = dialog.ListDialog(form.root, "Sequence variants", provider,
                            "Select one variant (prevResidue-residueNumber-postResidue)")
    form.setVar('selectVariant', dlg.values[0].get())

class AddMutationWizard(EmWizard):
  _targets = [(InsertVariants, ['addVariant'])]   #addVariant: LabelParam: 'Add variant to list of mutations'

  def show(self, form, *params):
    protocol = form.protocol
    form.setVar('toMutateList', protocol.toMutateList.get() +   #toMutateList: TextParam: List of variants to apply
                '{}\n'.format(protocol.selectVariant.get()))


