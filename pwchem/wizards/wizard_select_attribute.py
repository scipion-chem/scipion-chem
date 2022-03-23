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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
from pyworkflow.gui.tree import ListTreeProviderString
import pyworkflow.wizard as pwizard
from pyworkflow.gui import dialog
import pyworkflow.object as pwobj

import pwchem.protocols as chemprot

class SelectAttributeWizardBase(pwizard.Wizard):
  """Base wizard for selecting an attribute from those contained in the items of a given input
  inputParam: Name of the input parameter where the items are stored
  outputParam:
  """
  _targets = []
  def getFirstItem(self, form, inputParam):
      inputPointer = getattr(form.protocol, inputParam)
      if issubclass(inputPointer.__class__, pwobj.PointerList):
          inputPointer = inputPointer[0]

      inputSet = inputPointer.get()
      if issubclass(inputSet.__class__, pwobj.Set):
          item = inputSet.getFirstItem()
      elif issubclass(inputSet.__class__, pwobj.Object):
          item = inputSet
      return item

  def getInputAttributes(self, form, inputParam):
    attrNames = []
    item = self.getFirstItem(form, inputParam)
    for key, attr in item.getAttributesToStore():
      attrNames.append(key)
    return attrNames

  def show(self, form, inputParam, outputParam, *params):
    attrsList = self.getInputAttributes(form, inputParam)
    finalAttrsList = []
    for i in attrsList:
      finalAttrsList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalAttrsList)
    dlg = dialog.ListDialog(form.root, "Filter set", provider,
                            "Select one of the attributes")
    form.setVar(outputParam, dlg.values[0].get())


class SelectAttributeWizard(SelectAttributeWizardBase):
  """Assist in the creation of python formula to be evaluated. In Steps"""
  _targets = [(chemprot.ProtocolConsensusDocking, ['action']),
              (chemprot.ProtocolScoreDocking, ['corrAttribute'])]
  _inputs = {chemprot.ProtocolConsensusDocking: 'inputMoleculesSets',
             chemprot.ProtocolScoreDocking: 'inputSmallMolecules'}

  def show(self, form, *params):
      outParam = ''
      # Retreiving the output parameter of the protocol used (the one that was clicked)
      # Determining the name of the input parameter whose attributes will be displayed
      for target in self._targets:
          if form.protocol.__class__ == target[0]:
              inParam = self._inputs[target[0]]
              outParam = target[1][0]   #Only first param in the list of targets

      if outParam:
          super().show(form, inputParam=inParam, outputParam=outParam, *params)
