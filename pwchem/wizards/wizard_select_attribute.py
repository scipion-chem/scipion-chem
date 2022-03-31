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

class SelectAttributeWizard(pwizard.Wizard):
  """Assist in the creation of python formula to be evaluated. In Steps"""
  _targets = [(chemprot.ProtocolConsensusDocking, ['action'])]

  def getInputAttributes(self, form):
    attrNames = []
    item = form.protocol.inputMoleculesSets[0].get().getFirstItem()
    for key, attr in item.getAttributesToStore():
      attrNames.append(key)
    return attrNames

  def show(self, form, *params):
    attrsList = self.getInputAttributes(form)
    finalAttrsList = []
    for i in attrsList:
      finalAttrsList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalAttrsList)

    dlg = dialog.ListDialog(form.root, "Filter set", provider,
                            "Select one of the attributes")
    form.setVar('action', dlg.values[0].get())