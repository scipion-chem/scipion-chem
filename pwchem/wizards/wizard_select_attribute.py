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
from pwem.wizards import SelectAttributeWizard
import pwchem.protocols as chemprot


class SelectAttributeWizardChem(SelectAttributeWizard):
    def getInputAttributes(self, form, inputParam):
      attrNames = ['_objId']
      item = self.getFirstItem(form, inputParam[0])
      for key, attr in item.getAttributesToStore():
        attrNames.append(key)
      return attrNames

SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtAddAttribute,
                                      targets=['mapKey'],
                                      inputs=['inputSet'],
                                      outputs=['mapKey'])
SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtocolConsensusDocking,
                                      targets=['repAttr'],
                                      inputs=['inputMoleculesSets'],
                                      outputs=['repAttr'])
SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtocolScoreDocking,
                                      targets=['corrAttribute'],
                                      inputs=['inputMoleculesSets'],
                                      outputs=['corrAttribute'])



