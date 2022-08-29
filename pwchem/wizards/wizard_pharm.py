# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
Wizard srelated to pharmacophore protocols
"""

# Imports
from pwem.wizards import EmWizard, SelectChainWizard, SelectResidueWizard, VariableWizard

from pwchem.protocols import *
from pwchem.constants import *

######################## Variable wizards ####################

class WatchPharmacophoreAttributes(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getFeatureObj(self, inputPharm, featStr):
        featStr = featStr.split('|')[1].strip()
        for feat in inputPharm:
            if featStr == feat.__str__():
                return feat

    def getTypeIdx(self, inputFeat):
        typeStr = inputFeat.getType()
        choices = FEATURE_LABELS_SIMPLE + FEATURE_LABELS_ADVANCED
        return choices.index(typeStr)

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, outputParam = self.getInputOutput(form)

        inputFeatStr = getattr(protocol, inputParam[1]).get()
        pharmId = int(inputFeatStr.split('|')[0].split()[1])
        inputPharm = getattr(protocol, inputParam[0])[pharmId].get()

        inputFeat = self.getFeatureObj(inputPharm, inputFeatStr)

        form.setVar(outputParam[0], self.getTypeIdx(inputFeat))
        coords = inputFeat.getCoords()
        for i in range(3):
            form.setVar(outputParam[1+i], round(coords[i], 2))
        form.setVar(outputParam[4], round(inputFeat.getRadius(), 2))


WatchPharmacophoreAttributes().addTarget(protocol=ProtocolPharmacophoreModification,
                                         targets=['showCurrent'],
                                         inputs=['inputPharmacophores', 'currentFeatures'],
                                         outputs=['featType', 'featX', 'featY', 'featZ', 'featRadius'])


class AddPharmOperationWizard(VariableWizard):
  """Add an operation to perform on the pharmacophore"""
  _targets, _inputs, _outputs = [], {}, {}
  ADD, REM, MOD = 0, 1, 2

  def getPrevList(self, inputParam, protocol):
      prevList = getattr(protocol, inputParam[-1]).get()
      if not prevList:
        prevList = ''
      elif not prevList.endswith('\n'):
        prevList += '\n'
      return prevList

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    operation = getattr(protocol, inputParam[0]).get()
    curFeat = getattr(protocol, inputParam[1]).get()
    featType = protocol.getEnumText(inputParam[2])
    x, y, z = getattr(protocol, inputParam[3]).get(), getattr(protocol, inputParam[4]).get(), \
              getattr(protocol, inputParam[5]).get()
    radius = getattr(protocol, inputParam[6]).get()

    prevList = self.getPrevList(inputParam, protocol)

    if operation == self.ADD:
        addStr = 'ADD | {"Type": "%s", "Coords": "(%s, %s, %s)", "Radius": %s}' % (featType, x, y, z, radius)
        form.setVar(outputParam[0], prevList + '{}\n'.format(addStr.strip()))

    elif operation == self.REM:
        remStr = 'REMOVE | ' + curFeat
        form.setVar(outputParam[0], prevList + '{}\n'.format(remStr.strip()))

    elif operation == self.MOD:
      featId = ' '.join(curFeat.split()[:5])

      modStr = 'MODIFY | ' + featId + ' TO: ' + '{"Type": "%s", "Coords": "(%s, %s, %s)", "Radius": %s}' % \
               (featType, x, y, z, radius)
      form.setVar(outputParam[0], prevList + '{}\n'.format(modStr.strip()))


AddPharmOperationWizard().addTarget(protocol=ProtocolPharmacophoreModification,
                                    targets=['addOperation'],
                                    inputs=['operation', 'currentFeatures', 'featType',
                                            'featX', 'featY', 'featZ', 'featRadius',
                                            'operationList'],
                                    outputs=['operationList'])





























