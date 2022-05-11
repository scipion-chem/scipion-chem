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
The targets, inputs and outputs of this wizard can be easily accessed and modified
"""

# Imports
import pyworkflow.wizard as pwizard

class VariableWizard(pwizard.Wizard):
    """Wizard base class object where input and output paramNames can be modified and added, so
    one wizard can be used in several protocols with parameters of different names"""
    #Variables _targets, _inputs and _outputs must be created in sons
    # _targets, _inputs, _outputs = [], {}, {}

    def addTarget(self, protocol, targets, inputs, outputs):
        '''Add a target to a wizard and the input and output parameters are stored in a dictionary with
        (protocol, targetParamName) as key. Only the first target is used.'''
        self._targets += [(protocol, targets)]
        self._inputs.update({(protocol, targets[0]): inputs})
        self._outputs.update({(protocol, targets[0]): outputs})

    def getInputOutput(self, form):
        '''Retrieving input and output paramNames corresponding to the protocol and target of the wizard clicked'''
        outParam = ''
        for target in self._targets:
            if form.wizParamName == target[1][0] and form.protocol.__class__ == target[0]:
                prot, target = (target[0], target[1][0])
                inParam, outParam = self._inputs[(prot, target)], self._outputs[(prot, target)]
        return inParam, outParam
