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
    """Add a step of the workflow in the defined position"""
    #Variables _targets, _inputs and _outputs must be created in sons
    # _targets, _inputs, _outputs = [], {}, {}

    def addTarget(self, protocol, targets, inputs, outputs):
        self._targets += [(protocol, targets)]
        self._inputs.update({protocol: inputs})
        self._outputs.update({protocol: outputs})

    def getInputOutput(self, form):
        '''Retrieving input and output corresponding to the protocol where the wizard is used'''
        outParam = ''
        for target in self._targets:
            if form.protocol.__class__ == target[0]:
                inParam, outParam = self._inputs[target[0]], self._outputs[target[0]]
        return inParam, outParam
