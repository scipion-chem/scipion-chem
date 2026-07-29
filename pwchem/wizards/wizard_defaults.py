# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
This wizards are used to set default values to some parameters, depending on other inputs
"""

# Imports
from pwem.wizards.wizard import VariableWizard
from pwchem.protocols import ProtocolPoseBusters

class SetDefaultValue(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, outputParam = self.getInputOutput(form)

        outputDefValue = getattr(protocol, inputParam[0])()

        form.setVar(outputParam[0], outputDefValue)


SetDefaultValue().addTarget(protocol=ProtocolPoseBusters,
                            targets=['testsPassed'],
                            inputs=['getNTests'],
                            outputs=['testsPassed'])

SetDefaultValue().addTarget(protocol=ProtocolPoseBusters,
                            targets=['specificTests'],
                            inputs=['getDefTests'],
                            outputs=['specificTests'])

SetDefaultValue().addTarget(protocol=ProtocolPoseBusters,
                            targets=['filterValue'],
                            inputs=['getFilterDefValue'],
                            outputs=['filterValue'])