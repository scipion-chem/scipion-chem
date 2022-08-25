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
import pwchem.protocols as chemprot
from pwchem.wizards import VariableWizard

class AddElementSummaryWizard(VariableWizard):
    """Add a step of the workflow in the defined position"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        numSteps = protocol.countSteps()

        if getattr(protocol, inputParam[0]).get().strip() != '':
            index = int(getattr(protocol, inputParam[0]).get())
        else:
            index = numSteps + 1
        msjDic = protocol.createMSJDic()
        if index > numSteps:
            prevStr = getattr(protocol, outputParam[0]).get() \
                if getattr(protocol, outputParam[0]).get() is not None else ''
            form.setVar(outputParam[0], prevStr + str(msjDic) + '\n')

            newSum = protocol.createSummary()
            form.setVar(outputParam[1], newSum)

        elif numSteps >= index > 0:
            workSteps = getattr(protocol, outputParam[0]).get().split('\n')
            workSteps.insert(index-1, str(msjDic))
            form.setVar(outputParam[0], '\n'.join(workSteps))

            newSum = protocol.createSummary()
            form.setVar(outputParam[1], newSum)


class DeleteElementWizard(VariableWizard):
    """Delete the step of the workflow defined by the index"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        try:
            index = int(getattr(protocol, inputParam[0]).get().strip())
            if protocol.countSteps() >= index > 0:
                workSteps = getattr(protocol, outputParam[0]).get().split('\n')
                del workSteps[index - 1]
                workSteps.remove('')

                if workSteps != []:
                    form.setVar(outputParam[0], '\n'.join(workSteps) + '\n')
                else:
                    form.setVar(outputParam[0], '')
                newSum = protocol.createSummary()
                form.setVar(outputParam[1], newSum)
        except:
            print('Incorrect index')


AddElementWizard().addTarget(protocol=chemprot.ProtocolScoreDocking,
                             targets=['insertStep'],
                             inputs=['insertStep'],
                             outputs=['workFlowSteps', 'summarySteps'])

DeleteElementWizard().addTarget(protocol=chemprot.ProtocolScoreDocking,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])