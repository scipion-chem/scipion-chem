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
from pwchem.protocols import *
from pwchem.wizards import VariableWizard

class AddElementWizard(VariableWizard):
    """Add the content of a parameter to another"""
    _targets, _inputs, _outputs = [], {}, {}

    def curePrevList(self, prevList):
        if not prevList:
            prevList = ''
        elif not prevList.endswith('\n'):
            prevList += '\n'
        return prevList

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        inParam = getattr(protocol, inputParam[0]).get()
        if inParam and inParam.strip() != '':
            prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
            form.setVar(outputParam[0], prevList + '{}\n'.format(inParam.strip()))

class Add_FilterExpression(AddElementWizard):
    """Add ID or keyword in NCBI fetch protocol to the list"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        keep = protocol.getEnumText(inputParam[0])
        subset = protocol.getEnumText(inputParam[1])

        if subset and subset.strip() != '':
            prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
            towrite = prevList + '{} if in {}\n'.format(keep, subset)
            form.setVar(outputParam[0], towrite)


subGroups = list(ProtChemZINCFilter.subGroups.keys())
subChoices = ['subset_{}'.format(sb) for sb in subGroups]
Add_FilterExpression().addTarget(protocol=ProtChemZINCFilter,
                                 targets=['addFilter'],
                                 inputs=['mode', {'subGroup': subChoices}],
                                 outputs=['filterList'])

class AddLigandFilterExpression(AddElementWizard):
    """Add filter expression in ligand filter protocol"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        keepStr = protocol.getEnumText(inputParam[0])
        filterStr, fValue = protocol.getEnumText(inputParam[1]), getattr(protocol, inputParam[2]).get()

        prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
        towrite = prevList + '{} molecule if {} '.format(keepStr, filterStr.replace('x', str(fValue)).lower())
        if getattr(protocol, inputParam[1]).get() == 0:
            towrite += getattr(protocol, inputParam[3]).get()
        towrite += '\n'
        form.setVar(outputParam[0], towrite)

AddLigandFilterExpression().addTarget(protocol=ProtocolGeneralLigandFiltering,
                                      targets=['addFilter'],
                                      inputs=['mode', 'filter', 'filterValue', 'atomTypeFilter'],
                                      outputs=['filterList'])


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


AddElementSummaryWizard().addTarget(protocol=ProtocolScoreDocking,
                             targets=['insertStep'],
                             inputs=['insertStep'],
                             outputs=['workFlowSteps', 'summarySteps'])

DeleteElementWizard().addTarget(protocol=ProtocolScoreDocking,
                                targets=['deleteStep'],
                                inputs=['deleteStep'],
                                outputs=['workFlowSteps', 'summarySteps'])