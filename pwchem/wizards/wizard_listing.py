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
import json

from pwchem.protocols import *
from pwchem.wizards import VariableWizard
from pwchem.utils import createMSJDic

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

        prevList, newLine = self.curePrevList(getattr(protocol, outputParam[0]).get()), ''
        if len(inputParam) > 0:
            inParam = getattr(protocol, inputParam[0]).get()
            if inParam and inParam.strip() != '':
                newLine = f'{inParam.strip()}\n'
        elif hasattr(protocol, 'createElementLine'):
            newLine = protocol.createElementLine()

        form.setVar(outputParam[0], prevList + newLine)

AddElementWizard().addTarget(protocol=ProtocolRankDocking,
                             targets=['defineSummary'],
                             inputs=[],
                             outputs=['defineSummary'])
AddElementWizard().addTarget(protocol=ProtocolRANXFuse,
                             targets=['addAttr'],
                             inputs=[],
                             outputs=['inAttrs'])
AddElementWizard().addTarget(protocol=ProtocolOperateLibrary,
                             targets=['addAttribute'],
                             inputs=['filterAttr'],
                             outputs=['removeList'])


class AddNumberedElementWizard(AddElementWizard):
    """Add the content of a parameter to another numbering (and labeling) the elements in the output
    Input[0]: name of the function in the protocol that builds the summary line
    Input[1]: optional, parameters to pass to this function
    """
    _targets, _inputs, _outputs = [], {}, {}

    def getNewElementNumber(self, prevList):
        return len(prevList.split('\n'))

    def getSumLine(self, protocol, inputParam, outputParam):
        prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
        elemNumber = self.getNewElementNumber(prevList)

        sumLine = ''
        if hasattr(protocol, inputParam[0]):
            buildLineFunc = getattr(protocol, inputParam[0])
            sumLineParams = inputParam[1] if len(inputParam) > 1 else []
            sumLine = buildLineFunc(*sumLineParams)

        if not sumLine:
            return ''
        return f'{elemNumber}) {sumLine}'

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        prevList = self.curePrevList(getattr(protocol, outputParam[0]).get())
        sumLine = self.getSumLine(protocol, inputParam, outputParam)
        if sumLine.strip():
            form.setVar(outputParam[0], prevList + sumLine.strip() + '\n')

AddNumberedElementWizard().addTarget(protocol=ProtDefineSeqROI,
                                      targets=['addROI'],
                                      inputs=['buildSumLine'],
                                      outputs=['inROIs'])

AddNumberedElementWizard().addTarget(protocol=ProtChemGenerateVariants,
                                 targets=['addVariant'],
                                 inputs=['buildSumLine'],
                                 outputs=['toMutateList'])

AddNumberedElementWizard().addTarget(protocol=ProtDefineMultiEpitope,
                                      targets=['addROI'],
                                      inputs=['buildSumLine'],
                                      outputs=['multiSummary'])

AddNumberedElementWizard().addTarget(protocol=ProtModifyMultiEpitope,
                                      targets=['addMod'],
                                      inputs=['buildSumLine'],
                                      outputs=['modSummary'])

AddNumberedElementWizard().addTarget(protocol=ProtCombineScoresSeqROI,
                                      targets=['addFilter'],
                                      inputs=['buildSumLine'],
                                      outputs=['filtSummary'])
AddNumberedElementWizard().addTarget(protocol=ProtCombineScoresSeqROI,
                                      targets=['addCondFilter'],
                                      inputs=['buildCondSumLine'],
                                      outputs=['condSummary'])

AddNumberedElementWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                      targets=['addScore'],
                                      inputs=['buildScoreSumLine'],
                                      outputs=['scoreSummary'])
AddNumberedElementWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                      targets=['addScoreDef'],
                                      inputs=['buildScoreSumLineDef'],
                                      outputs=['scoreSummaryDef'])

AddNumberedElementWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                      targets=['addLinker'],
                                      inputs=['buildLinkerSumLine'],
                                      outputs=['linkerSummary'])
AddNumberedElementWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                      targets=['addEval'],
                                      inputs=['buildEvalSumLine'],
                                      outputs=['evalSummary'])


AddElementWizard().addTarget(protocol=ProtChemZINCFilter,
                             targets=['addFilter'],
                             inputs=[],
                             outputs=['filterList'])

AddElementWizard().addTarget(protocol=ProtocolGeneralLigandFiltering,
                             targets=['addFilter'],
                             inputs=[],
                             outputs=['filterList'])

AddElementWizard().addTarget(protocol=ProtocolOperateLibrary,
                             targets=['addFilter'],
                             inputs=[],
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
        msjDic = createMSJDic(protocol)
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

class WatchElementWizard(VariableWizard):
    """Watch the parameters of the step of the workflow defined by the index"""
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        try:
            index = int(protocol.watchStep.get().strip())
            if protocol.countSteps() >= index > 0:
                workSteps = protocol.workFlowSteps.get().split('\n')
                msjDic = eval(workSteps[index - 1])
                for pName in msjDic:
                    if pName in protocol.getStageParamsDic(type='Normal').keys():
                        val = eval(msjDic[pName]) if msjDic[pName] in ['True', 'False'] else msjDic[pName]
                        form.setVar(pName, val)
                    elif pName in protocol.getStageParamsDic(type='Enum').keys():
                        enumParam = protocol.getParam(pName)
                        idx = enumParam.choices.index(msjDic[pName])
                        form.setVar(pName, idx)
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

class AddResidueWizard(AddElementWizard):
  """Adds a selected resiude(s) in a chain+residues congifuration parameters to a list of residues
  """
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    prevList, newLine = self.curePrevList(getattr(protocol, outputParam[0]).get()), ''
    if len(inputParam) > 0:
      inChain, inResidues = getattr(protocol, inputParam[0]).get(), getattr(protocol, inputParam[1]).get()
      chainDic, resDic = json.loads(inChain), json.loads(inResidues)
      newLine = f'{{"model": {chainDic["model"]}, "chain": "{chainDic["chain"]}", ' \
               f'"index": "{resDic["index"]}", "residues": "{resDic["residues"]}"}}\n'
    elif hasattr(protocol, 'createElementLine'):
      newLine = protocol.createElementLine()

    form.setVar(outputParam[0], prevList + newLine)