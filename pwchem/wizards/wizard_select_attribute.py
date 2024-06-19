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
import pyworkflow.object as pwobj
from pyworkflow.gui.tree import ListTreeProviderString
from pyworkflow.gui import dialog

from pwem.wizards import SelectAttributeWizard
from pwem.objects import String

from pwchem.wizards import VariableWizard
import pwchem.protocols as chemprot


class SelectAttributeWizardChem(SelectAttributeWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getInputPointer(self, multiPointer, pointerStr):
      for i, pointer in enumerate(multiPointer):
        if str(i).strip() == pointerStr.split('//')[0]:
          break
      return pointer

    def getInputSet(self, form, inputParam, inputStr=None):
      inputPointer = getattr(form.protocol, inputParam)
      if issubclass(inputPointer.__class__, pwobj.PointerList):
        inputPointer = inputPointer[0] if not inputStr else self.getInputPointer(inputPointer, inputStr)
      return inputPointer.get()

    def getFirstItem(self, form, inputParam, inputStr=None):
        inputSet = self.getInputSet(form, inputParam, inputStr)
        if issubclass(inputSet.__class__, pwobj.Set):
            item = inputSet.getFirstItem()
        elif issubclass(inputSet.__class__, pwobj.Object):
            item = inputSet
        return item

    def getInputAttributes(self, form, inputParam):
      attrNames = ['_objId']
      inputStr = getattr(form.protocol, inputParam[1]).get() if len(inputParam) > 1 else None
      item = self.getFirstItem(form, inputParam[0], inputStr)
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

SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtocolRankDocking,
                                      targets=['defineScore'],
                                      inputs=['inputMoleculesSets', 'defineInput'],
                                      outputs=['defineScore'])
SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtOptimizeMultiEpitopeGrape,
                                      targets=['inScore'],
                                      inputs=['inputROISets', 'inSet'],
                                      outputs=['inScore'])

for label in ['1', '2', '1_ratio', '2_ratio']:
  SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtScoreCorrelation,
                                        targets=[f'inputID_{label}'],
                                        inputs=[f'inputSet_{label}'],
                                        outputs=[f'inputID_{label}'])

  SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtScoreCorrelation,
                                        targets=[f'inputScore_{label}'],
                                        inputs=[f'inputSet_{label}'],
                                        outputs=[f'inputScore_{label}'])


class SelectMultiAttributeWizardChem(SelectAttributeWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    attrsList = self.getInputAttributes(form, inputParam)
    finalAttrsList = []
    for i in attrsList:
      finalAttrsList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalAttrsList)
    dlg = dialog.ListDialog(form.root, "Filter set", provider,
                            "Select one of the attributes")
    form.setVar(outputParam[0], ';'.join([val.get() for val in dlg.values]))

SelectMultiAttributeWizardChem().addTarget(protocol=chemprot.ProtChemOperateSet,
                                           targets=['remColumns'],
                                           inputs=['inputSet'],
                                           outputs=['remColumns'])


class SelectAttributeWizardListOperate(SelectAttributeWizardChem):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    operation = getattr(form.protocol, inputParam[0]).get()
    if not operation in [1, 2]:
      inputParam = [inputParam[1]]
    else:
      inputParam = [inputParam[2]]

    attrsList = self.getInputAttributes(form, inputParam)
    finalAttrsList = []
    for i in attrsList:
      finalAttrsList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalAttrsList)
    dlg = dialog.ListDialog(form.root, "Filter set", provider,
                            "Select one of the attributes")
    form.setVar(outputParam[0], dlg.values[0].get())

SelectAttributeWizardListOperate().addTarget(protocol=chemprot.ProtChemOperateSet,
                                      targets=['refColumn'],
                                      inputs=['operation', 'inputSet', 'inputMultiSet'],
                                      outputs=['refColumn'])
SelectAttributeWizardListOperate().addTarget(protocol=chemprot.ProtChemOperateSet,
                                      targets=['filterColumn'],
                                      inputs=['operation', 'inputSet', 'inputMultiSet'],
                                      outputs=['filterColumn'])
SelectAttributeWizardListOperate().addTarget(protocol=chemprot.ProtOperateSeqROI,
                                      targets=['bestAttribute'],
                                      inputs=['operation', 'inputROIsSets', 'inputMultiSet'],
                                      outputs=['bestAttribute'])


########################## Sequence Attributes ####################################

class SelectSequenceAttribute(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    protocol = form.protocol
    inputParams, outputParam = self.getInputOutput(form)
    inSeq = getattr(protocol, inputParams[0]).get()

    attrNames = list(map(String, inSeq.getAttributesDic()))

    provider = ListTreeProviderString(attrNames)
    dlg = dialog.ListDialog(form.root, "Sequence attributes", provider,
                            "Select one attribute")
    form.setVar(outputParam[0], dlg.values[0].get())


SelectSequenceAttribute().addTarget(protocol=chemprot.ProtExtractSeqsROI,
                                 targets=['inputAttribute'],
                                 inputs=['inputSequence'],
                                 outputs=['inputAttribute'])

class CheckSequencesAttribute(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getProtAttributeValues(self, protocol):
      return protocol.getAttributeValues()

  def show(self, form, *params):
    from pwchem.viewers.viewer_structure_attributes import plotSequenceAttribute
    protocol = form.protocol
    inputParams, outputParam = self.getInputOutput(form)
    thres = getattr(protocol, inputParams[0]).get()
    attrName = getattr(protocol, inputParams[1]).get()

    attrValues = self.getProtAttributeValues(protocol)
    plotSequenceAttribute(attrValues, attrName, thres)

CheckSequencesAttribute().addTarget(protocol=chemprot.ProtExtractSeqsROI,
                                    targets=['thres'],
                                    inputs=['thres', 'inputAttribute'],
                                    outputs=[])

class SelectMoleculesSubGroup(VariableWizard):
  """Select a molecules subgroup label """
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    molSet = getattr(protocol, inputParam[0]).get()
    vType = protocol.getEnumText(inputParam[1])
    vType = protocol.typeLabels[vType]

    ligDic = molSet.getGroupIndexes()[vType]
    outputLabels = list(ligDic.keys())

    finalOutputLabels = []
    for i in outputLabels:
      finalOutputLabels.append(pwobj.String(i))
    provider = ListTreeProviderString(finalOutputLabels)
    dlg = dialog.ListDialog(form.root, "Select molecule subgroup", provider,
                            "Select one of the subgroups")
    form.setVar(outputParam[0], dlg.values[0].get())


SelectMoleculesSubGroup().addTarget(protocol=chemprot.ProtDefineContactStructROIs,
                                    targets=['ligandSelection'],
                                    inputs=['inputSmallMols', 'selectionType'],
                                    outputs=['ligandSelection'])

class SelectEvaluationOrigin(VariableWizard):
  """Select a molecules subgroup label """
  _targets, _inputs, _outputs = [], {}, {}

  def importDDGProtocols(self):
    iedbProts = []
    try:
      from ddg.protocols import ProtDDGEvaluations
      iedbProts = ['DDG']
    except:
      print('No DDG protocols detected')
    return iedbProts

  def importIIITDProtocols(self):
    iedbProts = []
    try:
      from immuno.protocols import ProtIIITDEvaluations
      iedbProts = ['IIITD']
    except:
      print('No IIITD protocols detected')
    return iedbProts

  def importIEDBProtocols(self):
    iedbProts = []
    try:
      from iedb.protocols import ProtMHCIIPopulationCoverage
      iedbProts = ['IEDB']
    except:
      print('No IEDB protocols detected')
    return iedbProts

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)

    evalProts = []
    evalProts += self.importIEDBProtocols()
    evalProts += self.importIIITDProtocols()
    evalProts += self.importDDGProtocols()

    finalOutputLabels = []
    for i in evalProts:
      finalOutputLabels.append(pwobj.String(i))
    provider = ListTreeProviderString(finalOutputLabels)
    dlg = dialog.ListDialog(form.root, "Select evaluation origin", provider,
                            "Select one of the protocols")
    form.setVar(outputParam[0], dlg.values[0].get())


# SelectEvaluationOrigin().addTarget(protocol=chemprot.ProtOptimizeMultiEpitope,
#                                     targets=['multiEval'],
#                                     inputs=[],
#                                     outputs=['multiEval'])