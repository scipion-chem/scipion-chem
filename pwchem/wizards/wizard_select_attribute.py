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
from pwchem.viewers import SmallMoleculesLibraryViewer

SELECT_STR = "Select one of the attributes"

class SelectFromListWizard(VariableWizard):
    '''This wizard let the user select from a list that comes from a function i the target protocol,
    which is set in the input param'''
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParam, outputParam = self.getInputOutput(form)

        attrsList = getattr(protocol, inputParam[0])()
        finalAttrsList = []
        for i in attrsList:
          finalAttrsList.append(pwobj.String(i))
        provider = ListTreeProviderString(finalAttrsList)
        dlg = dialog.ListDialog(form.root, "Attribute selection", provider,
                                SELECT_STR)
        form.setVar(outputParam[0], dlg.values[0].get())

SelectFromListWizard().addTarget(protocol=chemprot.ProtMapAttributeToSeqROIs,
                                 targets=['attrName'],
                                 inputs=['getInputAttributes'],
                                 outputs=['attrName'])

SelectFromListWizard().addTarget(protocol=chemprot.ProtocolGeneralLigandFiltering,
                                 targets=['scoreFilter'],
                                 inputs=['getInputAttributes'],
                                 outputs=['scoreFilter'])

SelectFromListWizard().addTarget(protocol=chemprot.ProtocolOperateLibrary,
                                 targets=['filterAttr'],
                                 inputs=['getInputAttributes'],
                                 outputs=['filterAttr'])

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
SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtOptimizeMultiEpitope,
                                      targets=['inScore'],
                                      inputs=['inputROISets', 'inSet'],
                                      outputs=['inScore'])

class SelectMultiPointerAttributeWizard(VariableWizard):
  # todo: generalize in  SelectAttributeWizardChem
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    protocol = form.protocol
    _, outputParam = self.getInputOutput(form)

    attrsList = protocol.getAllInputScores()
    finalAttrsList = []
    for i in attrsList:
      finalAttrsList.append(pwobj.String(i))
    provider = ListTreeProviderString(finalAttrsList)
    dlg = dialog.ListDialog(form.root, "Attribute selection", provider,
                            SELECT_STR)
    form.setVar(outputParam[0], dlg.values[0].get())

SelectMultiPointerAttributeWizard().addTarget(protocol=chemprot.ProtOptimizeMultiEpitope,
                                              targets=['inScoreDef'],
                                              inputs=['inputROISets'],
                                              outputs=['inScoreDef'])

SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtocolRANXFuse,
                                      targets=['inAttrName'],
                                      inputs=['inputSets', 'inSetID'],
                                      outputs=['inAttrName'])
SelectAttributeWizardChem().addTarget(protocol=chemprot.ProtocolRANXFuse,
                                      targets=['inAttrVal'],
                                      inputs=['inputSets', 'inSetID'],
                                      outputs=['inAttrVal'])

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
                            SELECT_STR)
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
                            SELECT_STR)
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
      from iedb.protocols import ProtMHCPopulationCoverage
      iedbProts = ['IEDB']
    except:
      print('No IEDB protocols detected')
    return iedbProts

  def show(self, form, *params):
    _, outputParam = self.getInputOutput(form)

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

class SelectHeadersLibrary(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getInputHeaders(self, form):
    viewer = form.protocol
    return viewer.getLibrary().getHeaders()

  def show(self, form, *params):
    _, outputParam = self.getInputOutput(form)
    headers = self.getInputHeaders(form)

    finalOutputLabels = []
    for i in headers:
      finalOutputLabels.append(pwobj.String(i))
    provider = ListTreeProviderString(finalOutputLabels)
    dlg = dialog.ListDialog(form.root, "Select evaluation origin", provider,
                            "Select one of the protocols")
    form.setVar(outputParam[0], dlg.values[0].get())

for i in range(1, 3):
  SelectHeadersLibrary().addTarget(protocol=SmallMoleculesLibraryViewer,
                                   targets=[f'chooseHeader{i}'],
                                   inputs=[],
                                   outputs=[f'chooseHeader{i}'])

class SetParamValue(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getParamValue(self, form, inFunction, scoreIdx):
    viewer = form.protocol
    return getattr(viewer, inFunction)(scoreIdx)


  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)

    scoreIdx = outputParam[0][-1]
    paramValue = self.getParamValue(form, inputParam[0], scoreIdx)
    form.setVar(outputParam[0], paramValue)

for i in range(1, 3):
  SetParamValue().addTarget(protocol=SmallMoleculesLibraryViewer,
                            targets=[f'trueMin{i}'],
                            inputs=['getMinValue'],
                            outputs=[f'trueMin{i}'])

  SetParamValue().addTarget(protocol=SmallMoleculesLibraryViewer,
                            targets=[f'trueMax{i}'],
                            inputs=['getMaxValue'],
                            outputs=[f'trueMax{i}'])

class SelectZINCSubsetWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}


  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    prot = form.protocol
    zincSubsets = getattr(prot, inputParam[1])

    subset = prot.getEnumText(inputParam[0])
    hRange, logpRange = zincSubsets[subset][0], zincSubsets[subset][1]
    form.setVar(outputParam[0], hRange[0]), form.setVar(outputParam[1], hRange[1])
    form.setVar(outputParam[2], logpRange[0]), form.setVar(outputParam[3], logpRange[1])

SelectZINCSubsetWizard().addTarget(protocol=chemprot.ProtChemImportSmallMolecules,
                                   targets=['setRanges20'],
                                   inputs=['zinc20Subset', 'zinc20Subsets'],
                                   outputs=['minSize20', 'maxSize20', 'minLogP20', 'maxLogP20'])

SelectZINCSubsetWizard().addTarget(protocol=chemprot.ProtChemImportMoleculesLibrary,
                                   targets=['setRanges20'],
                                   inputs=['zinc20Subset', 'zinc20Subsets'],
                                   outputs=['minSize20', 'maxSize20', 'minLogP20', 'maxLogP20'])

SelectZINCSubsetWizard().addTarget(protocol=chemprot.ProtChemImportMoleculesLibrary,
                                   targets=['setRanges22'],
                                   inputs=['zinc22Subset', 'zinc22Subsets'],
                                   outputs=['minSize22', 'maxSize22', 'minLogP22', 'maxLogP22'])
