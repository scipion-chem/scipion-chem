# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
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
This wizard will extract the chains from a atomic structure (pdb) file in
order to select it in the protocol.
Then, it will load the structure and will take all chain related
information such as name and number of residues.
"""

# Imports
import os, requests, json

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pwem.wizards import EmWizard
from pwem.convert import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence

from pwchem.protocols import ProtDefinePockets, ProtChemPairWiseAlignment, ProtDefineSeqROI
from pwchem.wizards.wizard_base import VariableWizard
from pwchem.utils import RESIDUES3TO1, RESIDUES1TO3

class AddResidueWizard(EmWizard):
    _targets = [(ProtDefinePockets, ['addResidue'])]

    def show(self, form, *params):
        protocol = form.protocol
        chainDic, resDic = json.loads(protocol.chain_name.get()), json.loads(protocol.resPosition.get())
        form.setVar('inResidues', protocol.inResidues.get() +
                    '{"model": %s, "chain": "%s", "index": "%s", "residues": "%s"}\n' %
                    (chainDic['model'], chainDic['chain'], resDic['index'], resDic['residues']))


######################## Variable wizards ####################

class SelectChainWizard(VariableWizard):
    '''Opens the input AtomStruct and allows you to select one of the present chains'''
    _targets, _inputs, _outputs = [], {}, {}

    @classmethod
    def getModelsChainsStep(cls, protocol, atomStructName='inputAtomStruct'):
      """ Returns (1) list with the information
         {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
         (2) list with residues, position and chain (modelsFirstResidue)"""
      structureHandler = AtomicStructHandler()
      fileName = ""
      AS = getattr(protocol, atomStructName).get()
      if AS is not None:
        fileName = os.path.abspath(AS.getFileName())
        if str(type(AS).__name__) == 'SchrodingerAtomStruct':
          fileName = os.path.abspath(AS.convert2PDB())
        else:
          fileName = os.path.abspath(AS.getFileName())

      structureHandler.read(fileName)
      structureHandler.getStructure()
      return structureHandler.getModelsChains()

    def editionListOfChains(self, listOfChains):
      self.chainList = []
      for model, chainDic in listOfChains.items():
        for chainID, lenResidues in chainDic.items():
          self.chainList.append(
            '{"model": %d, "chain": "%s", "residues": %d}' %
            (model, str(chainID), lenResidues))

    def show(self, form, *params):
      inputParam, outputParam = self.getInputOutput(form)
      protocol = form.protocol
      try:
        print('InputParam: ', inputParam[0])
        listOfChains, listOfResidues = self.getModelsChainsStep(protocol, inputParam[0])
      except Exception as e:
        print("ERROR: ", e)
        return

      self.editionListOfChains(listOfChains)
      finalChainList = []
      for i in self.chainList:
        finalChainList.append(String(i))
      provider = ListTreeProviderString(finalChainList)
      dlg = dialog.ListDialog(form.root, "Model chains", provider,
                              "Select one of the chains (model, chain, "
                              "number of chain residues)")
      form.setVar(outputParam[0], dlg.values[0].get())


class SelectChainWizard2(SelectChainWizard):
  _targets, _inputs, _outputs = [], {}, {}


SelectChainWizard().addTarget(protocol=ProtDefinePockets,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])

SelectChainWizard().addTarget(protocol=ProtChemPairWiseAlignment,
                              targets=['chain_name1'],
                              inputs=['inputAtomStruct1'],
                              outputs=['chain_name1'])

#There cannot be 2 VariableWizards of same type in the same protocol
#because we cannot know which target was clicked
SelectChainWizard2().addTarget(protocol=ProtChemPairWiseAlignment,
                               targets=['chain_name2'],
                               inputs=['inputAtomStruct2'],
                               outputs=['chain_name2'])


class SelectResidueWizard(SelectChainWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def editionListOfResidues(self, modelsFirstResidue, model, chain):
      residueList = []
      for modelID, chainDic in modelsFirstResidue.items():
        if int(model) == modelID:
          for chainID, seq_number in chainDic.items():
            if chain == chainID:
              for i in seq_number:
                residueList.append('{"index": %d, "residue": "%s"}' % (i[0], str(i[1])))
      return residueList

    def getResidues(self, form, inputParam):
      protocol = form.protocol
      inputObj = getattr(protocol, inputParam[0]).get()
      if issubclass(type(inputObj), AtomStruct):
          try:
            modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol, inputParam[0])
          except Exception as e:
            print("ERROR: ", e)
            return
          selection = getattr(protocol, inputParam[1]).get()

          struct = json.loads(selection)  # From wizard dictionary
          chain, model = struct["chain"].upper().strip(), int(struct["model"])

          residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
          finalResiduesList = []
          for i in residueList:
            finalResiduesList.append(String(i))

      elif issubclass(type(inputObj), Sequence):
          finalResiduesList = []
          for i, res in enumerate(inputObj.getSequence()):
            stri = '{"index": %s, "residue": "%s"}' % (i + 1, RESIDUES1TO3[res])
            finalResiduesList.append(String(stri))

      return finalResiduesList

    def show(self, form, *params):
      inputParam, outputParam = self.getInputOutput(form)
      finalResiduesList = self.getResidues(form, inputParam)
      provider = ListTreeProviderString(finalResiduesList)
      dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                              "Select one residue (residue number, "
                              "residue name)")
      roiStr = ''
      idxs = [json.loads(dlg.values[0].get())['index'], json.loads(dlg.values[-1].get())['index']]
      for i in range(idxs[0] - 1, idxs[1]):
        roiStr += RESIDUES3TO1[json.loads(finalResiduesList[i].get())['residue']]

      intervalStr = '{"index": "%s-%s", "residues": "%s"}' % (idxs[0], idxs[1], roiStr)
      form.setVar(outputParam[0], intervalStr)


SelectResidueWizard().addTarget(protocol=ProtDefinePockets,
                                targets=['resPosition'],
                                inputs=['inputAtomStruct', 'chain_name'],
                                outputs=['resPosition'])

SelectResidueWizard().addTarget(protocol=ProtDefineSeqROI,
                                targets=['resPosition'],
                                inputs=['inputSequence'],
                                outputs=['resPosition'])


class AddSequenceResidueWizard(EmWizard):
  _targets = [(ProtDefineSeqROI, ['addResidue'])]

  def show(self, form, *params):
    protocol = form.protocol
    form.setVar('inResidues', protocol.inResidues.get() + '{}\n'.format(protocol.resPosition.get()))