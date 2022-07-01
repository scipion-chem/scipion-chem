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
from pwem.wizards import EmWizard, SelectChainWizard, SelectResidueWizard, VariableWizard
from pwem.convert import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence

from pwchem.protocols import ProtDefineStructROIs, ProtChemPairWiseAlignment, ProtDefineSeqROI, \
  ProtMapSequenceROI, ProtExtractSeqsROI
from pwchem.objects import SequenceVariants
from pwchem.viewers.viewers_sequences import SequenceAliViewer, SequenceAliView
from pwchem.utils import RESIDUES3TO1, RESIDUES1TO3, runOpenBabel
from pwchem.utils.utilsFasta import pairwiseAlign, calculateIdentity

class AddResidueWizard(EmWizard):
    _targets = [(ProtDefineStructROIs, ['addResidue'])]

    def show(self, form, *params):
        protocol = form.protocol
        chainDic, resDic = json.loads(protocol.chain_name.get()), json.loads(protocol.resPosition.get())
        form.setVar('inResidues', protocol.inResidues.get() +
                    '{"model": %s, "chain": "%s", "index": "%s", "residues": "%s"}\n' %
                    (chainDic['model'], chainDic['chain'], resDic['index'], resDic['residues']))


######################## Variable wizards ####################

class SelectChainWizardQT(SelectChainWizard):
  _targets, _inputs, _outputs = [], {}, {}
  @classmethod
  def getModelsChainsStep(cls, protocol, inputParamName):
    """ Returns (1) list with the information
       {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
       (2) list with residues, position and chain (modelsFirstResidue)"""
    structureHandler = AtomicStructHandler()
    fileName = ""
    AS = getattr(protocol, inputParamName).get()
    if type(AS) == str:
      if os.path.exists(AS):
        fileName = AS
      else:
        pdbID = AS
        url = "https://www.rcsb.org/structure/"
        URL = url + ("%s" % pdbID)
        try:
          response = requests.get(URL)
        except:
          raise Exception("Cannot connect to PDB server")
        if (response.status_code >= 400) and (response.status_code < 500):
          raise Exception("%s is a wrong PDB ID" % pdbID)
        fileName = structureHandler.readFromPDBDatabase(os.path.basename(pdbID), dir="/tmp/")

    elif str(type(AS).__name__) == 'SchrodingerAtomStruct':
        fileName = os.path.abspath(AS.convert2PDB())
    elif AS.getFileName().endswith('.pdbqt'):
      proteinFile = AS.getFileName()
      inName, inExt = os.path.splitext(os.path.basename(proteinFile))
      fileName = os.path.abspath(os.path.join(protocol.getProject().getPath(inName + '.pdb')))
      args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), fileName)
      runOpenBabel(protocol=protocol, args=args, popen=True)

    else:
        fileName = os.path.abspath(AS.getFileName())

    structureHandler.read(fileName)
    structureHandler.getStructure()
    return structureHandler.getModelsChains()

class SelectResidueWizardQT(SelectResidueWizard):
  _targets, _inputs, _outputs = [], {}, {}
  @classmethod
  def getModelsChainsStep(cls, protocol, inputParamName):
      return SelectChainWizardQT().getModelsChainsStep(protocol, inputParamName)

SelectChainWizardQT().addTarget(protocol=ProtDefineStructROIs,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtChemPairWiseAlignment,
                              targets=['chain_name1'],
                              inputs=['inputAtomStruct1'],
                              outputs=['chain_name1'])

SelectChainWizardQT().addTarget(protocol=ProtChemPairWiseAlignment,
                              targets=['chain_name2'],
                              inputs=['inputAtomStruct2'],
                              outputs=['chain_name2'])

SelectChainWizardQT().addTarget(protocol=ProtExtractSeqsROI,
                              targets=['chain_name'],
                              inputs=['inputAS'],
                              outputs=['chain_name'])


SelectResidueWizardQT().addTarget(protocol=ProtDefineStructROIs,
                                targets=['resPosition'],
                                inputs=['inputAtomStruct', 'chain_name'],
                                outputs=['resPosition'])

SelectResidueWizardQT().addTarget(protocol=ProtDefineSeqROI,
                                targets=['resPosition'],
                                inputs=[{'chooseInput': ['inputSequence', 'inputSequenceVariants']}],
                                outputs=['resPosition'])



class SelectChainPairwiseWizard(SelectChainWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def editionListOfChains(self, listOfChains, inputParam, protocol):
        chainList = []
        inAS = getattr(protocol, inputParam[0]).get()
        inROI = getattr(protocol, inputParam[1]).get()

        handler = AtomicStructHandler(inAS.getFileName())
        for model, chainDic in listOfChains.items():
          for chainID, lenResidues in chainDic.items():
            wDic = {"model": model, "chain": str(chainID), "residues": lenResidues}
            if inROI:
                inSeq = inROI.getSequence()
                seq = str(handler.getSequenceFromChain(modelID=model, chainID=chainID))

                alignFile = 'preAlign{}_chain{}.fa'.format(protocol.getObjId(), chainID)
                alignFile = os.path.abspath(protocol.getProject().getTmpPath(alignFile))
                if not os.path.exists(alignFile):
                    pairwiseAlign(inSeq, seq, alignFile, force=True)
                ident = calculateIdentity(alignFile)

                wDic["identity"] = ident

            chainList.append(str(wDic).replace("'", '"'))
        return chainList

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        try:
          listOfChains, listOfResidues = self.getModelsChainsStep(protocol, inputParam[0])
        except Exception as e:
          print("ERROR: ", e)
          return

        chainList = self.editionListOfChains(listOfChains, inputParam, protocol)
        finalChainList = []
        for i in chainList:
          finalChainList.append(String(i))
        provider = ListTreeProviderString(finalChainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains (model, chain, "
                                "number of chain residues)")
        form.setVar(outputParam[0], dlg.values[0].get())



SelectChainPairwiseWizard().addTarget(protocol=ProtMapSequenceROI,
                                      targets=['chain_name'],
                                      inputs=['inputAtomStruct', 'inputSequenceROIs'],
                                      outputs=['chain_name'])


class PreviewAlignmentWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)

    model = json.loads(getattr(protocol, inputParam[0]).get())['model']
    chainID = json.loads(getattr(protocol, inputParam[0]).get())['chain']
    alignFile = os.path.abspath(protocol._getPath('preAlign_chain{}.fa'.format(chainID)))
    if not os.path.exists(alignFile):
        inAS = getattr(protocol, inputParam[1]).get()
        inROI = getattr(protocol, inputParam[2]).get()
        inSeq = inROI.getSequence()
        handler = AtomicStructHandler(inAS.getFileName())
        seq = str(handler.getSequenceFromChain(modelID=model, chainID=chainID))
        pairwiseAlign(inSeq, seq, alignFile, force=True)

    SequenceAliView([alignFile], cwd=None).show()

PreviewAlignmentWizard().addTarget(protocol=ProtMapSequenceROI,
                                   targets=['preview'],
                                   inputs=['chain_name', 'inputAtomStruct', 'inputSequenceROIs'],
                                   outputs=[])
