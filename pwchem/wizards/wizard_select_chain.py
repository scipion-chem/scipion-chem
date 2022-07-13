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
from Bio.PDB import PDBParser, MMCIFParser

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pwem.wizards import EmWizard, SelectChainWizard, SelectResidueWizard, VariableWizard
from pwem.convert import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence

from pwchem.protocols import ProtDefineStructROIs, ProtChemPairWiseAlignment, ProtDefineSeqROI, ProtMapSequenceROI, \
  ProtDefineSetOfSequences
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



SelectChainWizardQT().addTarget(protocol=ProtMapSequenceROI,
                                targets=['chain_name'],
                                inputs=['inputAtomStruct', 'inputSequenceROIs'],
                                outputs=['chain_name'])


class PreviewAlignmentWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getPDBInputAS(self, AS, protocol):
      if str(type(AS).__name__) == 'SchrodingerAtomStruct':
        fileName = os.path.abspath(AS.convert2PDB())
      elif AS.getFileName().endswith('.pdbqt'):
        proteinFile = AS.getFileName()
        inName, inExt = os.path.splitext(os.path.basename(proteinFile))
        fileName = os.path.abspath(os.path.join(protocol.getProject().getPath(inName + '.pdb')))
        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), fileName)
        runOpenBabel(protocol=protocol, args=args, popen=True)

      else:
        fileName = os.path.abspath(AS.getFileName())
      return fileName

  def show(self, form, *params):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)

    model = json.loads(getattr(protocol, inputParam[0]).get())['model']
    chainID = json.loads(getattr(protocol, inputParam[0]).get())['chain']
    alignFile = os.path.abspath(protocol._getPath('preAlign_chain{}.fa'.format(chainID)))
    if not os.path.exists(alignFile):
        inASFile = self.getPDBInputAS(getattr(protocol, inputParam[1]).get(), form.protocol)
        inROI = getattr(protocol, inputParam[2]).get()
        inSeq = inROI.getSequence()
        handler = AtomicStructHandler(inASFile)
        seq = str(handler.getSequenceFromChain(modelID=model, chainID=chainID))
        pairwiseAlign(inSeq, seq, alignFile, force=True)

    SequenceAliView([alignFile], cwd=None).show()

PreviewAlignmentWizard().addTarget(protocol=ProtMapSequenceROI,
                                   targets=['preview'],
                                   inputs=['chain_name', 'inputAtomStruct', 'inputSequenceROIs'],
                                   outputs=[])


class SelectLigandWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def is_het(self, residue):
      res = residue.id[0]
      return res != " " and res != "W"

    def extract_ligands(self, ASFile):
        """ Extraction of the heteroatoms of .pdb files """
        molNames = []
        if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
            pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
            parser = PDBParser().get_structure(pdb_code, ASFile)
        elif ASFile.endswith('.cif'):
            pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
            parser = MMCIFParser().get_structure(pdb_code, ASFile)
        else:
            print('Unknown AtomStruct file format')
            return
        for model in parser:
            for chain in model:
                for residue in chain:
                    if self.is_het(residue):
                        if not residue.resname in molNames:
                            molNames.append(residue.resname)
        return molNames

    def show(self, form, *params):
      protocol = form.protocol
      inputParam, outputParam = self.getInputOutput(form)

      ASFile = getattr(protocol, inputParam[0]).get().getFileName()
      molNames = self.extract_ligands(ASFile)

      finalList = []
      for i in molNames:
        finalList.append(String(i))
      provider = ListTreeProviderString(finalList)
      dlg = dialog.ListDialog(form.root, "Ligand Names", provider,
                              "Select one of the ligands")
      form.setVar(outputParam[0], dlg.values[0].get())


SelectLigandWizard().addTarget(protocol=ProtDefineStructROIs,
                               targets=['molName'],
                               inputs=['inputAtomStruct'],
                               outputs=['molName'])


class SelectElementWizard(VariableWizard):
    """Lists the items in a SetOfX and choose one"""
    _targets, _inputs, _outputs = [], {}, {}

    def getListOfElements(self, protocol, inputParam):
      eleList = []
      scipionSet = getattr(protocol, inputParam[0]).get()
      if scipionSet is not None:
        for element in scipionSet:
            eleList.append(element.__str__())
      return eleList

    def show(self, form, *params):
      protocol = form.protocol
      inputParam, outputParam = self.getInputOutput(form)
      try:
        listOfElements = self.getListOfElements(protocol, inputParam)
      except Exception as e:
        print("ERROR: ", e)
        return

      finalList = []
      for i in listOfElements:
        finalList.append(String(i))
      provider = ListTreeProviderString(finalList)
      dlg = dialog.ListDialog(form.root, "Set items", provider,
                              "Select one of items in the set")
      form.setVar(outputParam[0], dlg.values[0].get())


SelectElementWizard().addTarget(protocol=ProtDefineStructROIs,
                               targets=['ligName'],
                               inputs=['inSmallMols'],
                               outputs=['ligName'])


SelectChainWizardQT().addTarget(protocol=ProtDefineSetOfSequences,
                                targets=['inpChain'],
                                inputs=[{'inputOrigin': ['inputSequence', #will not be displayed, but keep order of EnumParam
                                                         'inputAtomStruct', 'inputPDB']}],
                                outputs=['inpChain'])

SelectResidueWizardQT().addTarget(protocol=ProtDefineSetOfSequences,
                                  targets=['inpPositions'],
                                  inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct', 'inputPDB']},
                                          'inpChain'],
                                  outputs=['inpPositions'])


class AddSequenceWizard(SelectResidueWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParam = self.getInputOutput(form)

        # StructName
        inputObj = getattr(protocol, inputParams[0]).get()
        pdbFile, AS = '', False
        if issubclass(type(inputObj), str):
            outStr = [inputObj]
            AS = True
        elif issubclass(type(inputObj), AtomStruct):
            pdbFile = inputObj.getFileName()
            outStr = [os.path.splitext(os.path.basename(pdbFile))[0]]
            AS = True
        elif issubclass(type(inputObj), Sequence):
            outStr = [inputObj.getId()]

        if AS:
            # Chain
            chainJson = getattr(protocol, inputParams[1]).get()
            chainId = json.loads(chainJson)['chain']

        # Positions
        posJson = getattr(protocol, inputParams[2]).get()
        if posJson:
            posIdxs = json.loads(posJson)['index']
            seq = json.loads(posJson)['residues']
            outStr += [posIdxs]
        else:
            outStr += ['FIRST-LAST']
            finalResiduesList = self.getResidues(form, inputParams)
            idxs = [json.loads(finalResiduesList[0].get())['index'], json.loads(finalResiduesList[-1].get())['index']]
            seq = self.getSequence(finalResiduesList, idxs)

        chainStr = ''
        if AS:
          chainStr = ', "chain": "{}"'.format(chainId)

        prevStr = getattr(protocol, outputParam[0]).get()
        lenPrev = len(prevStr.strip().split('\n')) + 1
        if prevStr.strip() == '':
          lenPrev -= 1
        elif not prevStr.endswith('\n'):
          prevStr += '\n'

        seqFile = protocol.getProject().getTmpPath('{}_{}_{}.fa'.format(outStr[0], lenPrev, outStr[1]))
        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(outStr[0], seq))

        jsonStr = '%s) {"name": "%s"%s, "index": "%s", "seqFile": "%s"}\n' % \
                  (lenPrev, outStr[0], chainStr, outStr[1], seqFile)
        form.setVar(outputParam[0], prevStr + jsonStr)


AddSequenceWizard().addTarget(protocol=ProtDefineSetOfSequences,
                              targets=['addInput'],
                              inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct', 'inputPDB']},
                                      'inpChain', 'inpPositions'],
                              outputs=['inputList'])