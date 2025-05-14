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
from Bio.PDB import PDBParser, MMCIFParser, PDBIO

from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pwem.wizards import SelectChainWizard, SelectResidueWizard, VariableWizard
from pwem.convert import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence, Pointer

from pwchem.protocols import *

from pwchem.objects import SequenceVariants, SetOfStructROIs, SetOfSmallMolecules
from pwchem.viewers.viewers_sequences import SequenceAliView
from pwchem.utils import RESIDUES3TO1, RESIDUES1TO3, runOpenBabel, natural_sort, parseAtomStruct
from pwchem.utils.utilsFasta import pairwiseAlign, calculateIdentity


class SelectLigandAtom(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def extract_atoms(self, molFile, protocol):
    """ Extraction of the atoms in a ligand file """
    if not molFile.endswith('.pdb'):
        proj = protocol.getProject()

        inName, inExt = os.path.splitext(os.path.basename(molFile))
        oFile = os.path.abspath(proj.getTmpPath(inName + '.pdb'))

        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(molFile), oFile)
        runOpenBabel(protocol=protocol, args=args, cwd=proj.getTmpPath(), popen=True)
        molFile = oFile

    parser = PDBParser().get_structure(molFile, molFile)
    atomNames = []
    for model in parser:
      for atom in model.get_atoms():
        atomNames.append(atom.get_name())
    atomNames = natural_sort(atomNames)
    return atomNames

  def getMol(self, inSet, molName):
    myMol = None
    for mol in inSet:
      if mol.__str__() == molName:
        myMol = mol.clone()
        break
    if myMol == None:
      print('The input ligand is not found')
      return None
    else:
      return myMol

  def show(self, form, *params):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)

    inSet, molName = getattr(protocol, inputParam[0]).get(), getattr(protocol, inputParam[1]).get()
    mol = self.getMol(inSet, molName)
    molFile = mol.getPoseFile() if mol.getPoseFile() else mol.getFileName()
    atomNames = self.extract_atoms(molFile, protocol)

    finalList = []
    for i in atomNames:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Ligand atom names", provider,
                            "Select one of the ligand atoms")
    form.setVar(outputParam[0], dlg.values[0].get())

class SelectLigandWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def isHet(self, residue):
    res = residue.id[0]
    return res != " " and res != "W"

  def extract_ligands(self, ASFile, protocol):
    """ Extraction of the heteroatoms of .pdb files """
    if ASFile.endswith('.pdbqt'):
        proj = protocol.getProject()
        oFile = os.path.abspath(proj.getTmpPath('inputStructure.pdb'))
        args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(ASFile), oFile)
        runOpenBabel(protocol=protocol, args=args, cwd=proj.getTmpPath(), popen=True)
        ASFile = oFile

    if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
      pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
      parser = PDBParser().get_structure(pdb_code, ASFile)
    elif ASFile.endswith('.cif'):
      pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
      parser = MMCIFParser().get_structure(pdb_code, ASFile)
    else:
      print('Unknown AtomStruct file format')
      return

    molNames = []
    for model in parser:
      for chain in model:
        for residue in chain:
          if self.isHet(residue):
            if not residue.resname in molNames:
              molNames.append(residue.resname)
    return molNames

  def show(self, form, *params):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)

    ASFile = getattr(protocol, inputParam[0]).get().getFileName()
    molNames = self.extract_ligands(ASFile, protocol)

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

SelectLigandWizard().addTarget(protocol=ProtocolRMSDDocking,
                               targets=['refLigName'],
                               inputs=['refAtomStruct'],
                               outputs=['refLigName'])

class SelectMultiLigandWizard(SelectLigandWizard):
  def show(self, form, *params):
      protocol = form.protocol
      inputParam, outputParam = self.getInputOutput(form)

      ASFile = getattr(protocol, inputParam[0]).get().getFileName()
      molNames = self.extract_ligands(ASFile, protocol)

      finalList = []
      for i in molNames:
        finalList.append(String(i))
      provider = ListTreeProviderString(finalList)
      dlg = dialog.ListDialog(form.root, "Ligand Names", provider,
                              "Select one of the ligands")

      vals = [v.get() for v in dlg.values]
      form.setVar(outputParam[0], ', '.join(vals))

SelectMultiLigandWizard().addTarget(protocol=ProtChemPrepareReceptor,
                               targets=['het2keep'],
                               inputs=['inputAtomStruct'],
                               outputs=['het2keep'])


class AddROIWizard(VariableWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def getPrevList(self, protocol, outputParam):
    prevList = getattr(protocol, outputParam[0]).get()
    if not prevList or not prevList.strip():
      prevList = ''
    elif not prevList.endswith('\n'):
      prevList += '\n'
    lenPrev = len(prevList.split('\n'))
    return prevList, lenPrev

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    prevList, lenPrev = self.getPrevList(protocol, outputParam)
    roiDef = protocol.getDefinedROILine()
    if protocol.origin.get() == 2 and protocol.extLig.get():
      _, newPointers = protocol.getNewPointers()
      form.setVar(outputParam[1], newPointers)

    form.setVar(outputParam[0], f'{prevList}{lenPrev}) {roiDef}')

AddROIWizard().addTarget(protocol=ProtDefineStructROIs,
                         targets=['addROI'],
                         inputs=[],
                         outputs=['inROIs', 'inputPointers'])



######################## Variable wizards ####################

class SelectChainWizardQT(SelectChainWizard):
  _targets, _inputs, _outputs = [], {}, {}

  @classmethod
  def getInputFilename(cls, protocol, inputObj, structureHandler):
    if type(inputObj) == str:
      if os.path.exists(inputObj):
        fileName = inputObj
      else:
        pdbID = inputObj
        url = "https://www.rcsb.org/structure/"
        URL = url + ("%s" % pdbID)
        try:
          response = requests.get(URL)
        except:
          raise Exception("Cannot connect to PDB server")
        if (response.status_code >= 400) and (response.status_code < 500):
          raise Exception("%s is a wrong PDB ID" % pdbID)
        fileName = structureHandler.readFromPDBDatabase(os.path.basename(pdbID), dir="/tmp/")

    elif str(type(inputObj).__name__) == 'SchrodingerAtomStruct':
        fileName = os.path.abspath(inputObj.convert2PDB())

    elif str(type(inputObj).__name__) == 'SetOfStructROIs' or str(type(inputObj).__name__) == 'SetOfSmallMolecules':
        fileName = inputObj.getProteinFile()

    else:
        fileName = os.path.abspath(inputObj.getFileName())

    if fileName.endswith('.pdbqt'):
      inName, inExt = os.path.splitext(os.path.basename(fileName))
      pdbFile = os.path.abspath(os.path.join(protocol.getProject().getPath(inName + '.pdb')))
      args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(fileName), pdbFile)
      runOpenBabel(protocol=protocol, args=args, popen=True)
      fileName = pdbFile
    return fileName

  @classmethod
  def getModelsChainsStep(cls, protocol, inputObj):
    """ Returns (1) list with the information
       {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
       (2) list with residues, position and chain (modelsFirstResidue)"""
    structureHandler = AtomicStructHandler()
    fileName = cls.getInputFilename(protocol, inputObj, structureHandler)
    structureHandler.read(fileName)
    structureHandler.getStructure()
    return structureHandler.getModelsChains()

class SelectResidueWizardQT(SelectResidueWizard, SelectChainWizardQT):
  _targets, _inputs, _outputs = [], {}, {}
  @classmethod
  def getModelsChainsStep(cls, protocol, inputObj):
      return SelectChainWizardQT().getModelsChainsStep(protocol, inputObj)

  def checkNoPDBQT(self, form, fileName):
      if fileName.endswith('.pdbqt'):
        protocol = form.protocol
        inName, inExt = os.path.splitext(os.path.basename(fileName))
        pdbFile = os.path.abspath(os.path.join(protocol.getProject().getPath(inName + '.pdb')))
        args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(fileName), pdbFile)
        runOpenBabel(protocol=protocol, args=args, popen=True)
        fileName = pdbFile
      return fileName

  def getResidues(self, form, inputObj, chainStr):
    if issubclass(type(inputObj), Sequence) or str(type(inputObj).__name__) == 'SequenceVariants':
      finalResiduesList = []
      for i, res in enumerate(inputObj.getSequence()):
        if res in RESIDUES1TO3:
          res3 = RESIDUES1TO3[res]
        else:
          res3 = res
        stri = '{"index": %s, "residue": "%s"}' % (i + 1, res3)
        finalResiduesList.append(String(stri))

    else:
      structureHandler = AtomicStructHandler()
      fileName = self.getInputFilename(form.protocol, inputObj, structureHandler)

      structureHandler.read(fileName)
      structureHandler.getStructure()
      _, modelsFirstResidue = structureHandler.getModelsChains()

      struct = json.loads(chainStr)  # From wizard dictionary
      chain, model = struct["chain"].upper().strip(), int(struct["model"])

      residueList = self.editionListOfResidues(modelsFirstResidue, model, chain)
      finalResiduesList = []
      for i in residueList:
        finalResiduesList.append(String(i))

    return finalResiduesList


class SelectAtomWizardQT(SelectResidueWizardQT):
  _targets, _inputs, _outputs = [], {}, {}

  def getAtoms(self, form, inputObj, chainStr, resStr):
    protocol = form.protocol
    structureHandler = AtomicStructHandler()
    parser = parseAtomStruct(self.getInputFilename(protocol, inputObj, structureHandler))

    chainJson, residueJson = json.loads(chainStr), json.loads(resStr)  # From wizard dictionary
    chainID, modelID = chainJson["chain"].upper().strip(), int(chainJson["model"])
    residueID = int(residueJson['index'].split('-')[0])

    atomList = self.editionListOfAtoms(parser, modelID, chainID, residueID)
    finalAtomList = []
    for i in atomList:
      finalAtomList.append(String(i))
    return finalAtomList

  def editionListOfAtoms(self, parser, modelID, chainID, residueID):
    atomList = []
    for model in parser:
      if model.get_id() == modelID:
        for chain in model.get_chains():
          if chain.get_id() in chainID:
            for residue in chain.get_residues():
              if residue.get_id()[1] == residueID:
                for i, atom in enumerate(residue.get_atoms()):
                  atomID = atom.get_id()
                  atomList.append(f'{{"index": {i+1}, "atom": "{atomID}"}}')

    return atomList

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    inputObj = getattr(protocol, inputParams[0]).get()
    chainStr = getattr(protocol, inputParams[1]).get()
    residueStr = getattr(protocol, inputParams[2]).get()

    finalAtomsList = self.getAtoms(form, inputObj, chainStr, residueStr)

    provider = ListTreeProviderString(finalAtomsList)
    dlg = dialog.ListDialog(form.root, "Residue atoms", provider,
                            "Select one atom (atom number, "
                            "atom name)")

    idx, atomID = json.loads(dlg.values[0].get())['index'], json.loads(dlg.values[0].get())['atom']

    intervalStr = '{"index": "%s", "atom": "%s"}' % (idx, atomID)
    form.setVar(outputParam[0], intervalStr)



SelectChainWizardQT().addTarget(protocol=ProtDefineStructROIs,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])
SelectChainWizardQT().addTarget(protocol=ProtDefineStructROIs,
                              targets=['chain_name2'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name2'])

SelectChainWizardQT().addTarget(protocol=ProtChemPairWiseAlignment,
                              targets=['chain_name1'],
                              inputs=['inputAtomStruct1'],
                              outputs=['chain_name1'])

SelectChainWizardQT().addTarget(protocol=ProtChemPairWiseAlignment,
                              targets=['chain_name2'],
                              inputs=['inputAtomStruct2'],
                              outputs=['chain_name2'])

SelectChainWizardQT().addTarget(protocol=ProtSeqCalculateConservation,
                              targets=['chain_name'],
                              inputs=['inputAS'],
                              outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtExtractLigands,
                              targets=['chain_name'],
                              inputs=['inputStructure'],
                              outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtCalculateSASA,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtCalculateSASA,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
                              outputs=['chain_name'])

SelectChainWizardQT().addTarget(protocol=ProtMapAttributeToSeqROIs,
                              targets=['chain_name'],
                              inputs=['inputAtomStruct'],
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

    def editionListOfChains(self, listOfChains, inAS, inROI, protocol):
        chainList = []

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
        inputObj = getattr(protocol, inputParam[0]).get()
        try:
            listOfChains, listOfResidues = self.getModelsChainsStep(protocol, inputObj)
        except Exception as e:
            print("ERROR: ", e)
            return

        inROI = getattr(protocol, inputParam[1]).get()
        chainList = self.editionListOfChains(listOfChains, inputObj, inROI, protocol)
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


class SelectMultiChainWizard(SelectChainWizardQT):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
      inputParams, outputParam = self.getInputOutput(form)
      protocol = form.protocol
      inputObj = getattr(protocol, inputParams[0]).get()
      try:
          listOfChains, listOfResidues = self.getModelsChainsStep(protocol, inputObj)
      except Exception as e:
          print("ERROR: ", e)
          return

      chainList = self.editionListOfChains(listOfChains)
      finalChainList = []
      for i in chainList:
        finalChainList.append(String(i))
      provider = ListTreeProviderString(finalChainList)
      dlg = dialog.ListDialog(form.root, "Model chains", provider,
                              "Select one of the chains (model, chain, "
                              "number of chain residues)")
      if len(dlg.values) > 1:
        chains = []
        for selChain in dlg.values:
          chains += ['{}-{}'.format(json.loads(selChain.get())['model'], json.loads(selChain.get())['chain'])]
        chainInfo = '{"model-chain": "%s"}' % (', '.join(chains))
      else:
        chainInfo = dlg.values[0].get()
      form.setVar(outputParam[0], chainInfo)

SelectMultiChainWizard().addTarget(protocol=ProtDefineStructROIs,
                                   targets=['keep_chain_name'],
                                   inputs=['inputAtomStruct'],
                                   outputs=['keep_chain_name'])

SelectMultiChainWizard().addTarget(protocol=ProtChemPrepareReceptor,
                                   targets=['chain_name'],
                                   inputs=['inputAtomStruct'],
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
    inputAS = getattr(protocol, inputParam[1]).get()

    alignFile = os.path.abspath(protocol._getPath('preAlign_chain{}.fa'.format(chainID)))
    if not os.path.exists(alignFile):
        inASFile = self.getPDBInputAS(inputAS, form.protocol)
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


class SelectElementWizard(VariableWizard):
  """Lists the items in a SetOfX and choose one"""
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, protocol, scipionSet):
    eleList = []
    if scipionSet is not None:
      for element in scipionSet:
          eleList.append(element.__str__())
    return eleList

  def displayDialog(self, form, inputParam):
    protocol = form.protocol
    try:
      scipionSet = getattr(protocol, inputParam[0]).get()
      listOfElements = self.getListOfElements(protocol, scipionSet)
    except Exception as e:
      print("ERROR: ", e)
      return

    finalList = []
    for i in listOfElements:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Set items", provider,
                            "Select one of items in the set")
    return dlg

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    dlg = self.displayDialog(form, inputParam)
    form.setVar(outputParam[0], dlg.values[0].get())
    
SelectElementWizard().addTarget(protocol=ProtocolLigandParametrization,
                               targets=['inputLigand'],
                               inputs=['inputSmallMolecules'],
                               outputs=['inputLigand'])


class SelectElementMultiPointWizard(SelectElementWizard):
  """Lists the items in a SetOfX selected from a multipointer"""
  _targets, _inputs, _outputs = [], {}, {}

  def getInputSet(self, multiPointer, setStr):
    inSet = None
    for inPointer in multiPointer:
      curSet = inPointer.get()
      if curSet.__str__() == setStr:
        inSet = curSet
    return inSet

  def displayDialog(self, form, inputParam):
    protocol = form.protocol
    try:
      multiPointer = getattr(protocol, inputParam[0])
      inputSetStr = getattr(protocol, inputParam[1]).get()

      inSet = self.getInputSet(multiPointer, inputSetStr)
      listOfElements = self.getListOfElements(protocol, inSet)
    except Exception as e:
      print("ERROR: ", e)
      return

    finalList = []
    for i in listOfElements:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Set items", provider,
                            "Select one of items in the set")
    return dlg

SelectElementMultiPointWizard().addTarget(protocol=ProtDefineMultiEpitope,
                                 targets=['inROI'],
                                 inputs=['inputROIsSets', 'inSet'],
                                 outputs=['inROI'])

class SelectInputSetWizard(VariableWizard):
  '''Select a set form a MultiPointer param'''
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, multipointer):
    eleList = []
    if multipointer is not None:
      for inPointer in multipointer:
        eleList.append(inPointer.get())
    return eleList

  def displayDialog(self, form, inputParam):
    protocol = form.protocol
    try:
      multiPointer = getattr(protocol, inputParam[0])
      listOfElements = self.getListOfElements(multiPointer)
    except Exception as e:
      print("ERROR: ", e)
      return

    finalList = []
    for i in listOfElements:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "Select set", provider,
                            "Select one of the sets in the input")
    return dlg

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    dlg = self.displayDialog(form, inputParam)
    form.setVar(outputParam[0], dlg.values[0].get())

SelectInputSetWizard().addTarget(protocol=ProtDefineMultiEpitope,
                                 targets=['inSet'],
                                 inputs=['inputROIsSets'],
                                 outputs=['inSet'])
SelectInputSetWizard().addTarget(protocol=ProtCombineScoresSeqROI,
                                 targets=['inSet'],
                                 inputs=['conditionalROIs'],
                                 outputs=['inSet'])


class SelectMultiElementWizard(SelectElementWizard):
  """Lists the items in a SetOfX and choose one or several"""
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    dlg = self.displayDialog(form, inputParam)
    values = [val.get().strip() for val in dlg.values]
    form.setVar(outputParam[0], ', '.join(values))

SelectMultiElementWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                     targets=['seleLinker'],
                                     inputs=['inLinkerSet'],
                                     outputs=['seleLinker'])

class SelectMultiSeqWizard(SelectMultiElementWizard):
  """Lists the items in a SetOfSequences and choose several"""
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, protocol, scipionSet):
    eleList = []
    if scipionSet is not None:
      for element in scipionSet:
        eleList.append(element.getSeqName())
    return ['All'] + eleList

class SelectMultiMolWizard(SelectMultiElementWizard):
  """Lists the interacting mols in a SetOfSequences and choose several"""
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, protocol, seqSet):
    return ['All'] + seqSet.getInteractMolNames()

SelectMultiSeqWizard().addTarget(protocol=ProtExtractInteractingMols,
                                 targets=['chooseSeq'],
                                 inputs=['inputSequences'],
                                 outputs=['chooseSeq'])

SelectMultiMolWizard().addTarget(protocol=ProtExtractInteractingMols,
                                 targets=['chooseMol'],
                                 inputs=['inputSequences'],
                                 outputs=['chooseMol'])

class SelectMultiEpitopeElementWizard(SelectElementWizard):
  """Lists the items in a MultiEpitope and choose several"""
  _targets, _inputs, _outputs = [], {}, {}

  def getListOfElements(self, protocol, scipionSet):
    eleList = []
    if scipionSet is not None:
      for item in scipionSet:
        elemType = 'Linker' if hasattr(item, '_type') and getattr(item, '_type') == 'Linker' else 'Epitope'
        elem = f'{elemType} (ID {item.getObjId()}): {item.getROISequence()}'
        eleList.append(elem)
    return eleList

SelectMultiEpitopeElementWizard().addTarget(protocol=ProtModifyMultiEpitope,
                                 targets=['inROI'],
                                 inputs=['inputMultiEpitope'],
                                 outputs=['inROI'])

class SelectSetMultiPointerWizard(SelectElementWizard):
  """Lists the items in a multipointer of SetOfX and choose one"""
  _targets, _inputs, _outputs = [], {}, {}

  def displayDialog(self, form, inputParam):
    protocol = form.protocol
    try:
      scipionSet = getattr(protocol, inputParam[0])
      listOfElements = self.getListOfElements(protocol, scipionSet)
    except Exception as e:
      print("ERROR: ", e)
      return

    finalList = []
    for i in listOfElements:
      finalList.append(String(i))
    provider = ListTreeProviderString(finalList)
    dlg = dialog.ListDialog(form.root, "MultiPointer sets", provider,
                            "Select one of the sets in the multipointer")
    return dlg

  def getListOfElements(self, protocol, scipionMultiPointer):
    eleList = []
    if scipionMultiPointer is not None:
      for i, element in enumerate(scipionMultiPointer):
        eleList.append(f'{i}//{element.get().__str__()}')
    return eleList

SelectSetMultiPointerWizard().addTarget(protocol=ProtocolRankDocking,
                                        targets=['defineInput'],
                                        inputs=['inputMoleculesSets'],
                                        outputs=['defineInput'])

SelectSetMultiPointerWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                        targets=['inSet'],
                                        inputs=['inputROISets'],
                                        outputs=['inSet'])
SelectSetMultiPointerWizard().addTarget(protocol=ProtOptimizeMultiEpitope,
                                        targets=['linkProtSet'],
                                        inputs=['inputROISets'],
                                        outputs=['linkProtSet'])


SelectSetMultiPointerWizard().addTarget(protocol=ProtocolRANXFuse,
                                        targets=['inSetID'],
                                        inputs=['inputSets'],
                                        outputs=['inSetID'])

class SelectElementMultiPointerWizard(SelectElementWizard):
    """Lists the items in a multipointer of SetOfX and choose one"""
    _targets, _inputs, _outputs = [], {}, {}

    def getListOfElements(self, protocol, scipionSet, i):
      eleList = []
      if scipionSet is not None:
        for element in scipionSet:
            eleList.append('Set {} | {}'.format(i, element.__str__()))
      return eleList

    def show(self, form, *params):
      protocol = form.protocol
      inputParam, outputParam = self.getInputOutput(form)
      try:
        listOfElements = []
        scipionMultiSet = getattr(protocol, inputParam[0])
        for i, scipionPointer in enumerate(scipionMultiSet):
            listOfElements += self.getListOfElements(protocol, scipionPointer.get(), i)
      except Exception as e:
        print("ERROR: ", e)
        return

      finalList = []
      for i in listOfElements:
        finalList.append(String(i))
      provider = ListTreeProviderString(finalList)
      dlg = dialog.ListDialog(form.root, "Sets items", provider,
                              "Select one of items in the sets")
      form.setVar(outputParam[0], dlg.values[0].get())


SelectElementWizard().addTarget(protocol=ProtDefineStructROIs,
                               targets=['ligName'],
                               inputs=['inSmallMols'],
                               outputs=['ligName'])

SelectElementWizard().addTarget(protocol=ProtocolShapeDistancesFiltering,
                               targets=['inputReferenceMolecule'],
                               inputs=['inputRefSmallMolecules'],
                               outputs=['inputReferenceMolecule'])

SelectElementWizard().addTarget(protocol=ProtocolFingerprintFiltering,
                               targets=['inputReferenceMolecule'],
                               inputs=['inputRefSmallMolecules'],
                               outputs=['inputReferenceMolecule'])

SelectElementWizard().addTarget(protocol=ProtocolRMSDDocking,
                                targets=['refMolName'],
                                inputs=['refSmallMolecules'],
                                outputs=['refMolName'])

SelectElementWizard().addTarget(protocol=ProtSeqCalculateConservation,
                               targets=['outSeq'],
                               inputs=['inputSequences'],
                               outputs=['outSeq'])

SelectElementMultiPointerWizard().addTarget(protocol=ProtocolPharmacophoreModification,
                                           targets=['currentFeatures'],
                                           inputs=['inputPharmacophores'],
                                           outputs=['currentFeatures'])


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
        pdbFile, AS, addPointer = '', False, True
        if issubclass(type(inputObj), str):
            outStr = [inputObj]
            AS, addPointer = True, False
        elif issubclass(type(inputObj), AtomStruct):
            pdbFile = inputObj.getFileName()
            outStr = [os.path.splitext(os.path.basename(pdbFile))[0]]
            AS = True
        elif issubclass(type(inputObj), Sequence):
            outStr = [inputObj.getId().replace("|", "_")]

        if AS:
            # Chain
            chainJson = getattr(protocol, inputParams[1]).get()
            chainId = json.loads(chainJson)['chain']
        else:
            chainJson = ''

        # Positions
        posJson = getattr(protocol, inputParams[2]).get()
        if posJson:
            posIdxs = json.loads(posJson)['index']
            seq = json.loads(posJson)['residues']
            outStr += [posIdxs]
        else:
            outStr += ['FIRST-LAST']
            finalResiduesList = self.getResidues(form, inputObj, chainJson)
            idxs = [json.loads(finalResiduesList[0].get())['index'], json.loads(finalResiduesList[-1].get())['index']]
            seq = self.getSequence(finalResiduesList, idxs)

        chainStr, chainFileId = '', ''
        if AS:
          chainStr = ', "chain": "{}"'.format(chainId)
          chainFileId = '_{}'.format(chainId)

        prevStr = getattr(protocol, outputParam[0]).get()
        lenPrev = len(prevStr.strip().split('\n')) + 1
        if prevStr.strip() == '':
          lenPrev -= 1
        elif not prevStr.endswith('\n'):
          prevStr += '\n'

        seqFile = protocol.getProject().getTmpPath('{}{}_{}.fa'.format(outStr[0], chainFileId, outStr[1]))
        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(outStr[0], seq))

        jsonStr = '%s) {"name": "%s"%s, "index": "%s", "seqFile": "%s"}\n' % \
                  (lenPrev, outStr[0], chainStr, outStr[1], seqFile)
        form.setVar(outputParam[0], prevStr + jsonStr)

        if addPointer:
            outPointers = outputParam[1]
            prevPointers = getattr(protocol, outPointers)
            prevPointers.append(Pointer(inputObj))
            form.setVar(outPointers, prevPointers)



AddSequenceWizard().addTarget(protocol=ProtDefineSetOfSequences,
                              targets=['addInput'],
                              inputs=[{'inputOrigin': ['inputSequence', 'inputAtomStruct', 'inputPDB']},
                                      'inpChain', 'inpPositions'],
                              outputs=['inputList', 'inputPointers'])