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
from pwem.wizards import SelectChainWizard, SelectResidueWizard, VariableWizard
from pwem.convert import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence, Pointer

from pwchem.protocols import *

from pwchem.objects import SequenceVariants, SetOfStructROIs
from pwchem.viewers.viewers_sequences import SequenceAliView
from pwchem.utils import RESIDUES3TO1, RESIDUES1TO3, runOpenBabel, natural_sort
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

class AddCoordinatesWizard(AddROIWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    prevList, lenPrev = self.getPrevList(protocol, outputParam)

    coords = []
    for i in range(3):
        coords.append(getattr(protocol, inputParams[i]).get())

    roiDef = '%s) Coordinate: {"X": %s, "Y": %s, "Z": %s}\n' % \
             (lenPrev, coords[0], coords[1], coords[2])

    form.setVar(outputParam[0], prevList + roiDef)

AddCoordinatesWizard().addTarget(protocol=ProtDefineStructROIs,
                                 targets=['addCoordinate'],
                                 inputs=['coordX', 'coordY', 'coordZ'],
                                 outputs=['inROIs'])

class AddResidueWizard(AddROIWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol
        chainDic, resDic = json.loads(getattr(protocol, inputParams[0]).get()), \
                           json.loads(getattr(protocol, inputParams[1]).get())

        prevList, lenPrev = self.getPrevList(protocol, outputParam)

        roiDef = '%s) Residues: {"model": %s, "chain": "%s", "index": "%s", "residues": "%s"}\n' % \
                 (lenPrev, chainDic['model'], chainDic['chain'], resDic['index'], resDic['residues'])

        form.setVar(outputParam[0], prevList + roiDef)

class SetResidueWizard(AddResidueWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    chainDic, resDic = json.loads(getattr(protocol, inputParams[0]).get()), \
                       json.loads(getattr(protocol, inputParams[1]).get())

    resDef = '{"model": %s, "chain": "%s", "index": "%s", "residues": "%s"}' % \
             (chainDic['model'], chainDic['chain'], resDic['index'], resDic['residues'])
    form.setVar(outputParam[0], resDef)

AddResidueWizard().addTarget(protocol=ProtDefineStructROIs,
                             targets=['addResidue'],
                             inputs=['chain_name', 'resPosition'],
                             outputs=['inROIs'])

class AddLigandWizard(AddROIWizard):
  _targets, _inputs, _outputs = [], {}, {}
  
  def getPrevPointersIds(self, prevPointers):
      ids = []
      for p in prevPointers:
          ids.append(p.get().getObjId())
      return ids
  
  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    prevList, lenPrev = self.getPrevList(protocol, outputParam)
    
    if not getattr(protocol, inputParams[0]):
        roiDef = '%s) Ligand: {"molName": "%s"}\n' % \
                 (lenPrev, getattr(protocol, inputParams[3]).get())
    else:
        prevPointers = getattr(protocol, outputParam[1])
        prevIds = self.getPrevPointersIds(prevPointers)
        newObj = getattr(protocol, inputParams[1]).get()
        newId = newObj.getObjId()

        if not newId in prevIds:
            newIndex = len(prevPointers)
            prevPointers.append(Pointer(newObj))
        else:
            newIndex = prevIds.index(newId)

        roiDef = '%s) Ext-Ligand: {"pointerIdx": "%s", "ligName": "%s"}\n' % \
                 (lenPrev, newIndex, getattr(protocol, inputParams[2]).get())
        form.setVar(outputParam[1], prevPointers)

    form.setVar(outputParam[0], prevList + roiDef)



AddLigandWizard().addTarget(protocol=ProtDefineStructROIs,
                            targets=['addLigand'],
                            inputs=['extLig', 'inSmallMols', 'ligName', 'molName'],
                            outputs=['inROIs', 'inputPointers'])

class AddPPIsWizard(AddROIWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol

    prevList, lenPrev = self.getPrevList(protocol, outputParam)

    chain1Dic, chain2Dic = json.loads(getattr(protocol, inputParams[0]).get()), \
                           json.loads(getattr(protocol, inputParams[1]).get())
    iDist = getattr(protocol, inputParams[2]).get()

    c1, c2 = '{}-{}'.format(chain1Dic['model'], chain1Dic['chain']), \
             '{}-{}'.format(chain2Dic['model'], chain2Dic['chain'])

    roiDef = '%s) PPI: {"chain1": "%s", "chain2": "%s", "interDist": "%s"}\n' % \
             (lenPrev, c1, c2, iDist)

    form.setVar(outputParam[0], prevList + roiDef)

AddPPIsWizard().addTarget(protocol=ProtDefineStructROIs,
                          targets=['addPPI'],
                          inputs=['chain_name', 'chain_name2', 'ppiDistance'],
                          outputs=['inROIs'])

class AddNResidueWizard(AddROIWizard):
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParams, outputParam = self.getInputOutput(form)
    protocol = form.protocol
    resTypes, resDist = getattr(protocol, inputParams[0]).get(), \
                        getattr(protocol, inputParams[1]).get()
    resLink = protocol.getEnumText(inputParams[2])

    prevList, lenPrev = self.getPrevList(protocol, outputParam)

    roiDef = '%s) Near_Residues: {"residues": "%s", "distance": "%s", "linkage": "%s"}\n' % \
             (lenPrev, resTypes, resDist, resLink)

    form.setVar(outputParam[0], prevList + roiDef)


AddNResidueWizard().addTarget(protocol=ProtDefineStructROIs,
                              targets=['addNRes'],
                              inputs=['resNRes', 'resDistance', 'linkNRes'],
                              outputs=['inROIs'])


######################## Variable wizards ####################

class SelectChainWizardQT(SelectChainWizard):
  _targets, _inputs, _outputs = [], {}, {}
  @classmethod
  def getModelsChainsStep(cls, protocol, inputObj):
    """ Returns (1) list with the information
       {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
       (2) list with residues, position and chain (modelsFirstResidue)"""
    structureHandler = AtomicStructHandler()
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

    elif str(type(inputObj).__name__) == 'SetOfStructROIs':
        fileName = inputObj.getProteinFile()

    else:
        fileName = os.path.abspath(inputObj.getFileName())

    if fileName.endswith('.pdbqt'):
      inName, inExt = os.path.splitext(os.path.basename(fileName))
      pdbFile = os.path.abspath(os.path.join(protocol.getProject().getPath(inName + '.pdb')))
      args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(fileName), pdbFile)
      runOpenBabel(protocol=protocol, args=args, popen=True)
      fileName = pdbFile

    structureHandler.read(fileName)
    structureHandler.getStructure()
    return structureHandler.getModelsChains()

class SelectResidueWizardQT(SelectResidueWizard):
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
    if issubclass(type(inputObj), SetOfStructROIs):
        inputObj = inputObj.getProteinFile()
        inputObj = self.checkNoPDBQT(form, inputObj)
    elif issubclass(type(inputObj), AtomStruct):
        inputObj = inputObj.getFileName()
        inputObj = self.checkNoPDBQT(form, inputObj)

    return super().getResidues(form, inputObj, chainStr)


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

  def displayDialog(self, form):
    protocol = form.protocol
    inputParam, outputParam = self.getInputOutput(form)
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
    dlg = self.displayDialog(form)
    form.setVar(outputParam[0], dlg.values[0].get())

class SelectMultiElementWizard(SelectElementWizard):
  """Lists the items in a SetOfX and choose one or several"""
  _targets, _inputs, _outputs = [], {}, {}

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    dlg = self.displayDialog(form)
    values = [val.get() for val in dlg.values]
    form.setVar(outputParam[0], ','.join(values))


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