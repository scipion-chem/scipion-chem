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
from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pyworkflow.wizard import Wizard
from pwem.wizards import GetStructureChainsWizard, EmWizard
from pwem.convert import AtomicStructHandler
import os, requests
from pwchem.protocols import ProtDefinePockets, ProtChemPdbFastaAlignment

class GetChainsWizard(Wizard):
    """
    This wizard will extract the chains from a atomic structure (pdb) file in
    order to select it in the protocol.
    Then, it will load the structure and will take all chain related
    information such as name and number of residues.
    """

    # list with tuples to target protocol parameters
    _targets = []


    def getChain_Information(self, protocol):
        """ Function that return the information about the different
            chains in a atomic structure using the class AtomicStructHandler.
        """

        # Help about AtomicStructHandler:
        # Given an atomic structure returns two dictionaries:
        #             (1) for all models and respective chains (chainID and length of residues)
        #             (2) for each chain list of residues

        structureHandler = AtomicStructHandler()

        try:
            filename = os.path.abspath(protocol.inputAtomStruct.get().getFileName())
            structureHandler.read(filename)
            structureHandler.getStructure()

            chains, residues = structureHandler.getModelsChains()

            chainInf = []  # List to save the chain information in a "list"
            n_chains = len(list(chains.items())[0][1])
            for model, chainDic in chains.items():
                for ID, resnumber in chainDic.items():
                    chainInf.append(
                        '{"Chain": "%s", "Number of residues": %d, "Number of chains": %d}' %
                        (str(ID), resnumber, n_chains))
            return chainInf
        except:
          print('Unable to read the input atom structure')

    def show(self, form, *params):
        """Show and select a chain"""

        protocol = form.protocol
        try:
            chainInf = self.getChain_Information(protocol)
        except Exception as e:
            print("ERROR: ", "A pdb file was not entered in the Atomic structure field. Please enter it.", e)
            return

        # This are the chains:
        chainList = []
        for chain in chainInf:
            chainList.append(String(chain))

        # Get a data provider from the greetings to be used in the tree (dialog)
        provider = ListTreeProviderString(chainList)
        dlg = dialog.ListDialog(form.root, "Model chains", provider,
                                "Select one of the chains")
        form.setVar('chain_name', dlg.values[0].get())


class SelectChainWizard(GetStructureChainsWizard):

      _targets = [(ProtDefinePockets, ['chain_name']),
                  (ProtChemPdbFastaAlignment, ['chain_name'])]

      @classmethod
      def getModelsChainsStep(cls, protocol, atomStructName='inputAtomStruct'):
        """ Returns (1) list with the information
           {"model": %d, "chain": "%s", "residues": %d} (modelsLength)
           (2) list with residues, position and chain (modelsFirstResidue)"""
        structureHandler = AtomicStructHandler()
        fileName = ""
        if hasattr(protocol, 'pdbId'):
          if protocol.pdbId.get() is not None:
            pdbID = protocol.pdbId.get()
            url = "https://www.rcsb.org/structure/"
            URL = url + ("%s" % pdbID)
            try:
              response = requests.get(URL)
            except:
              raise Exception("Cannot connect to PDB server")
            if (response.status_code >= 400) and (response.status_code < 500):
              raise Exception("%s is a wrong PDB ID" % pdbID)
            fileName = structureHandler.readFromPDBDatabase(
              os.path.basename(pdbID), dir="/tmp/")
          else:
            fileName = protocol.pdbFile.get()
        else:
          AS = getattr(protocol, atomStructName).get()
          if AS is not None:
            fileName = os.path.abspath(AS.getFileName())

        structureHandler.read(fileName)
        structureHandler.getStructure()
        # listOfChains, listOfResidues = structureHandler.getModelsChains()
        return structureHandler.getModelsChains()

      def show(self, form, *params):
        protocol = form.protocol
        try:
          listOfChains, listOfResidues = self.getModelsChainsStep(protocol)
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
        form.setVar('chain_name', dlg.values[0].get())


class SelectChainWizard2(SelectChainWizard):
    _targets = [(ProtChemPdbFastaAlignment, ['chain_name2'])]

    def show(self, form, *params):
        protocol = form.protocol
        try:
            listOfChains, listOfResidues = self.getModelsChainsStep(protocol, atomStructName='inputAtomStruct2')
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
        form.setVar('chain_name2', dlg.values[0].get())


class SelectResidueWizard(SelectChainWizard):
      _targets = [(ProtDefinePockets, ['resPosition'])]

      def editionListOfResidues(self, modelsFirstResidue, model, chain):
        self.residueList = []
        for modelID, chainDic in modelsFirstResidue.items():
          if int(model) == modelID:
            for chainID, seq_number in chainDic.items():
              if chain == chainID:
                for i in seq_number:
                  self.residueList.append(
                    '{"residue": %d, "%s"}' % (i[0], str(i[1])))

      def getResidues(self, form):
        protocol = form.protocol
        try:
          modelsLength, modelsFirstResidue = self.getModelsChainsStep(protocol)
        except Exception as e:
          print("ERROR: ", e)
          return
        selection = protocol.chain_name.get()

        model = selection.split(',')[0].split(':')[1].strip()
        chain = selection.split(',')[1].split(':')[1].split('"')[1]
        self.editionListOfResidues(modelsFirstResidue, model, chain)
        finalResiduesList = []
        for i in self.residueList:
          finalResiduesList.append(String(i))
        return finalResiduesList

      def show(self, form, *params):
        finalResiduesList = self.getResidues(form)
        provider = ListTreeProviderString(finalResiduesList)
        dlg = dialog.ListDialog(form.root, "Chain residues", provider,
                                "Select one residue (residue number, "
                                "residue name)")
        form.setVar('resPosition', dlg.values[0].get())


class AddResidueWizard(EmWizard):
    _targets = [(ProtDefinePockets, ['addResidue'])]

    def show(self, form, *params):
        protocol = form.protocol
        chain, pos = protocol.chain_name.get(), protocol.resPosition.get()
        form.setVar('inResidues', protocol.inResidues.get() +
                    '{} | {}\n'.format(chain, pos))

