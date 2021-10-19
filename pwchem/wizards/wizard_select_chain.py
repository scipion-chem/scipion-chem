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

from pwem.convert import AtomicStructHandler
import os

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