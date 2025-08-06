# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import numpy as np
import os
import json

from pyworkflow.protocol import params
from pyworkflow.object import Float
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RDKIT_DIC
from pwchem.objects import SetOfSmallMolecules, SmallMoleculesLibrary
from pwchem.constants import DESCRIPTOR_CATEGORIES

class ProtocolLigandCharacterization(EMProtocol):
    _label = 'Ligand characterization'

    def _defineParams(self, form):
        form.addSection(label='Parameters')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules, SmallMoleculesLibrary', allowsNull=False,
                      label="Input Small Molecules",
                      help='Set of molecules to process.')
        # Name all parameters using booleans
        form.addParam('useConstitutional', params.BooleanParam, default=True,
                    label='Calculate constitutional descriptors')
        form.addParam('useElectronic', params.BooleanParam, default=True,
                    label='Calculate electronic descriptors')
        form.addParam('useTopological', params.BooleanParam, default=True,
                    label='Calculate topological descriptors')
        form.addParam('useGeometrical', params.BooleanParam, default=True,
                    label='Calculate geometrical descriptors')
        form.addParam('useRing', params.BooleanParam, default=True,
                    label='Calculate ring descriptors')
        form.addParam('useFragment', params.BooleanParam, default=True,
                    label='Calculate fragment descriptors')
        form.addParam('useOther', params.BooleanParam, default=True,
                    label='Calculate other descriptors')

    def _insertAllSteps(self):
        self._insertFunctionStep('runDescriptorCalc')
        self._insertFunctionStep('createOutputStep')

    def runDescriptorCalc(self):
        mol_set = self.inputSmallMolecules.get()

        # Build molecule dictionary: {molName: path}
        mol_dict = {mol.molName.get(): mol.getFileName() for mol in mol_set}
        input_json = self._getExtraPath("input_mols.json")
        with open(input_json, 'w') as f:
            json.dump(mol_dict, f)

        # Build descriptor flags
        flags = {cat: getattr(self, f"use{cat.capitalize()}") for cat in DESCRIPTOR_CATEGORIES}
        flags_json = json.dumps(flags)

        # Save flags as file (optional, not used directly by script)
        flags_path = self._getExtraPath("descriptor_flags.json")
        with open(flags_path, 'w') as f:
            f.write(flags_json)

        # Output file path
        output_json = self._getExtraPath("output_results.json")

        # Call the RDKit script in the correct environment
        Plugin.runScript(self, 'ligand_descriptor_calc.py',
                        f"{input_json} {output_json} '{flags_json}'",
                        env=RDKIT_DIC,
                        cwd=self._getPath())

        # Load the results
        with open(output_json, 'r') as f:
            data = json.load(f)

        self.header = data['header']
        self.propertyDict = data['property_dict']

    def createOutputStep(self):
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            if mol.molName.get() not in self.propertyDict:
                continue
            row = self.propertyDict[mol.molName.get()]
            
            for i in range(len(row)):
                prop_name = self.header[i]
                value = row[i]
                # Find the descriptor category
                category = 'uncategorized'
                for cat, prop_list in DESCRIPTOR_CATEGORIES.items():
                    if prop_name in prop_list:
                        category = cat
                        break
                attr_name = f"_Property_{category}_{prop_name}"
                setattr(mol, attr_name, Float(value))

            newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    def _summary(self):
        pass

    def getPluginCategory(self):
        return "Virtual Screening"
