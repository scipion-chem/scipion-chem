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

import os
import json

from pyworkflow.protocol import params
from pyworkflow.object import Float
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RDKIT_DIC
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import DESCRIPTOR_CATEGORIES


class ProtocolLigandCharacterization(EMProtocol):
    _label = 'Ligand characterization'

    def _defineParams(self, form):
        form.addSection(label='Parameters')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules",
                      help='Set of molecules to process.')

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
        mol_input = self.inputSmallMolecules.get()
        mol_dict = {mol.molName.get(): os.path.abspath(mol.getFileName()) for mol in mol_input}

        input_json = self._getExtraPath("input_mols.json")
        with open(input_json, 'w') as f:
            json.dump(mol_dict, f)

        flags = {cat: bool(getattr(self, f"use{cat.capitalize()}").get()) for cat in DESCRIPTOR_CATEGORIES}
        flags_path = self._getExtraPath("descriptor_flags.json")
        with open(flags_path, 'w') as f:
            json.dump(flags, f)

        output_json = self._getExtraPath("output_results.json")

        scripts_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
        if not os.path.isdir(scripts_dir):
            self.error(f"Scripts dir not found: {scripts_dir}")

        env = dict(RDKIT_DIC) if isinstance(RDKIT_DIC, dict) else {}
        prev_py = env.get('PYTHONPATH', '')
        env['PYTHONPATH'] = os.pathsep.join([p for p in [scripts_dir, prev_py] if p])

        Plugin.runScript(
            self,
            'ligand_descriptor_calc.py',
            f"{os.path.abspath(input_json)} {os.path.abspath(output_json)} {os.path.abspath(flags_path)}",
            env=env,
            cwd=scripts_dir
        )

        if not os.path.exists(output_json) or os.path.getsize(output_json) == 0:
            self.error("Descriptor script produced no output.")

    def createOutputStep(self):
        mol_input = self.inputSmallMolecules.get()

        output_json = self._getExtraPath("output_results.json")
        with open(output_json, 'r') as f:
            data = json.load(f)

        header = data.get('header', [])
        property_dict = data.get('property_dict', {})

        new_mols = SetOfSmallMolecules.createCopy(mol_input, self._getPath(), copyInfo=True)

        for mol in mol_input:
            mol_copy = mol.clone()
            name = mol_copy.molName.get()

            if name not in property_dict:
                new_mols.append(mol_copy)
                continue

            row = property_dict[name]
            for i, prop_name in enumerate(header):
                value = row[i]
                category = next(
                    (cat for cat, props in DESCRIPTOR_CATEGORIES.items() if prop_name in props),
                    'uncategorized'
                )
                attr_name = f"_Property_{category}_{prop_name}"
                if not (isinstance(value, float) and value != value):  # skip nan
                    setattr(mol_copy, attr_name, Float(value))

            new_mols.append(mol_copy)

        new_mols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=new_mols)

    def _summary(self):
        msgs = []
        out_json = self._getExtraPath("output_results.json")
        if os.path.exists(out_json):
            msgs.append(f"Descriptors stored in {os.path.basename(out_json)}")
        else:
            msgs.append("No descriptor output found.")
        return msgs

    def getPluginCategory(self):
        return "Virtual Screening"