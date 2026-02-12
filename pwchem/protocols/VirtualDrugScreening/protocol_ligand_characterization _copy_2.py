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
import glob

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
        form.addSection(label='Parameters COSS')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules, SmallMoleculesLibrary', allowsNull=False,
                      label="Input Small Molecules",
                      help='Set of molecules to process.')
        form.addParam('useLibrary', params.BooleanParam, default=False,
                      label="Use Small Molecules Library")

        # Boolean flags for descriptor categories
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

    # --------------------------------------------------------
    # Workflow
    # --------------------------------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runDescriptorCalc')
        self._insertFunctionStep('createOutputStep')

    # --------------------------------------------------------
    # Step 1: Run external RDKit descriptor calculator
    # --------------------------------------------------------
    def runDescriptorCalc(self):
        mol_input = self.inputSmallMolecules.get()
        is_lib = isinstance(mol_input, SmallMoleculesLibrary) or self.useLibrary.get()

        # Build {molName: filePath}
        mol_dict = {}
        if is_lib:
            inDir = os.path.abspath(self._getTmpPath())
            try:
                ligFiles = mol_input.splitInFiles(inDir)
            except Exception as e:
                self.info(f"splitInFiles() failed ({type(e).__name__}: {e}). Trying directory glob.")
                try:
                    lib_src = mol_input.getFileName()
                except Exception:
                    raise
                ligFiles = []
                if os.path.isdir(lib_src):
                    exts = ('*.sdf', '*.mol', '*.mol2', '*.smi', '*.smiles')
                    for pat in exts:
                        ligFiles.extend(glob.glob(os.path.join(lib_src, pat)))
                    ligFiles = [os.path.abspath(p) for p in ligFiles]
                    if not ligFiles:
                        self.error(f"No ligand files found in directory: {lib_src}")
                else:
                    self.error(f"Library source is not a directory and splitInFiles() failed: {lib_src}")

            for f in ligFiles:
                base = os.path.basename(f)
                name, _ = os.path.splitext(base)
                i, orig = 2, name
                while name in mol_dict:
                    name = f"{orig}_{i}"
                    i += 1
                mol_dict[name] = os.path.abspath(f)
        else:
            mol_dict = {mol.molName.get(): mol.getFileName() for mol in mol_input}

        # Serialize input molecules
        input_json = self._getExtraPath("input_mols.json")
        with open(input_json, 'w') as f:
            json.dump(mol_dict, f)

        # Descriptor flags
        flags = {cat: bool(getattr(self, f"use{cat.capitalize()}").get()) for cat in DESCRIPTOR_CATEGORIES}
        flags_path = self._getExtraPath("descriptor_flags.json")
        with open(flags_path, 'w') as f:
            json.dump(flags, f)

        # Output path
        output_json = self._getExtraPath("output_results.json")

        # Prepare environment
        scripts_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
        if not os.path.isdir(scripts_dir):
            self.error(f"Scripts dir not found: {scripts_dir}")

        env = dict(RDKIT_DIC) if isinstance(RDKIT_DIC, dict) else {}
        prev_py = env.get('PYTHONPATH', '')
        env['PYTHONPATH'] = os.pathsep.join([p for p in [scripts_dir, prev_py] if p])

        # Run descriptor calculation script
        Plugin.runScript(
            self,
            'ligand_descriptor_calc.py',
            f"{os.path.abspath(input_json)} {os.path.abspath(output_json)} {os.path.abspath(flags_path)}",
            env=env,
            cwd=scripts_dir
        )

        if not os.path.exists(output_json) or os.path.getsize(output_json) == 0:
            self.error("Descriptor script produced no output.")

    # --------------------------------------------------------
    # Step 2: Create output with descriptors visible in Scipion
    # --------------------------------------------------------
    def createOutputStep(self):
        inp = self.inputSmallMolecules.get()

        #if isinstance(inp, SmallMoleculesLibrary) or self.useLibrary.get():
        #    self.info("Input is a SmallMoleculesLibrary; skipping property attachment. "
        #              "Descriptor results are available in 'output_results.json'.")
        #   return

        newMols = SetOfSmallMolecules.createCopy(inp, self._getPath(), copyInfo=True)

        # Load descriptor data
        output_json = self._getExtraPath("output_results.json")
        with open(output_json, 'r') as f:
            data = json.load(f)

        header = data.get('header', [])
        property_dict = data.get('property_dict', {})

        for mol in inp:
            mol_copy = mol.clone()
            name = mol_copy.molName.get()

            if name not in property_dict:
                newMols.append(mol_copy)
                continue

            row = property_dict[name]
            for i, prop_name in enumerate(header):
                value = row[i]
                category = next(
                    (cat for cat, props in DESCRIPTOR_CATEGORIES.items() if prop_name in props),
                    'uncategorized'
                )
                attr_name = f"_Property_{category}_{prop_name}"
                setattr(mol_copy, attr_name, Float(value))

            newMols.append(mol_copy)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------------------------------------
    # Summary and plugin metadata
    # --------------------------------------------------------
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
