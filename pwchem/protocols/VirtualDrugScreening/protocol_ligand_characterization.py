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
# * e-mail address 'scipion@cnb.csic.es'
# *

import os
import glob
import json

from pyworkflow.protocol import params
from pyworkflow.object import Float
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RDKIT_DIC
from pwchem.objects import SetOfSmallMolecules, SmallMoleculesLibrary
from pwchem.constants import DESCRIPTOR_CATEGORIES

OUTPUT_JSON = "output_results.json"


class ProtocolLigandCharacterization(EMProtocol):
    _label = 'Ligand characterization'

    def _defineParams(self, form):
        form.addSection(label='Parameters')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules, SmallMoleculesLibrary',
                      allowsNull=False,
                      label="Input Small Molecules",
                      help='Set of molecules or library to process.')

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

    def _buildMolDict(self):
        """Build {molName: filePath} dict handling both SetOfSmallMolecules and SmallMoleculesLibrary."""
        molInput = self.inputSmallMolecules.get()
        molDict = {}

        if isinstance(molInput, SmallMoleculesLibrary):
            inDir = os.path.abspath(self._getTmpPath())
            try:
                ligFiles = molInput.splitInFiles(inDir)
            except Exception as e:
                self.info(f"splitInFiles() failed ({type(e).__name__}: {e}). Trying directory glob.")
                libSrc = molInput.getFileName()
                ligFiles = []
                if os.path.isdir(libSrc):
                    for pat in ('*.sdf', '*.mol', '*.mol2', '*.smi', '*.smiles'):
                        ligFiles.extend(glob.glob(os.path.join(libSrc, pat)))
                    ligFiles = [os.path.abspath(p) for p in ligFiles]
                    if not ligFiles:
                        self.error(f"No ligand files found in directory: {libSrc}")
                else:
                    self.error(f"Library source is not a directory and splitInFiles() failed: {libSrc}")

            for f in ligFiles:
                base = os.path.basename(f)
                name, _ = os.path.splitext(base)
                i, orig = 2, name
                while name in molDict:
                    name = f"{orig}_{i}"
                    i += 1
                molDict[name] = os.path.abspath(f)
        else:
            molDict = {mol.molName.get(): os.path.join(self.getProject().getPath(), mol.getFileName())
                       for mol in molInput}

        return molDict

    def runDescriptorCalc(self):
        molDict = self._buildMolDict()

        inputJson = self._getExtraPath("input_mols.json")
        with open(inputJson, 'w') as f:
            json.dump(molDict, f)

        flags = {cat: bool(getattr(self, f"use{cat.capitalize()}").get()) for cat in DESCRIPTOR_CATEGORIES}
        flagsPath = self._getExtraPath("descriptor_flags.json")
        with open(flagsPath, 'w') as f:
            json.dump(flags, f)

        outputJson = self._getExtraPath(OUTPUT_JSON)
        scriptsDir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..', 'scripts'))
        if not os.path.isdir(scriptsDir):
            self.error(f"Scripts dir not found: {scriptsDir}")

        env = dict(RDKIT_DIC) if isinstance(RDKIT_DIC, dict) else {}
        prevPy = env.get('PYTHONPATH', '')
        env['PYTHONPATH'] = os.pathsep.join([p for p in [scriptsDir, prevPy] if p])

        Plugin.runScript(
            self,
            'ligand_descriptor_calc.py',
            f"{os.path.abspath(inputJson)} {os.path.abspath(outputJson)} {os.path.abspath(flagsPath)}",
            env=env,
            cwd=scriptsDir
        )

        if not os.path.exists(outputJson) or os.path.getsize(outputJson) == 0:
            self.error("Descriptor script produced no output.")

    def createOutputStep(self):
        molInput = self.inputSmallMolecules.get()
        outputJson = self._getExtraPath(OUTPUT_JSON)

        with open(outputJson, 'r') as f:
            data = json.load(f)

        header = data.get('header', [])
        propertyDict = data.get('property_dict', {})

        if isinstance(molInput, SmallMoleculesLibrary):
            newMols = SetOfSmallMolecules(filename=self._getPath('smallMolecules.sqlite'))
            for molName, molPath in self._buildMolDict().items():
                molCopy = SetOfSmallMolecules.ITEM_TYPE()
                molCopy.setFileName(molPath)
                molCopy.molName.set(molName)

                if molName in propertyDict:
                    row = propertyDict[molName]
                    for i, propName in enumerate(header):
                        value = row[i]
                        category = next(
                            (cat for cat, props in DESCRIPTOR_CATEGORIES.items() if propName in props),
                            'uncategorized'
                        )
                        attrName = f"_Property_{category}_{propName}"
                        setattr(molCopy, attrName, Float(value))

                newMols.append(molCopy)
        else:
            newMols = SetOfSmallMolecules.createCopy(molInput, self._getPath(), copyInfo=True)
            for mol in molInput:
                molCopy = mol.clone()
                name = molCopy.molName.get()

                if name not in propertyDict:
                    newMols.append(molCopy)
                    continue

                row = propertyDict[name]
                for i, propName in enumerate(header):
                    value = row[i]
                    category = next(
                        (cat for cat, props in DESCRIPTOR_CATEGORIES.items() if propName in props),
                        'uncategorized'
                    )
                    attrName = f"_Property_{category}_{propName}"
                    setattr(molCopy, attrName, Float(value))

                newMols.append(molCopy)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    def _summary(self):
        msgs = []
        outJson = self._getExtraPath(OUTPUT_JSON)
        if os.path.exists(outJson):
            msgs.append(f"Descriptors stored in {os.path.basename(outJson)}")
        else:
            msgs.append("No descriptor output found.")
        return msgs

    def getPluginCategory(self):
        return "Virtual Screening"