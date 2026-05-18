"""
Script to calculate selected RDKit molecular descriptors for a set of input molecules.
Designed to be called from a Scipion protocol using a separate RDKit environment.
"""

import os
import sys
import json
import importlib.util

from rdkit import Chem
from rdkit.Chem import Descriptors

HERE = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
for p in (REPO_ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

from rdkit_IO import getMolsFromFile

spec = importlib.util.spec_from_file_location(
    'constants',
    os.path.join(REPO_ROOT, 'pwchem', 'constants.py')
)
_constants = importlib.util.module_from_spec(spec)
spec.loader.exec_module(_constants)
DESCRIPTOR_CATEGORIES = _constants.DESCRIPTOR_CATEGORIES


def loadMolecule(path):
    """Load first molecule from file, handling smi files safely."""
    ext = os.path.splitext(path)[1].lower()
    if ext in ('.smi', '.smiles'):
        with open(path) as f:
            line = f.readline().strip().split()[0]
        return Chem.MolFromSmiles(line)
    else:
        mols, _ = getMolsFromFile(path)
        return mols[0] if mols else None


def getEnabledDescriptors(enabledCategories):
    """Return descList and set of enabled descriptor names."""
    descList = Descriptors.descList
    enabledSet = set()
    for name, _ in descList:
        for cat, props in DESCRIPTOR_CATEGORIES.items():
            if enabledCategories.get(cat, False) and name in props:
                enabledSet.add(name)
                break
    return descList, enabledSet


def calcMolDescriptors(mol, descList, enabledSet):
    """Calculate descriptors for a single molecule."""
    row = []
    for name, func in descList:
        if name in enabledSet:
            try:
                value = func(mol)
            except Exception:
                value = float('nan')
        else:
            value = float('nan')
        row.append(value)
    return row


def main():
    if len(sys.argv) < 4:
        raise SystemExit("Usage: python ligand_descriptor_calc.py input.json output.json flags.json")

    inputJson = sys.argv[1]
    outputJson = sys.argv[2]
    flagsJson = sys.argv[3]

    with open(flagsJson, 'r') as f:
        enabledCategories = json.load(f)

    with open(inputJson, 'r') as f:
        moleculePaths = json.load(f)

    descList, enabledSet = getEnabledDescriptors(enabledCategories)
    header = [name for name, _ in descList]
    propertyDict = {}

    for molName, path in moleculePaths.items():
        mol = loadMolecule(path)

        if mol is None:
            propertyDict[molName] = [float('nan')] * len(header)
            continue

        propertyDict[molName] = calcMolDescriptors(mol, descList, enabledSet)

    with open(outputJson, 'w') as f:
        json.dump({'header': header, 'property_dict': propertyDict}, f)


if __name__ == "__main__":
    main()