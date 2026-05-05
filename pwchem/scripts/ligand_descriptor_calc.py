# Script to calculate selected RDKit molecular descriptors for a set of input molecules
# Designed to be called from a Scipion protocol using a separate RDKit environment

# --- PATH bootstrap: make sure 'pwchem' is importable regardless of cwd/conda ---

# Script to calculate selected RDKit molecular descriptors for a set of input molecules
# Designed to be called from a Scipion protocol using a separate RDKit environment

import os, sys, json
from rdkit import Chem
from rdkit.Chem import Descriptors

HERE = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
for p in (REPO_ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)

try:
    from pwchem.constants import DESCRIPTOR_CATEGORIES
except ModuleNotFoundError:
    import importlib
    constants = importlib.import_module('constants')
    DESCRIPTOR_CATEGORIES = getattr(constants, 'descriptor_categories')

def load_molecule(path):
    path = os.path.abspath(path)
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"[WARNING] Skipping missing or empty file: {path}")
        return None
    try:
        if path.endswith(('.smi', '.smiles')):
            with open(path) as f:
                line = f.readline().strip()
            return Chem.MolFromSmiles(line)
        else:
            mol = Chem.MolFromMolFile(path, sanitize=True)
            if mol is None:
                print(f"[WARNING] RDKit could not read molecule: {path}")
            return mol
    except Exception as e:
        print(f"[ERROR] Failed to load {path}: {e}")
        return None

def main():
    if len(sys.argv) < 4:
        raise SystemExit("Usage: python ligand_descriptor_calc.py input.json output.json flags.json")

    input_json  = sys.argv[1]
    output_json = sys.argv[2]
    flags_json  = sys.argv[3]
    
    with open(input_json, 'r') as f:
        molecule_paths = json.load(f)

    with open(flags_json, 'r') as f:
        enabled_descriptors = set(json.load(f))

    descriptor_list = Descriptors.descList
    header = [name for name, _ in descriptor_list]
    property_dict = {}

    for mol_name, path in molecule_paths.items():
        mol = load_molecule(path)
        if mol is None:
            property_dict[mol_name] = [float('nan')] * len(header)
            continue

        row = []
        for name, desc_func in descriptor_list:
            if name not in enabled_descriptors:
                row.append(float('nan'))
                continue
            try:
                value = desc_func(mol)
            except Exception:
                value = float('nan')
            row.append(value)

        property_dict[mol_name] = row

    with open(output_json, 'w') as f:
        json.dump({'header': header, 'property_dict': property_dict}, f)

if __name__ == "__main__":
    main()