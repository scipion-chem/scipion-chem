# Script to calculate selected RDKit molecular descriptors for a set of input molecules
# Designed to be called from a Scipion protocol using a separate RDKit environment

# --- PATH bootstrap: make sure 'pwchem' is importable regardless of cwd/conda ---

# Script to calculate selected RDKit molecular descriptors for a set of input molecules
# Designed to be called from a Scipion protocol using a separate RDKit environment

import os, sys, json
from rdkit import Chem
from rdkit.Chem import Descriptors

# --- PATH bootstrap ---
HERE = os.path.abspath(os.path.dirname(__file__))
REPO_ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
for p in (REPO_ROOT, HERE):
    if p not in sys.path:
        sys.path.insert(0, p)
# ----------------------

# Try to import descriptor_categories (fallback to constants)
try:
    from pwchem.scripts.constants import descriptor_categories
except ModuleNotFoundError:
    import importlib
    constants = importlib.import_module('constants')
    descriptor_categories = getattr(constants, 'descriptor_categories')

def load_molecule(path):
    path = os.path.abspath(path)

    if not os.path.exists(path) or os.path.getsize(path) == 0:
        print(f"[WARNING] missing file: {path}")
        return None

    try:
        if path.endswith(('.smi', '.smiles')):
            with open(path) as f:
                line = f.readline().strip()
            return Chem.MolFromSmiles(line)

        elif path.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(path, sanitize=True)
            mol = suppl[0] if len(suppl) > 0 else None
            if mol is None:
                print(f"[WARNING] Could not read SDF: {path}")
            return mol

        else:
            mol = Chem.MolFromMolFile(path, sanitize=True)
            if mol is None:
                print(f"[WARNING] MolFromMolFile failed: {path}")
            return mol

    except Exception as e:
        print(f"[ERROR] Failed to load {path}: {e}")
        return None



def main():
    if len(sys.argv) < 4:
        raise SystemExit("Usage: python ligand_descriptor_calc.py input.json output.json flags.json")

    input_json = sys.argv[1]
    output_json = sys.argv[2]
    flags_json = sys.argv[3]

    with open(flags_json, 'r') as f:
        enabled_categories = json.load(f)

    with open(input_json, 'r') as f:
        molecule_paths = json.load(f)

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
            evaluate = any(
                enabled_categories.get(cat, False) and name in props
                for cat, props in descriptor_categories.items()
            )

            try:
                value = desc_func(mol) if evaluate else float('nan')
            except Exception:
                value = float('nan')
            row.append(value)

        property_dict[mol_name] = row

    with open(output_json, 'w') as f:
        json.dump({'header': header, 'property_dict': property_dict}, f)


if __name__ == "__main__":
    main()
