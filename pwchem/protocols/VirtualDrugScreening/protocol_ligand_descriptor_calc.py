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

# Script to calculate selected RDKit molecular descriptors for a set of input molecules
# Designed to be called from a Scipion protocol using a separate RDKit environment

import json
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from constants import descriptor_categories  # Si no funciona lo pego directamente

def load_molecule(path):
    return Chem.MolFromMolFile(path, sanitize=True)

def main():
    import sys
    input_json = sys.argv[1]
    output_json = sys.argv[2]
    enabled_categories = json.loads(sys.argv[3])  # Dict that indicates which descriptor categories are enabled
    
    with open(input_json, 'r') as f:
        molecule_paths = json.load(f)

    descriptor_list = Descriptors.descList
    header = [name for name, _ in descriptor_list]
    property_dict = {}

    for mol_name, path in molecule_paths.items():
        rdkit_mol = load_molecule(path)
        if rdkit_mol is None:
            continue

        row = []
        for name, desc_func in descriptor_list:
            evaluate = False
            for cat in descriptor_categories:
                if enabled_categories.get(cat, False) and name in descriptor_categories[cat]:
                    evaluate = True
                    break

            if evaluate:
                try:
                    value = desc_func(rdkit_mol)
                except Exception:
                    value = float('nan')
            else:
                value = float('nan')

            row.append(value)

        property_dict[mol_name] = row

    with open(output_json, 'w') as f:
        json.dump({'header': header, 'property_dict': property_dict}, f)

if __name__ == "__main__":
    main()
