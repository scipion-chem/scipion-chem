#!/usr/bin/env python


import os
import sys

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor

from utils import parseMoleculeFile

# Fixed output name; the wizard looks for this file in the cwd it sets.
OUTPUT_IMAGE_NAME = "ligand_atom_numbers.png"



def main():
    if len(sys.argv) != 2:
        raise SystemExit(
            "Usage: rdkit_ligand_atom_numbers.py <ligand file (sdf, mol, mol2, pdb)>"
        )

    ligandFile = os.path.realpath(sys.argv[1])

    if not os.path.isfile(ligandFile):
        raise FileNotFoundError(
            "Ligand file not found: {}".format(ligandFile)
        )

    molecule = parseMoleculeFile(ligandFile)
    if molecule is None:
        raise ValueError("Could not read molecule from {}".format(ligandFile))

    # Assign numbers in original file atom order, skipping H and D.
    # These are 1-based heavy-atom numbers, not RDKit or PyMOL indices.
    atomLabels = {}
    heavyAtomNumber = 0

    for atom in molecule.GetAtoms():
        element = atom.GetSymbol().upper()

        # Exclude hydrogen and deuterium, including isotope-labelled H.
        if element in ("H", "D"):
            continue

        heavyAtomNumber += 1

        atomLabels[atom.GetIdx()] = "{}:{}".format(
            atom.GetSymbol(),
            heavyAtomNumber
        )

    if heavyAtomNumber == 0:
        raise ValueError("The ligand contains no heavy atoms.")

    # Generate a 2D depiction without changing the original atom ordering.
    molecule2d = Chem.Mol(molecule)
    rdDepictor.Compute2DCoords(molecule2d)

    options = Draw.MolDrawOptions()
    options.addAtomIndices = False
    options.bondLineWidth = 2.0
    options.fixedFontSize = 18
    options.padding = 0.12

    for atomIndex, label in atomLabels.items():
        options.atomLabels[atomIndex] = label

    # Saved in the working directory chosen by the caller through cwd.
    outputFile = os.path.join(os.getcwd(), OUTPUT_IMAGE_NAME)

    Draw.MolToFile(
        molecule2d,
        outputFile,
        size=(1600, 1200),
        options=options
    )

if __name__ == "__main__":
    main()