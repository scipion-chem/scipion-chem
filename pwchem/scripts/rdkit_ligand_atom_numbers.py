#!/usr/bin/env python

import os
import sys

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import rdDepictor


def readMolecule(filename):
    """Read the first molecule while preserving its original atom order."""

    extension = os.path.splitext(filename)[1].lower()

    if extension == ".sdf":
        supplier = Chem.SDMolSupplier(
            filename,
            removeHs=False
        )

        for molecule in supplier:
            if molecule is not None:
                return molecule

        raise ValueError(
            "No valid molecule could be read from {}".format(filename)
        )

    if extension == ".mol":
        molecule = Chem.MolFromMolFile(
            filename,
            removeHs=False
        )

        if molecule is None:
            raise ValueError(
                "Could not read molecule from {}".format(filename)
            )

        return molecule

    raise ValueError(
        "Unsupported ligand format: {}".format(extension)
    )


def main():
    if len(sys.argv) != 3:
        print(
            "Usage: rdkit_ligand_atom_numbers.py "
            "<input_ligand.sdf|mol> <output.png>"
        )
        sys.exit(1)

    ligandFile = os.path.abspath(sys.argv[1])
    outputFile = os.path.abspath(sys.argv[2])

    if not os.path.isfile(ligandFile):
        raise FileNotFoundError(
            "Ligand file not found: {}".format(ligandFile)
        )

    molecule = readMolecule(ligandFile)

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

    outputDirectory = os.path.dirname(outputFile)

    if outputDirectory:
        os.makedirs(outputDirectory, exist_ok=True)

    Draw.MolToFile(
        molecule2d,
        outputFile,
        size=(1600, 1200),
        options=options
    )


if __name__ == "__main__":
    main()