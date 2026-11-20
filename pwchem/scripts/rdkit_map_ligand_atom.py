import argparse
import json

from rdkit import Chem


def loadMolecule(filename):
    """Load the first molecule from an SDF or MOL file."""

    if filename.lower().endswith('.sdf'):
        supplier = Chem.SDMolSupplier(filename, removeHs=False)
        mol = next((m for m in supplier if m is not None), None)
    else:
        mol = Chem.MolFromMolFile(filename, removeHs=False)

    if mol is None:
        raise ValueError('Could not read molecule: {}'.format(filename))

    if mol.GetNumConformers() == 0:
        raise ValueError(
            'Molecule has no coordinates: {}'.format(filename)
        )

    return mol


def heavyAtomIndices(mol):
    """Return heavy-atom indices in original atom order."""

    return [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 1
    ]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--reference', required=True)
    parser.add_argument('--pose', required=True)
    parser.add_argument('--atom-index', required=True, type=int)
    parser.add_argument('--output', required=True)
    args = parser.parse_args()

    reference = loadMolecule(args.reference)
    pose = loadMolecule(args.pose)

    reference_heavy = heavyAtomIndices(reference)

    if not 1 <= args.atom_index <= len(reference_heavy):
        raise ValueError(
            'Atom index {} is outside the reference range '
            '(1-{}).'.format(args.atom_index, len(reference_heavy))
        )

    # Convert the wizard's 1-based heavy-atom index to an RDKit index.
    reference_idx = reference_heavy[args.atom_index - 1]
    reference_atom = reference.GetAtomWithIdx(reference_idx)

    # Map reference atoms to the pose using molecular connectivity.
    match = pose.GetSubstructMatch(reference)

    if not match:
        raise ValueError(
            'Could not map reference molecule onto pose: {}'.format(
                args.pose
            )
        )

    poseIdx = match[reference_idx]
    poseAtom = pose.GetAtomWithIdx(poseIdx)

    if reference_atom.GetAtomicNum() != poseAtom.GetAtomicNum():
        raise ValueError(
            'Mapped atom element mismatch: reference={}, pose={}'.format(
                reference_atom.GetSymbol(),
                poseAtom.GetSymbol()
            )
        )

    position = pose.GetConformer().GetAtomPosition(poseIdx)

    result = {
        'coordinates': [position.x, position.y, position.z],
        'element': poseAtom.GetSymbol(),
        'pose_atom_index': poseIdx + 1,
        'reference_atom_index': args.atom_index
    }

    with open(args.output, 'w') as handle:
        json.dump(result, handle)


if __name__ == '__main__':
    main()