#######################################################################
# This script contains functions from PyBioMed described in the       #
# published paper: Dong et al. J Cheminform  (2018) 10:16             #
# https://doi.org/10.1186/s13321-018-0270-2. It is intended to        #
# replicate their calculations of Kier Flexibilities of small         #
# molecules to provide features for our scoring functions.            #
#                                                                     #
# Format conversion function at the end of the script added by        #
# @milesmcgibbon                                                      #
#                                                                     #
#######################################################################


from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from openbabel import openbabel as ob
import subprocess
# from oddt import toolkit,shape
# from molent.molent import entropy, binary_similarity, atomic_smiles

def CalculateKappaAlapha1(mol):
    """
    #################################################################
    Calculation of molecular shape index for one bonded fragment
    with Alapha
    ---->kappam1
    Usage:
        result=CalculateKappaAlapha1(mol)
        Input: mol is a molecule object.
        Output: result is a numeric value.
    #################################################################
    """
    P1 = mol.GetNumBonds(onlyHeavy=1)
    A = mol.GetNumHeavyAtoms()
    alpha = rdMolDescriptors.CalcHallKierAlpha(mol)
    denom = P1 + alpha
    if denom:
        kappa = (A + alpha) * (A + alpha - 1) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateKappaAlapha2(mol):
    """
    #################################################################
    Calculation of molecular shape index for two bonded fragment
    with Alapha
    ---->kappam2
    Usage:
        result=CalculateKappaAlapha2(mol)
        Input: mol is a molecule object.
        Output: result is a numeric value.
    #################################################################
    """
    P2 = len(Chem.FindAllPathsOfLengthN(mol, 2))
    A = mol.GetNumHeavyAtoms()
    alpha = rdMolDescriptors.CalcHallKierAlpha(mol)
    denom = P2 + alpha
    if denom:
        kappa = (A + alpha - 1) * (A + alpha - 2) ** 2 / denom ** 2
    else:
        kappa = 0.0
    return round(kappa, 3)

def CalculateFlexibility(mol):
    """
    #################################################################
    Calculation of Kier molecular flexibility index
    ---->phi
    Usage:
        result = CalculateFlexibility(mol)
        Input: mol is a molecule object.
        Output: result is a numeric value.`
    #################################################################
    """
    kappa1 = CalculateKappaAlapha1(mol)
    kappa2 = CalculateKappaAlapha2(mol)
    A = mol.GetNumHeavyAtoms()
    phi = kappa1 * kappa2 / (A + 0.0)
    return phi


def SmilePrep(file, molecule_encoder=None):
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("pdbqt", "pdb")
    mol = ob.OBMol()

    # Read and convert the molecule
    if not obConversion.ReadString(mol, file):
        raise ValueError("Failed to read the PDBQT file with Open Babel.")

    outMDL = obConversion.WriteString(mol)

    # Convert to RDKit molecule
    refmol = Chem.MolFromPDBBlock(outMDL, removeHs=False, sanitize=False)
    if refmol is None:
        raise ValueError("Failed to convert PDB block to RDKit Mol.")

    # Sanitize the molecule if necessary
    # try:
    #     Chem.SanitizeMol(refmol)
    #
    # except Exception as e:
    #     raise ValueError(f"Sanitization failed: {str(e)}")

    # Check and use the molecule encoder if provided
    if molecule_encoder:
        mol_feature = molecule_encoder(refmol)
        return refmol, mol_feature
    else:
        return refmol




