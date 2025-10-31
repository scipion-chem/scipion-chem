from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys, os

from utils import getMolFilesDic, parseParams, getBaseName, writeMol

DEFAULT_VALENCES = {
    'H': 1,
    'C': 4,
    'N': 3,
    'O': 2,
    'F': 1,
    'P': 3,
    'S': 2,
    'Cl': 1,
    'Br': 1,
    'I': 1,
    'B': 3,
}

def fixLigand(mol):
    """
    Add explicit hydrogens to carbon atoms and ensure tetravalent nitrogens are assigned a +1 charge
    """
    chemProblems = Chem.DetectChemistryProblems(mol)
    for chemProblem in chemProblems:
        if chemProblem.GetType() == 'AtomValenceException':
            atom = mol.GetAtomWithIdx(chemProblem.GetAtomIdx())
            symbol = atom.GetSymbol()
            valence = atom.GetTotalValence()  # total number of bonds (including to H)

            if symbol in DEFAULT_VALENCES:
                defValence = DEFAULT_VALENCES[symbol]
                if valence != defValence and atom.GetFormalCharge() == 0:
                    # Adjust formal charge to match valence difference
                    charge = defValence - valence
                    atom.SetFormalCharge(charge)
    Chem.SanitizeMol(mol)
    for a in mol.GetAtoms():
        rad = a.GetNumRadicalElectrons()
        if rad:
            a.SetNumExplicitHs(rad)
            a.SetNumRadicalElectrons(0)
    return Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)

#################################################################################################################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']
    outFormat = paramsDic['outputFormat'] if 'outputFormat' in paramsDic else '.sdf'
    outFormat = '.' + outFormat if not outFormat.startswith('.') else outFormat

    sanitize = eval(paramsDic['sanitize']) if 'sanitize' in paramsDic else True
    molsDict, _ = getMolFilesDic(molFiles, sanitize)
    for mol in molsDict:
        molBase = getBaseName(molsDict[mol])
        outBase = os.path.join(outDir, molBase)

        if paramsDic['doHydrogens']:
            mol = fixLigand(mol)
        if paramsDic['doGasteiger']:
            AllChem.ComputeGasteigerCharges(mol)

        setMolName = not mol.HasProp('_Name') or not mol.GetProp('_Name')
        if mol:
            outFile = outBase + outFormat
            writeMol(mol, outFile, setName=setMolName)
