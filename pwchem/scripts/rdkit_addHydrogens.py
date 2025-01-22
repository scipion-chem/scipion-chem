from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys, os

def parseMoleculeFile(molFile):
    if molFile.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molFile)
    elif molFile.endswith('.mol'):
        mol = Chem.MolFromMolFile(molFile)
    elif molFile.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molFile)
    elif molFile.endswith('.smi'):
        f = open(molFile, "r")
        firstline = next(f)
        mol = Chem.MolFromSmiles(str(firstline))
    elif molFile.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(molFile)
        for mol in suppl:
            break
    else:
        mol = Chem.MolFromSmiles(molFile)

    return mol

def getMolFilesDic(molFiles):
    mols_dict = {}
    for molFile in molFiles:
        m = parseMoleculeFile(molFile)
        mols_dict[m] = molFile

    mols = list(mols_dict.keys())
    return mols_dict, mols


def parseParams(paramsFile):
    paramsDic = {}
    with open(paramsFile) as f:
        for line in f:
            key, value = line.strip().split(':')
            if key == 'ligandFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic

def fixLigand(mol):
    """
    Add explicit hydrogens to carbon atoms and ensure tetravalent nitrogens are assigned a +1 charge
    """
    chem_problems = Chem.DetectChemistryProblems(mol)
    for chem_problem in chem_problems:
        if chem_problem.GetType() == 'AtomValenceException':
            atom = mol.GetAtomWithIdx(chem_problem.GetAtomIdx())
            if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0 and atom.GetExplicitValence() == 4:
                atom.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    return Chem.AddHs(mol, addCoords=True)


#################################################################################################################

def writeMol(mol, outFile, cid=-1, setName=False):
    w = Chem.SDWriter(outFile)
    molName = os.path.split(os.path.splitext(outFile)[0])[-1]
    if setName:
        mol.SetProp('_Name', molName)
    w.write(mol, cid)
    w.close()

###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']

    mols_dict, _ = getMolFilesDic(molFiles)
    for mol in mols_dict:
        molBase = os.path.splitext(os.path.basename(mols_dict[mol]))[0]
        outBase = os.path.join(outDir, molBase)

        if paramsDic['doHydrogens']:
            mol = fixLigand(mol)
        if paramsDic['doGasteiger']:
            AllChem.ComputeGasteigerCharges(mol)

        setMolName = not mol.HasProp('_Name') or not mol.GetProp('_Name')
        if mol:
            outFile = outBase + '.sdf'
            writeMol(mol, outFile, setName=setMolName)
