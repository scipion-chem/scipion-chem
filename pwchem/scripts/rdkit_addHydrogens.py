from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys, os

from pwchem.utils import getBaseName
from pwchem.utils.rdkitUtils import writeMol, getMolFilesDic
from pwchem.utils.scriptUtils import parseParams

def fixLigand(mol):
    """
    Add explicit hydrogens to carbon atoms and ensure tetravalent nitrogens are assigned a +1 charge
    """
    chemProblems = Chem.DetectChemistryProblems(mol)
    for chem_problem in chemProblems:
        if chem_problem.GetType() == 'AtomValenceException':
            atom = mol.GetAtomWithIdx(chem_problem.GetAtomIdx())
            if atom.GetSymbol() == 'N' and atom.GetFormalCharge() == 0 and atom.GetExplicitValence() == 4:
                atom.SetFormalCharge(1)
    Chem.SanitizeMol(mol)
    return Chem.AddHs(mol, addCoords=True)


#################################################################################################################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']

    molsDict, _ = getMolFilesDic(molFiles)
    for mol in molsDict:
        molBase = getBaseName(molsDict[mol])
        outBase = os.path.join(outDir, molBase)

        if paramsDic['doHydrogens']:
            mol = fixLigand(mol)
        if paramsDic['doGasteiger']:
            AllChem.ComputeGasteigerCharges(mol)

        setMolName = not mol.HasProp('_Name') or not mol.GetProp('_Name')
        if mol:
            outFile = outBase + '.sdf'
            writeMol(mol, outFile, setName=setMolName)
