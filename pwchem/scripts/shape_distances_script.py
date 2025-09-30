from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys

from utils import getMolFilesDic, parseParams, parseMoleculeFile

#################################################################################################################
def distanceCalculation(molsDict, ref, distance, ignore, prealign, permuts):
    distanceDic = {}
    for mol, molFile in molsDict.items():
        if molFile.split('.')[-1] in ['smi', 'smile', 'smiles']:
            AllChem.EmbedMolecule(mol)
        if prealign:
            try:
                if permuts:
                    rmsd = AllChem.GetBestRMS(mol, ref)
                else:
                    rmsd = rdMolAlign.AlignMol(mol, ref)
            except:
                rmsd = 1000
                print('No substructure found for ', molFile)

        if distance == "Tanimoto":
            tanimoto = rdShapeHelpers.ShapeTanimotoDist(ref, mol, ignoreHs=ignore)
            distanceDic[molFile] = tanimoto
        elif distance == "Protrude":
            protude = rdShapeHelpers.ShapeProtrudeDist(ref, mol, ignoreHs=ignore)
            distanceDic[molFile] = protude

        elif distance == "RMSD":
            distanceDic[molFile] = rmsd

    return distanceDic

def writeFinalfile(filename, finalDict, distance):
    with open(filename, 'w') as f:
        f.write(f'MoleculeName\t{distance}\n')
        for molecule, coefficient in finalDict.items():
            f.write(str(molecule) + "\t" + str(coefficient) + "\n")


###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    molFiles = paramsDic['ligandFiles']
    molsDict, mols = getMolFilesDic(molFiles)
    objectiveFile = str(paramsDic["referenceFile"])
    objective = parseMoleculeFile(objectiveFile)
    distance = paramsDic["distanceType"]
    prealign, permuts = eval(paramsDic["prealign"]) if distance != "RMSD" else True, eval(paramsDic["prealignOrder"])
    ignore = eval(paramsDic["ignoreHydrogen"])

    distanceDic = distanceCalculation(molsDict, objective, distance, ignore, prealign, permuts)
    writeFinalfile(paramsDic['outputPath'], distanceDic, distance)

