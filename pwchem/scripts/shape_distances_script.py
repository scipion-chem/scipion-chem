from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys

def preprocessLigands(ligandsFiles):
    mols_dict = {}

    for molFile in ligandsFiles:
        if molFile.endswith('.mol2'):
            with open(molFile) as fIn:
                m = Chem.MolFromMol2Block(fIn.read(), sanitize=False)
            mols_dict[m] = molFile

        elif molFile.endswith('.mol'):
            m = Chem.MolFromMolFile(molFile)
            mols_dict[m] = molFile

        elif molFile.endswith('.pdb'):
            m = Chem.MolFromPDBFile(molFile)
            mols_dict[m] = molFile

        elif molFile.endswith('.smi'):
            f = open(molFile, "r")
            firstline = next(f)
            m = Chem.MolFromSmiles(str(firstline))
            mols_dict[m] = molFile

        elif molFile.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(molFile)
            for mol in suppl:
                mols_dict[mol] = molFile

    mols = list(mols_dict.keys())

    return mols_dict, mols

def preprocessObjective(objective_file):
    mol = " "
    if objective_file.endswith('.mol2'):
        mol = Chem.MolFromMol2File(objective_file)

    elif objective_file.endswith('.mol'):
        mol = Chem.MolFromMolFile(objective_file)

    elif objective_file.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(objective_file)

    elif objective_file.endswith('.smi'):
        f = open(objective_file, "r")
        firstline = next(f)
        mol = Chem.MolFromSmiles(str(firstline))

    elif objective_file.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(objective_file)
        for mol in suppl:
            mol = mol

    else:
        mol = Chem.MolFromSmiles(objective_file)


    return mol


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


#################################################################################################################
def distance_filter(mols_dict, objective, distance, ignore, prealign, permuts):
    distance_dict = {}
    ref = Chem.AddHs(objective)
    AllChem.EmbedMolecule(ref)
    for mol, molFile in mols_dict.items():
        mol = Chem.AddHs(mol)
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

        if distance == "Tanimoto Distance":
            tanimoto = rdShapeHelpers.ShapeTanimotoDist(ref, mol, ignoreHs=ignore)
            distance_dict[molFile] = tanimoto
        elif distance == "Protrude Distance":
            protude = rdShapeHelpers.ShapeProtrudeDist(ref, mol, ignoreHs=ignore)
            distance_dict[molFile] = protude

        elif distance == "RMSD":
            distance_dict[molFile] = rmsd

    return(distance_dict)

def filter(distance_dict, cut):
    final_dict = {}
    for chembl, number in distance_dict.items():
        if float(number) <= cut:
            final_dict[chembl] = number
    return(final_dict)

def write_finalfile(filename, final_dict):
    with open(filename, 'w') as f:
        f.write("# Molecules that have passed the shape filtering: \n")
        for molecule, coefficient in final_dict.items():
            f.write(str(molecule) + "\t" + str(coefficient) + "\n")


###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    molFiles = paramsDic['ligandFiles']
    mols_dict, mols = preprocessLigands(molFiles)
    objectiveFile = str(paramsDic["referenceFile"])
    objective = preprocessObjective(objectiveFile)
    distance = paramsDic["distanceType"]
    prealign, permuts = bool(paramsDic["prealign"]), bool(paramsDic["prealignOrder"])
    ignore = bool(paramsDic["ignoreHydrogen"])
    cut = float(paramsDic["cut-off"])

    distance_dict = distance_filter(mols_dict, objective, distance, ignore, prealign, permuts)
    final_dict = filter(distance_dict, cut)

    write_finalfile(paramsDic['outputPath'], final_dict)

