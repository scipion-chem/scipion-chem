from __future__ import print_function
import sys, os
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolAlign
import csv

###Funcion para procesar los archivos de input
def preprocessLigands(ligandsFiles):
    mols_dict = {}

    for molFile in ligandsFiles:
        if molFile.endswith('.mol2'):
            m = Chem.MolFromMol2File(molFile)
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

####Funcion para filtrar que se empleara en los dos filtros
def calculate_RMSD(paramsDict, reference_mol, filter_type):
    molFiles = paramsDic['ligandFiles']
    mols_dict, mols = preprocessLigands(molFiles)
    reference_dict, reference = preprocessLigands(reference_mol)

    result_dict = {}

    for mol in mols:

        try:
            if filter_type == "BestRMSD":
                # AllChem.GetBestRMS does a substructure mapping between the
                # reference and probe molecules so that it can try all symmetry-equivalent
                # mappings to get the best RMS value.

                matches = mol.GetSubstructMatches(mol, uniquify=False)
                maps = [list(enumerate(match)) for match in matches]
                result = AllChem.GetBestRMS(mol, reference, maps)

            else:
                result = rdkit.Chem.rdMolAlign.AlignMol(mol, reference)


        except:
            result = 'No sub-structure match found between the probe (' + str(mols_dict[mol]) + ') ' + 'and query mol'

        result_dict[mols_dict[mol]] = result

    return result_dict


def filter(num, result_dict):
    filtered_mol = []

    for mol, rmsd in result_dict.items():
        if str(rmsd).startswith("No sub-stucture"):
            pass
        else:
            if float(rmsd) <= num:
                filtered_mol.append(mol)

    return filtered_mol


def write_finalfile(filename, filtered_mol):
    with open(filename, 'w') as f:
        f.write("#Molecules that have passed the filter: \n")
        for molecule in filtered_mol:
            f.write(str(molecule) + "\n")


if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    filter_type = paramsDic['filterChoice']
    reference_mol = paramsDic['referenceFile']
    cut = float(paramsDic['cut-off'])

    result_dict = calculate_RMSD(paramsDic, reference_mol, filter_type)
    filtered_mol = filter(cut, result_dict)
    write_finalfile(paramsDic['outputPath'], filtered_mol)
