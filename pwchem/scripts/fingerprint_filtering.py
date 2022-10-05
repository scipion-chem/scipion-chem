from __future__ import print_function
import sys
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import (
    PandasTools,
    Draw,
    Descriptors,
    MACCSkeys,
    rdFingerprintGenerator,
)
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

####Funcion para filtrar que se empleara en los dos filtros

def init(molFiles):
    mols_dict, mols = preprocessLigands(molFiles)
    dfObj = pd.DataFrame(mols_dict.items(), columns=['ROMol', 'ChEMBL'])
    return dfObj


def calculate_coeficient(df, maccs_fp_query, maccs_fp_list, morgan_fp_query, morgan_fp_list):
    df["tanimoto_maccs"] = DataStructs.BulkTanimotoSimilarity(maccs_fp_query, maccs_fp_list)
    df["tanimoto_morgan"] = DataStructs.BulkTanimotoSimilarity(morgan_fp_query, morgan_fp_list)

    df["dice_maccs"] = DataStructs.BulkDiceSimilarity(maccs_fp_query, maccs_fp_list)
    df["dice_morgan"] = DataStructs.BulkDiceSimilarity(morgan_fp_query, morgan_fp_list)
    final_df = df[df.columns[1:]]
    final_df.to_csv("all_distances.csv")
    return df


def filter(df, fingerprint, coefficient, cut_off):
    column = (str(coefficient).lower() + "_" + str(fingerprint).lower())

    df_final = df[['ChEMBL', str(column)]]
    dictionary = dict(df_final.values)

    filtered_fp = {}

    for chembl, coef in dictionary.items():
        if float(coef) >= float(cut_off):
            filtered_fp[chembl] = coef

    return filtered_fp


if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''

    paramsDic = parseParams(sys.argv[1])
    fingerprint_type = paramsDic['fingerprint']
    coefficient_type = paramsDic['coefficient']
    mol_files = paramsDic["ligandFiles"]
    cut = float(paramsDic['cut-off'])
    objective = paramsDic["referenceFile"]
    df = init(mol_files)
    query = preprocessObjective(objective)
    ########################################################################################
    maccs_fp_query = MACCSkeys.GenMACCSKeys(query)
    maccs_fp_list = df["ROMol"].apply(MACCSkeys.GenMACCSKeys).tolist()
    morgan_fp_query = rdFingerprintGenerator.GetCountFPs([query])[0]
    morgan_fp_list = rdFingerprintGenerator.GetCountFPs(df["ROMol"].tolist())
    ######################################################################################
    df2 = calculate_coeficient(df, maccs_fp_query, maccs_fp_list, morgan_fp_query, morgan_fp_list)

    filtered_fp = filter(df2, fingerprint_type, coefficient_type, cut)

    ############################### WE have a dict with the molecules that have passed the filter

    f = open(paramsDic["outputPath"], "w")
    f.write("#The following molecules have passed the fingerprint filter: " + "\n")
    for key, item in filtered_fp.items():
        f.write(key + "," + str(item) + "\n")

    f.close()



