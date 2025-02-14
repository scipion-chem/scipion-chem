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

from utils import getMolFilesDic, parseParams, parseMoleculeFile

####Funcion para filtrar que se empleara en los dos filtros

def init(molFiles):
    molsDict, _ = getMolFilesDic(molFiles)
    dfObj = pd.DataFrame(molsDict.items(), columns=['ROMol', 'ChEMBL'])
    return dfObj


def calculateCoeficient(df, maccs_fp_query, maccs_fp_list, morgan_fp_query, morgan_fp_list):
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

    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    fingerprint_type = paramsDic['fingerprint']
    coefficient_type = paramsDic['coefficient']
    mol_files = paramsDic["ligandFiles"]
    cut = float(paramsDic['cut-off'])
    objective = paramsDic["referenceFile"]
    df = init(mol_files)
    query = parseMoleculeFile(objective)
    ########################################################################################
    maccs_fp_query = MACCSkeys.GenMACCSKeys(query)
    maccs_fp_list = df["ROMol"].apply(MACCSkeys.GenMACCSKeys).tolist()
    morgan_fp_query = rdFingerprintGenerator.GetCountFPs([query])[0]
    morgan_fp_list = rdFingerprintGenerator.GetCountFPs(df["ROMol"].tolist())
    ######################################################################################
    df2 = calculateCoeficient(df, maccs_fp_query, maccs_fp_list, morgan_fp_query, morgan_fp_list)

    filtered_fp = filter(df2, fingerprint_type, coefficient_type, cut)

    ############################### WE have a dict with the molecules that have passed the filter

    f = open(paramsDic["outputPath"], "w")
    f.write("#The following molecules have passed the fingerprint filter: " + "\n")
    for key, item in filtered_fp.items():
        f.write("{}\t{}\n".format(key, item))

    f.close()



