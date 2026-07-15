from __future__ import print_function
import sys
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from utils import getMolFilesDic, parseParams, parseMoleculeFile

TVERSKY_REF_A, TVERSKY_REF_B = 1.0, 0.0
TVERSKY_DB_A, TVERSKY_DB_B = 0.0, 1.0

####Funcion para filtrar que se empleara en los dos filtros

def getMolsDF(molFiles):
    molsDict, _ = getMolFilesDic(molFiles)
    dfObj = pd.DataFrame(molsDict.items(), columns=['ROMol', 'ChEMBL'])
    return dfObj

def performAnalysis(query, molsDF, fingerprintType, coefficientType, outPath):
    maccsFpQuery = MACCSkeys.GenMACCSKeys(query)
    maccsFpMols = molsDF["ROMol"].apply(MACCSkeys.GenMACCSKeys).tolist()
    morganFpQuery = rdFingerprintGenerator.GetCountFPs([query])[0]
    morganFpMols = rdFingerprintGenerator.GetCountFPs(molsDF["ROMol"].tolist())

    # 'All' runs every coefficient. 'Both' is kept for backwards compatibility
    # with older saved workflows (it means Tanimoto + Dice).
    doTanimoto = coefficientType in ['Tanimoto', 'Both', 'All']
    doDice     = coefficientType in ['Dice', 'Both', 'All']
    doTversky  = coefficientType in ['Tversky_Ref', 'All']
    doTverskyDB = coefficientType in ['Tversky_DB', 'All']

    if fingerprintType in ['MACCS', 'Both']:
        if doTanimoto:
            molsDF["tanimoto_maccs"] = DataStructs.BulkTanimotoSimilarity(maccsFpQuery, maccsFpMols)
        if doDice:
            molsDF["dice_maccs"] = DataStructs.BulkDiceSimilarity(maccsFpQuery, maccsFpMols)
        if doTversky:
            # IMPORTANT: reference (query) is the first arg -> this is Tversky_REF.
            molsDF["tversky_ref_maccs"] = DataStructs.BulkTverskySimilarity(
                maccsFpQuery, maccsFpMols, TVERSKY_REF_A, TVERSKY_REF_B)
        if doTverskyDB:
            molsDF["tversky_db_maccs"] = DataStructs.BulkTverskySimilarity(
                maccsFpQuery, maccsFpMols, TVERSKY_DB_A, TVERSKY_DB_B)

    if fingerprintType in ['Morgan', 'Both']:
        if doTanimoto:
            molsDF["tanimoto_morgan"] = DataStructs.BulkTanimotoSimilarity(morganFpQuery, morganFpMols)
        if doDice:
            molsDF["dice_morgan"] = DataStructs.BulkDiceSimilarity(morganFpQuery, morganFpMols)
        if doTversky:
            molsDF["tversky_ref_morgan"] = DataStructs.BulkTverskySimilarity(
                morganFpQuery, morganFpMols, TVERSKY_REF_A, TVERSKY_REF_B)
        if doTverskyDB:
            molsDF["tversky_db_morgan"] = DataStructs.BulkTverskySimilarity(
                morganFpQuery, morganFpMols, TVERSKY_DB_A, TVERSKY_DB_B)

    saveResults(molsDF, outPath)

def saveResults(df, outPath):
    finalDf = df[df.columns[1:]]
    finalDf.to_csv(outPath, sep='\t', index=False)

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''

    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    fingerprintType, coefficientType = paramsDic['fingerprint'], paramsDic['coefficient']
    molFiles, objFile = paramsDic["ligandFiles"], paramsDic["referenceFile"]

    molsDF = getMolsDF(molFiles)
    query = parseMoleculeFile(objFile)
    ########################################################################################
    performAnalysis(query, molsDF, fingerprintType, coefficientType, paramsDic["outputPath"])