from __future__ import print_function
import sys
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys, rdFingerprintGenerator

from utils import getMolFilesDic, parseParams, parseMoleculeFile

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

    if fingerprintType in ['MACCS', 'Both']:
        if coefficientType in ['Tanimoto', 'Both']:
            molsDF["tanimoto_maccs"] = DataStructs.BulkTanimotoSimilarity(maccsFpQuery, maccsFpMols)
        if coefficientType in ['Dice', 'Both']:
            molsDF["dice_maccs"] = DataStructs.BulkDiceSimilarity(maccsFpQuery, maccsFpMols)

    if fingerprintType in ['Morgan', 'Both']:
        if coefficientType in ['Tanimoto', 'Both']:
            molsDF["tanimoto_morgan"] = DataStructs.BulkTanimotoSimilarity(morganFpQuery, morganFpMols)
        if coefficientType in ['Dice', 'Both']:
            molsDF["dice_morgan"] = DataStructs.BulkDiceSimilarity(morganFpQuery, morganFpMols)
    
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
