from __future__ import print_function
from rdkit.Chem import Descriptors, Draw, PandasTools
import sys

from utils import getMolFilesDic, parseParams

#################################################################################################################

####Funcion para filtrar que se empleara en los dos filtros
def filterExpression(expression, mols):
    """Filter molecule by a generic expression, such as `mol.logp > 1`."""

    out = []
    for mol in mols:
        if isinstance(expression, list):
            fail = 0
            for e in expression:
                if eval(e):
                    fail += 1

            if fail >= 3:
                out.append(mol)

        else:
            if eval(expression):
                out.append(mol)
    return out


def ro5Filter(mols):
    filtered_mols = []
    out = filterExpression(mols=mols, expression=['Descriptors.ExactMolWt(mol) < 500',
                                                  'Descriptors.NumHAcceptors(mol) <= 10',
                                                  'Descriptors.NumHDonors(mol) <= 5',
                                                  'Descriptors.MolLogP(mol) <= 5'])

    filtered_mols.append(out)
    return filtered_mols


def ro3Filter(mols):
    filtered_mols = []
    out = filterExpression(mols=mols, expression=['Descriptors.ExactMolWt(mol) < 300',
                                                  'Descriptors.NumHAcceptors(mol) <= 3',
                                                  'Descriptors.NumHDonors(mol) <= 3',
                                                  'Descriptors.MolLogP(mol) <= 3'])

    filtered_mols.append(out)
    return filtered_mols

###############################################################################
def writeDetails(moleculas, dict):
    with open("details.txt", 'w') as a:
        for molecula in moleculas:
            molecular_weight = Descriptors.ExactMolWt(molecula)
            n_hba = Descriptors.NumHAcceptors(molecula)
            n_hbd = Descriptors.NumHDonors(molecula)
            logp = Descriptors.MolLogP(molecula)
            file_m = dict[molecula]

            a.write(file_m + " descriptors:" + "\n" )
            a.write("molecular weight: " + str(molecular_weight) + "\n")
            a.write("n_hba: " + str(n_hba) + "\n")
            a.write("n_hbd: " + str(n_hbd) + "\n")
            a.write("logp: " + str(logp) + "\n")
            a.write("\n")


if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    rule = paramsDic['rule']
    molFiles = paramsDic['ligandFiles']
    molsDict, mols = getMolFilesDic(molFiles)
    if rule in ['l5', 'ro5']:
        results = ro5Filter(mols)
    elif 'ro3' in rule:
        results = ro3Filter(mols)

    with open(paramsDic['outputPath'], 'w') as f:
        f.write("#The following molecules have been passed the " + str(rule) + " filter:" + "\n")
        results = sum(results, [])
        for molecule in results:
            file_f = molsDict[molecule]
            f.write(str(file_f) + "\n")

    writeDetails(results, molsDict)

