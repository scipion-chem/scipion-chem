from __future__ import print_function
import sys, os
import csv
from os.path import dirname, isfile, join
from multiprocessing import Pool
from itertools import chain
from functools import partial
import warnings
import oddt
from oddt.utils import is_molecule, compose_iter, chunker, method_caller
from oddt.scoring import scorer
from oddt.fingerprints import (InteractionFingerprint,
                               SimpleInteractionFingerprint,
                               dice)
from oddt.shape import usr, usr_cat, electroshape, usr_similarity
from oddt.scoring.descriptors import (autodock_vina_descriptor, fingerprints, oddt_vina_descriptor)
from oddt.scoring.functions import rfscore, nnscore, PLECscore
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw, PandasTools
import sys
import itertools


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
def _filter(expression, mols, soft_fail=0):
    """Filter molecule by a generic expression, such as `mol.logp > 1`."""

    # lista con las moleculas que pasan el filtro
    out = []
    for mol in mols:
        # saber si la expresion es una lista
        if isinstance(expression, list):
            fail = 0
            # si lo es entramos en el bucle
            for e in expression:
                # si la expresión es falsa
                if eval(e):
                    fail += 1

            if fail >= 3:
                out.append(mol)

        # este es en el caso de que la expresion no sea una lista (unica expresion)
        else:
            if eval(expression):
                out.append(mol)
    return out


def ro5_filter(paramsDic):
    filtered_mols = []
    molFiles = paramsDic['ligandFiles']
    mols_dict, mols = preprocessLigands(molFiles)

    out = _filter(mols=mols, expression=['Descriptors.ExactMolWt(mol) < 500',
                                         'Descriptors.NumHAcceptors(mol) <= 10',
                                         'Descriptors.NumHDonors(mol) <= 5',
                                         'Descriptors.MolLogP(mol) <= 5'], soft_fail=0)

    filtered_mols.append(out)
    return filtered_mols, mols_dict, mols


def ro3_filter(paramsDic):
    filtered_mols = []
    molFiles = paramsDic['ligandFiles']
    mols_dict, mols = preprocessLigands(molFiles)

    out = _filter(mols=mols, expression=['Descriptors.ExactMolWt(mol) < 300',
                                         'Descriptors.NumHAcceptors(mol) <= 3',
                                         'Descriptors.NumHDonors(mol) <= 3',
                                         'Descriptors.MolLogP(mol) <= 3'], soft_fail=0)

    filtered_mols.append(out)
    return filtered_mols, mols_dict, mols

###############################################################################
def comprobar(moleculas, dict):
    with open("comprobar.txt", 'w') as a:
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
    paramsDic = parseParams(sys.argv[1])
    rule = paramsDic['rule']
    if rule in ['l5', 'ro5']:
        results, mols_dict, mols = ro5_filter(paramsDic)
    elif 'ro3' in rule:
        results, mols_dict, mols = ro3_filter(paramsDic)

    print("mols")
    print(mols_dict)

    with open(paramsDic['outputPath'], 'w') as f:
        f.write("#The following molecules have been passed the " + str(rule) + " filter:" + "\n")
        results = sum(results, [])
        for molecule in results:
            file_f = mols_dict[molecule]
            f.write(str(file_f) + "\n")

    comprobar(results, mols_dict)

