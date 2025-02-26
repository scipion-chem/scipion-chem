# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

'''Script to score a SetOfSmallMolecules docked using the Open Drug Discovery Toolkit
(ODDT, https://github.com/oddt/oddt). It must be launch with the conda environment with rdkit and oddt'''

import sys, os, shutil

import oddt
from oddt.scoring.descriptors import (autodock_vina_descriptor, fingerprints, oddt_vina_descriptor)
from oddt.scoring.functions import rfscore, nnscore, PLECscore

from utils import parseParams

def oddt_vina_score(paramsDic):
    """Score molecules docking using the internal oddt vina score"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    vina_scores = ['vina_affinity', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen',
                   # intra-molecular interactions
                   'vina_intra_gauss1', 'vina_intra_gauss2', 'vina_intra_repulsion', 'vina_intra_hydrophobic',
                   'vina_intra_hydrogen', 'vina_num_rotors'
                   ]

    print('Running ODDT Vina scoring')
    oddt_vina_results = oddt_vina_descriptor(protein=rec, vina_scores=vina_scores).build(mols)
    #We return only the score (first column, affinity, and negative to resemble an score)
    return -oddt_vina_results[:,0] , molFiles


def oddt_rfscore_score(paramsDic, v):
    """Score a list of ligands docked to a protein using RFScore"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    #descs = rfscore(version=v, protein=rec).descriptor_generator.build(mols)
    #print('Descriptors: ', descs)

    if not 'model' in paramsDic:
        print('Training RFScore model: ', paramsDic['saveModel'])
        rfScoreModel = rfscore(version=v, spr=int(paramsDic['spr']))
        rfScoreModel.train(pdbbind_version=int(paramsDic['pdbbind']), sf_pickle=paramsDic['saveModel'])
    else:
        print('Loading RFScore model: ', paramsDic['model'])
        rfScoreModel = rfscore().load(paramsDic['model'])
    rfScoreModel.set_protein(rec)
    return rfScoreModel.predict(mols), molFiles


def oddt_nnscore_score(paramsDic):
    """Score a list of ligands docked to a protein using NNScore"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    #descs = nnscore(protein=rec).descriptor_generator.build(mols)
    # print('Descriptors: ', descs)


    if not 'model' in paramsDic:
        print('Training NNScore model: ', paramsDic['saveModel'])
        nnscore().train(pdbbind_version=int(paramsDic['pdbbind']), sf_pickle=paramsDic['saveModel'])
        nnScoreModel = nnscore().load(paramsDic['saveModel'])
    else:
        print('Loading NNScore model: ', paramsDic['model'])
        nnScoreModel = nnscore().load(paramsDic['model'])
    nnScoreModel.set_protein(rec)
    return nnScoreModel.predict(mols), molFiles


def oddt_PLECscore_score(paramsDic, v):
    """Score a list of ligands docked to a protein using PLECScore"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    #descs = PLECscore(version=v, protein=rec).descriptor_generator.build(mols)
    # print('Descriptors: ', descs)

    if not 'model' in paramsDic:
        print('Training PLECScore model: ', paramsDic['saveModel'])
        PLECScoreModel = PLECscore(version=v, depth_protein=int(paramsDic['depthProt']),
                                   depth_ligand=int(paramsDic['depthLig']), size=int(paramsDic['fingerSize']))
        PLECScoreModel.train(pdbbind_version=int(paramsDic['pdbbind']))
        PLECScoreModel.save(paramsDic['saveModel'])
    else:
        print('Loading PLECScore model: ', paramsDic['model'])
        PLECScoreModel = PLECscore().load(paramsDic['model'])

    PLECScoreModel.set_protein(rec)
    return PLECScoreModel.predict(mols), molFiles

def cleanPDB(pdbFile):
    auxFile = pdbFile.replace('.pdb', '_aux.pdb')
    with open(auxFile, 'w') as f:
        with open(pdbFile) as fIn:
            for line in fIn:
                if line.strip() != '' and \
                        line.split()[0] in ['HEADER', 'EXPDTA', 'REMARK', 'ATOM', 'HETATM', 'TER', 'CONECT']:
                    f.write(line)
    shutil.move(auxFile, pdbFile)
    return pdbFile

def preprocessLigands(ligandsFiles):
    mols = []
    for molFile in ligandsFiles:
        molExt = os.path.splitext(molFile)[1][1:]
        if molExt == 'pdb':
            molFile = cleanPDB(molFile)
        mols += list(oddt.toolkit.readfile(molExt, molFile))
    list(map(lambda x: x.addh(), mols))
    return mols

def preprocessReceptor(receptorFile):
    recExt = os.path.splitext(receptorFile)[1][1:]
    rec = next(oddt.toolkit.readfile(recExt, receptorFile))
    rec.protein = True
    rec.addh()
    return rec

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <Function_version> <PDBBind_train> <outputPath> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    function = paramsDic['function']
    if 'vina' in function.lower():
        results, mols = oddt_vina_score(paramsDic)
    elif 'rfscore' in function.lower():
        version = function.split('_')[1]
        results, mols = oddt_rfscore_score(paramsDic, int(version))
    elif 'nnscore' in function.lower():
        results, mols = oddt_nnscore_score(paramsDic)
    elif 'plecscore' in function.lower():
        version = function.split('_')[1]
        results, mols = oddt_PLECscore_score(paramsDic, version)

    with open(paramsDic['outputPath'], 'w') as f:
        for molFile, score in zip(mols, results):
            f.write('{}\t{}\n'.format(molFile, score))

