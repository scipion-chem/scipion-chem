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

'''Script to describe a SetOfSmallMolecules docked using the Open Drug Discovery Toolkit
(ODDT, https://github.com/oddt/oddt). It must be launch with the conda environment with rdkit and oddt'''

import sys, os
import numpy as np
from functools import partial

import oddt
from oddt.scoring.descriptors import (fingerprints, oddt_vina_descriptor,
                                      close_contacts_descriptor, universal_descriptor)
from oddt.scoring.functions import rfscore, nnscore, PLECscore
from oddt.fingerprints import *

from utils import parseParams

fpDic = {'PLEC': PLEC, 'SimpleInteractionFingerprint': SimpleInteractionFingerprint, 'ECFP': ECFP, 'SPLIF': SPLIF}

def oddt_vina_descriptors(paramsDic):
    """Return oddt vina descriptors"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    vina_scores = ['vina_affinity', 'vina_gauss1', 'vina_gauss2', 'vina_repulsion', 'vina_hydrophobic', 'vina_hydrogen',
                   'vina_intra_gauss1', 'vina_intra_gauss2', 'vina_intra_repulsion', 'vina_intra_hydrophobic',
                   'vina_intra_hydrogen', 'vina_num_rotors']

    descs = oddt_vina_descriptor(protein=rec, vina_scores=vina_scores)
    oddt_vina_results = descs.build(mols)
    oddt_vina_results = np.vstack((descs.titles, oddt_vina_results))
    return oddt_vina_results, molFiles


def oddt_rfscore_descriptors(paramsDic):
    """Return close contact descriptors"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])

    rfModel = rfscore(protein=rec, version=int(paramsDic['version']))
    descs = rfModel.descriptor_generator
    results = descs.build(mols)
    results = np.vstack((descs.titles, results))
    return results, molFiles

def oddt_nnscore_descriptors(paramsDic):
    """Score a list of ligands docked to a protein using PLECScore"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])

    nnModel = nnscore(protein=rec)
    descs = nnModel.descriptor_generator
    results = descs.build(mols)
    results = np.vstack((descs.titles, results))
    return results, molFiles

def oddt_fingerprints(paramsDic, fp='PLEC'):
    """Returns fingerprint descriptors"""
    molFiles = paramsDic['ligandFiles']
    mols = preprocessLigands(molFiles)
    rec = preprocessReceptor(paramsDic['receptorFile'])
    size = paramsDic['size'] if 'size' in paramsDic else None

    if fp == 'PLEC':
        fp_func = partial(fpDic[fp], depth_ligand=paramsDic['depth_ligand'], depth_protein=paramsDic['depth_protein'],
                          size=size, count_bits=True, sparse=True, ignore_hoh=True)
    elif 'SimpleInteractionFingerprint' in fp:
        size=168
        fp_func = partial(fpDic[fp], strict=True)

    # elif fp == 'SPLIF':
    #     fp_func = partial(fpDic[fp], depth=1, size=4096, distance_cutoff=4.5)
    #
    # elif fp == 'ECFP':
    #     fp_func = partial(fpDic[fp], depth=1, size=4096, count_bits=True, sparse=True,
    #                       use_pharm_features=False)

    descs = universal_descriptor(fp_func, protein=rec, shape=size, sparse=True)
    results = descs.build(mols)
    results = results.todense()
    results = np.vstack(([i for i in range(size)], results))
    return results, molFiles


############################### UTILS ############################

def preprocessLigands(ligandsFiles):
    mols = []
    for molFile in ligandsFiles:
        molExt = os.path.splitext(molFile)[1][1:]
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
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    descriptor = paramsDic['descriptor']
    if 'Vina' in descriptor:
        results, mols = oddt_vina_descriptors(paramsDic)
    elif 'RFScore' in descriptor:
        results, mols = oddt_rfscore_descriptors(paramsDic)
    elif 'NNScore' in descriptor:
        results, mols = oddt_nnscore_descriptors(paramsDic)
    elif 'PLECScore' in descriptor:
        newParams = {'depth_protein': 5, 'depth_ligand': 1, 'size': 65536}
        paramsDic.update(newParams)
        results, mols = oddt_fingerprints(paramsDic, fp='PLEC')
    elif 'Fingerprint' in descriptor:
        results, mols = oddt_fingerprints(paramsDic, fp=paramsDic['version'])

    mols = ['Mol file'] + mols
    with open(paramsDic['outputPath'], 'w') as f:
        for molFile, descriptors in zip(mols, results):
            if not descriptor in ['PLECScore', 'Fingerprint']:
                f.write('{}\t{}\n'.format(molFile, '\t'.join(str(x) for x in descriptors)))
            else:
                f.write('{}\t{}\n'.format(molFile, '\t'.join(str(x) for x in descriptors.flat)))

