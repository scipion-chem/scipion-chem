# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# # # *
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
import os, sys, copy
from operator import itemgetter

from rdkit import RDConfig, Chem, Geometry, DistanceGeometry
from rdkit.Chem import ChemicalFeatures, rdDistGeom, Draw, rdMolTransforms, AllChem
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
from rdkit.Numerics import rdAlignment

from utils import getMolFilesDic, parseParams, writeMol, getBaseName

ABSOLUTE, PROP_MOLS, PROP_FEATS = ['Absolute', 'Molecules proportion', 'Features proportion']
DBSCAN, KMEANS = ['DBSCAN', 'KMeans']


def buildPharmacophore(pharmDic):
    chemFeats, radii = [], []
    for featId in pharmDic:
        chemFeats.append(ChemicalFeatures.FreeChemicalFeature(pharmDic[featId]['type'],
                                                              Geometry.Point3D(*pharmDic[featId]['coords'])))
        radii.append(pharmDic[featId]['radius'])

    pcophore = Pharmacophore.Pharmacophore(chemFeats)
    for i in range(len(radii)):
        for j in range(i + 1, len(radii)):
            sumRadii = radii[i] + radii[j]
            pcophore.setLowerBound(i, j, max(pcophore.getLowerBound(i, j) - sumRadii, 0))
            pcophore.setUpperBound(i, j, pcophore.getUpperBound(i, j) + sumRadii)

    return pcophore

def getTransformMatrix(alignRef, confEmbed, atomMatch):
  alignProbe = []
  for matchIds in atomMatch:
    dummyPoint = Geometry.Point3D(0.0,0.0,0.0)
    for id in matchIds:
      dummyPoint += confEmbed.GetAtomPosition(id)
    dummyPoint /= len(matchIds)
    alignProbe.append(dummyPoint)
  return rdAlignment.GetAlignmentTransform(alignRef,alignProbe)

def transformEmbeddings(pcophore, embeddings, atomMatch):
  alignRef = [f.GetPos() for f in pcophore.getFeatures()]
  SSDs, confs = [], []
  for embedding in embeddings:
    conf = embedding.GetConformer()
    #print(conf.GetPositions())
    SSD, transformMatrix = getTransformMatrix(alignRef, conf, atomMatch)
    rdMolTransforms.TransformConformer(conf, transformMatrix)
    SSDs.append(SSD)
    confs.append(conf)
  return SSDs, confs

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> <outputDir>
    '''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles', 'moleculesFiles'], sep='::')
    it = sys.argv[1].split('_')[-1].split('.')[0]
    outDir = sys.argv[2]
    ligandFiles = paramsDic['ligandFiles']
    pharmDic = eval(paramsDic['pharmDic'])

    downSample = eval(paramsDic['downSample'])
    optimize = eval(paramsDic['optimize'])
    nAlignments = int(paramsDic['nAlignments'])
    maxSSD = float(paramsDic['maxSSD'])


#####################################################################
    # Load molecules and pharmacophore to RDKit objects
    dict_molecules = {}
    molFileDic, mols = getMolFilesDic(ligandFiles)
    if len(mols) > 0:
        pcophore = buildPharmacophore(pharmDic)

        ###################### Check for pharmacophore matches #######################
        fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
        featFactory = ChemicalFeatures.BuildFeatureFactory(fdefName)

        res = {}
        for i, mol in enumerate(mols):
            try:
                boundsMat = rdDistGeom.GetMoleculeBoundsMatrix(mol)
                canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(mol, featFactory, pcophore)
            except:
                pass
            if canMatch:
                # for (i, match) in enumerate(allMatches):
                #     for f in match:
                #         print("%d %s %s %s" % (i, f.GetFamily(), f.GetType(), f.GetAtomIds()))
                failed, boundsMatMatched, matched, matchDetails = \
                    EmbedLib.MatchPharmacophore(allMatches, boundsMat, pcophore, useDownsampling=downSample)
                if failed:
                    print("Couldn't embed molecule {}".format(molFileDic[mol]))
                    continue
            else:
                print("Couldn't match molecule {}".format(molFileDic[mol]))
                continue
            atomMatch = [list(x.GetAtomIds()) for x in matched]
            # for match in matched:
            #     print("%s %d" % (match.GetFamily(), match.GetAtomIds()[0]))

            try:
                molH = Chem.AddHs(mol, addCoords=True)
                Chem.AssignStereochemistry(molH, force=True, cleanIt=True)
                bm, embeddings, numFail = EmbedLib.EmbedPharmacophore(molH, atomMatch, pcophore, randomSeed=44, silent=True,
                                                                      targetNumber=nAlignments, count=nAlignments*10)
                # always yielding equal embeddings for all nAlignments
                if optimize:
                    optEmbeddings = []
                    for embMol in embeddings:
                        optMol = copy.deepcopy(embMol)
                        e1, e2 = EmbedLib.OptimizeMol(optMol, bm, atomMatch)
                        optEmbeddings.append(optMol)
                else:
                    optEmbeddings = embeddings

            except:
                print("Bounds smoothing failed for molecule {}".format(molFileDic[mol]))
                continue


            SSDs, confs = transformEmbeddings(pcophore, optEmbeddings, atomMatch)
            if SSDs:
                bestFitIndex = min(enumerate(SSDs), key=itemgetter(1))[0]

                #print('Best: ', bestFitIndex)
                #print('SSDs: ', SSDs)
                prevSSDs = []
                for i, conf in enumerate(confs):
                    if maxSSD >= SSDs[i] and not SSDs[i] in prevSSDs:
                        optMol = optEmbeddings[i]
                        outFile = os.path.join(outDir, getBaseName(molFileDic[mol]) + '_{}.sdf'.format(i))
                        writeMol(optMol, outFile, cid=confs[i].GetId())
                        prevSSDs.append(SSDs[i])

                        res[mol] = [outFile, SSDs[i]]

        with open(os.path.join(outDir, f'deviations_{it}.tsv'), 'w') as f:
            f.write('OriginalMolecule\tOutputMolecule\tDeviation\n')
            for mol in res:
                f.write('{}\t{}\t{}\n'.format(getBaseName(molFileDic[mol]), res[mol][0], res[mol][1]))

    else:
        print('None of the input molecules could be read by RDKit.\n'
              'Preparing the molecules with the RDKit ligand preparation protocol might solve this issue')


