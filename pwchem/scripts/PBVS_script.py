# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *
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
import os, sys, math, ast, pickle, operator
from collections import OrderedDict, Counter

#import pandas as pd
import numpy as np

from sklearn import datasets, cluster

from rdkit import RDConfig, Chem, Geometry, DistanceGeometry
from rdkit.Chem import ChemicalFeatures, rdDistGeom, Draw, rdMolTransforms, AllChem
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib


ABSOLUTE, PROP_MOLS, PROP_FEATS = ['Absolute', 'Molecules proportion', 'Features proportion']
DBSCAN, KMEANS = ['DBSCAN', 'KMeans']

def getBaseFileName(filename):
    return os.path.splitext(os.path.basename(filename))[0]

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
            key, value = line.strip().split('::')
            if key == 'ligandFiles' or key == 'moleculesFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic

def clusteringDBScan(feature_coord, dbEps):
    clus = None
    if feature_coord:
        model = cluster.DBSCAN(min_samples=1, eps=dbEps)
        clus = model.fit_predict(np.array(feature_coord))

    return clus

def clusteringKMeans(feature_coord, kNumber):
    k_means = None
    if feature_coord:
        if kNumber > len(feature_coord):
            kNumber = len(feature_coord)
        #print(f"Clustering: \nVariable k in k-means: {kNumber} of {len(feature_coord)} points\n")

        # Initialize k-means and compute clustering
        k_means = cluster.KMeans(n_clusters=kNumber, random_state=44)
        k_means.fit(feature_coord)

    return k_means

def get_clusters(feature_labels, min_cluster_size, top_cluster_number):
    if not feature_labels is None:
        feature_labels_count = Counter(feature_labels)

        feature_labels_count = sorted(
            feature_labels_count.items(), key=operator.itemgetter(1), reverse=True
        )
        # print(f"Clusters sorted by size: \n{feature_labels_count}\n")

        # Get number of the largest clusters, which are larger then the threshold (selected clusters)
        cluster_indices_sel = []

        for cluster_index, cluster_size in feature_labels_count:
            if cluster_size >= min_cluster_size and \
                    (top_cluster_number < 0 or  len(cluster_indices_sel) < top_cluster_number):
                cluster_indices_sel.append(cluster_index)

        # print(f"Cluster indices of selected clusters: \n{cluster_indices_sel}\n")

    else:
        cluster_indices_sel = None

    return cluster_indices_sel

def get_selected_cluster_center_coords(feature_coords, cluster_labels, cluster_indices_sel):
    cluster_centers_sel = []
    if cluster_indices_sel:
        for label in cluster_indices_sel:
            points_of_cluster = np.array(feature_coords)[cluster_labels == label]
            centroid_of_cluster = np.mean(points_of_cluster, axis=0)
            cluster_centers_sel.append(centroid_of_cluster)
    else:
        pass

    return list(cluster_centers_sel)

def create_Pharmacophore(final_dict, feature):
    samplePh4Feats = []
    for key, coord in final_dict.items():
        if key == feature:
            samplePh4Feats.append(ChemicalFeatures.FreeChemicalFeature(str(key[:-1]), Geometry.Point3D(*coord)))
    pcophore = Pharmacophore.Pharmacophore(samplePh4Feats)
    return pcophore

def getMinClusterSize(paramsDic, nMols, nFeats):
    numberType = paramsDic['minType']
    number = float(paramsDic['minSize'])

    if numberType == PROP_MOLS:
        minClusterSize = int(nMols * number)
    elif numberType == ABSOLUTE:
        minClusterSize = number
    elif numberType == PROP_FEATS:
        minClusterSize = int(nFeats * number)
    return minClusterSize


if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> <outputDir>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    outDir = sys.argv[2]
    paths_ligands = paramsDic['ligandFiles']
    dict_smiles = ast.literal_eval(paramsDic['ligandSmiles'])

    selFeatures = paramsDic['features'].split()

    method = paramsDic['method']
    if method == KMEANS:
        kNumber = int(paramsDic['kNumber'])
    elif method == DBSCAN:
        dbEps = float(paramsDic['eps'])

    top_cluster_number = int(paramsDic['topClusters'])

#####################################################################
    # Load molecules and assign bonds order based on the SMILES
    dict_molecules = {}
    molsDic, mols = preprocessLigands(paths_ligands)
    for mol in molsDic:
        dict_molecules[getBaseFileName(molsDic[mol])] = mol

    with open(os.path.join(outDir, 'molecules.pkl'), 'wb') as outp:
        pickle.dump(mols, outp, pickle.HIGHEST_PROTOCOL)

    # Building dictionary: {ligBase: [ChemMol(from file), ChemMol(from SMILES)]}
    dict_mol_reference = {}
    for key, value in dict_molecules.items():
        if key in dict_smiles.keys():
            smiles_mol = Chem.MolFromSmiles(dict_smiles[key])
            dict_mol_reference[key] = [dict_molecules[key], smiles_mol]

    dict_mol_reference = OrderedDict(sorted(dict_mol_reference.items()))

    molecules = []
    for key in dict_mol_reference:
        molecule, reference_molecule = dict_mol_reference[key]
        try:
            molecules += [AllChem.AssignBondOrdersFromTemplate(reference_molecule, molecule)]
        except:
            print('AssignBondOrdersFromTemplate failed for ', key)

    ###################### BUILD RDKIT PHARMACOPHORE #######################
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    # molecule_feature_frequencies = []
    # for molecule2 in molecules_2:
    #     features = [feature.GetFamily() for feature in factory.GetFeaturesForMol(molecule2)]
    #     feature_frequency = Counter(features)
    #     molecule_feature_frequencies.append(feature_frequency)
    #
    # feature_frequencies_df = (
    #     pd.DataFrame(molecule_feature_frequencies, index=[f"Mol{i}" for i, _ in enumerate(molecules_2, 1)],).fillna(0).astype(int))
    # print(feature_frequencies_df.transpose())


    ########################################
    # Build features for each molecule and retrieve coordinates
    features = {feat: [] for feat in selFeatures}
    for feat in features:
        for mol in molecules:
            features[feat].append(factory.GetFeaturesForMol(mol, includeOnly=feat))
    #print('features: ', features)

    features_coord = {}
    for feat in features:
        features_coord[feat] = [list(item.GetPos()) for sublist in features[feat] for item in sublist]
    #print('features_coord: ', features_coord)


    ############################################
    ## Clustering the feature coordinates using KMeans or DBScan
    clustDic = {}
    for feat in features_coord:
        if method == KMEANS:
            clustDic[feat] = clusteringKMeans(features_coord[feat], kNumber).labels_
        elif method == DBSCAN:
            clustDic[feat] = clusteringDBScan(features_coord[feat], dbEps)
   # print('clustDic: ', clustDic)

    ## Filters top clusters with minimum size
    cluster_indices_sel = {feat: [] for feat in selFeatures}
    for feat in cluster_indices_sel:
        min_cluster_size = getMinClusterSize(paramsDic, len(molecules), len(clustDic[feat]))
        selClusters = get_clusters(clustDic[feat], min_cluster_size, top_cluster_number)
        if selClusters:
            cluster_indices_sel[feat] = selClusters

   # print('cluster_indices_sel: ', cluster_indices_sel)
    with open(os.path.join(outDir, 'cluster_index.pkl'), 'wb') as outp:
        pickle.dump(cluster_indices_sel, outp, pickle.HIGHEST_PROTOCOL)

    # Save coordinates of the cluster centers
    cluster_centers_sel = {}
    for feat in cluster_indices_sel:
        center2 = str(feat).capitalize()
        cluster_centers_sel[center2] = get_selected_cluster_center_coords(features_coord[feat], clustDic[feat],
                                                                          cluster_indices_sel[feat])

    #print('cluster_centers_sel: ', cluster_centers_sel)
    with open(os.path.join(outDir, 'cluster_centers.pkl'), 'wb') as outp:
        pickle.dump(cluster_centers_sel, outp, pickle.HIGHEST_PROTOCOL)

