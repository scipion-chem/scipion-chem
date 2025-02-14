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
import os, sys, ast, pickle, operator
from collections import OrderedDict, Counter

import numpy as np

from sklearn import cluster

from rdkit import RDConfig, Chem, Geometry, DistanceGeometry
from rdkit.Chem import ChemicalFeatures, rdDistGeom, Draw, rdMolTransforms, AllChem
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib

from utils import getMolFilesDic, parseParams, getBaseName

ABSOLUTE, PROP_MOLS, PROP_FEATS = ['Absolute', 'Molecules proportion', 'Features proportion']
DBSCAN, KMEANS = ['DBSCAN', 'KMeans']


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

def get_cluster_center_radii(feature_coords, cluster_labels, cluster_indices_sel, propRadii):
    cluster_centers, clusters_radii = [], []
    if cluster_indices_sel:
        for label in cluster_indices_sel:
            points_of_cluster = np.array(feature_coords)[cluster_labels == label]
            centroid_of_cluster = np.mean(points_of_cluster, axis=0)
            cluster_centers.append(centroid_of_cluster)

            dists = []
            for point in points_of_cluster:
                dists.append(np.sqrt(np.sum((centroid_of_cluster - point) ** 2, axis=0)))

            radii = sorted(dists)[int(propRadii*len(dists))]
            if radii < 1:
                radii = 1
            clusters_radii.append(radii)

    return cluster_centers, clusters_radii

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
    '''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles', 'moleculesFiles'], sep='::')
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

    propRadii = float(paramsDic['propRadii'])

#####################################################################
    # Load molecules and assign bonds order based on the SMILES
    dict_molecules = {}
    molsDic, mols = getMolFilesDic(paths_ligands)
    for mol in molsDic:
        dict_molecules[getBaseName(molsDic[mol])] = mol

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

    ## Filters top clusters with minimum size
    cluster_indices = {feat: [] for feat in selFeatures}
    for feat in cluster_indices:
        if clustDic[feat] is not None:
            min_cluster_size = getMinClusterSize(paramsDic, len(molecules), len(clustDic[feat]))
            selClusters = get_clusters(clustDic[feat], min_cluster_size, top_cluster_number)
            if selClusters:
                cluster_indices[feat] = selClusters

    # print('cluster_indices_sel: ', cluster_indices_sel)

    # Save coordinates of the cluster centers
    cluster_centers, cluster_radii = {}, {}
    for feat in cluster_indices:
        cluster_centers[feat], cluster_radii[feat] = get_cluster_center_radii(features_coord[feat], clustDic[feat],
                                                                              cluster_indices[feat], propRadii)

    #print('cluster_centers_sel: ', cluster_centers_sel)
    with open(os.path.join(outDir, 'cluster_centers.pkl'), 'wb') as outp:
        pickle.dump(cluster_centers, outp, pickle.HIGHEST_PROTOCOL)

    with open(os.path.join(outDir, 'cluster_radii.pkl'), 'wb') as outp:
        pickle.dump(cluster_radii, outp, pickle.HIGHEST_PROTOCOL)

