import operator, sys
import numpy as np
from collections import defaultdict

from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator, MACCSkeys
from rdkit.Avalon import pyAvalonTools

from rdkit.Chem import DataStructs

from rdkit.ML.Cluster import Butina
from rdkit.SimDivFilters import rdSimDivPickers
from scipy.spatial.distance import squareform
from sklearn.cluster import HDBSCAN, DBSCAN, Birch
from sklearn.metrics import pairwise_distances

from utils import getMolFilesDic, parseParams, typeParamsDic

CLUSTERINGS = ['butina', 'dbscan', 'hdbscan', 'birch', 'bitbirch', 'kmedoids']
PICKERS = ['leaderpick', 'maxminpick']

distDic = {'tanimoto': DataStructs.cDataStructs.BulkTanimotoSimilarity,
           'dice': DataStructs.cDataStructs.BulkDiceSimilarity,
           'cosine': DataStructs.cDataStructs.BulkCosineSimilarity}

def labelsToClusterIndex(labels):
    unique_labels = np.unique(labels)
    cluster_groups = {label: np.where(labels == label)[0].tolist() for label in unique_labels if label != -1}
    return list(cluster_groups.values())

def assignPointsToClusters(picks, fps, mols, distType='tanimoto'):
    clusters = defaultdict(list)
    for i, idx in enumerate(picks):
        clusters[i].append(mols[idx])
    sims = np.zeros((len(picks), len(fps)))
    for i in range(len(picks)):
        pick = picks[i]
        sims[i, :] = distDic[distType.lower()](fps[pick], fps)
        sims[i, i] = 0
    best = np.argmax(sims, axis=0)
    for i, idx in enumerate(best):
      if i not in picks:
        clusters[idx].append(mols[i])
    return list(clusters.values())

#################################################################################################################
def makeFingerprints(mols, fingerprint, paramsDic):
    fing = fingerprint.lower()

    if fing in ['morgan', 'rdkit', 'atompair', 'topologicaltorsion']:
        commonKwargs = {'fpSize': paramsDic["fingerSize"]}
        if fing == 'morgan':
            fingKwargs = {"radius": paramsDic["radius"], 'includeChirality': paramsDic["useChiralty"]}
            fingGen = rdFingerprintGenerator.GetMorganGenerator(**commonKwargs, **fingKwargs)

        elif fing == 'rdkit':
            fingKwargs = {"minPath": paramsDic["minPath"], 'maxPath': paramsDic["maxPath"],
                          'useHs': paramsDic["useHs"]}
            fingGen = rdFingerprintGenerator.GetRDKitFPGenerator(**commonKwargs, **fingKwargs)

        elif fing == 'atompair':
            fingKwargs = {"minDistance": paramsDic["minDistance"], 'use2D': paramsDic["minDistance"],
                          'includeChirality': paramsDic["useChiralty"]}
            if int(paramsDic["maxDistance"]) > 0:
                fingKwargs.update({"maxDistance": paramsDic["maxDistance"]})
            fingGen = rdFingerprintGenerator.GetAtomPairGenerator(**commonKwargs, **fingKwargs)

        elif fing == 'topologicaltorsion':
            fingKwargs = {'includeChirality': paramsDic["useChiralty"]}
            fingGen = rdFingerprintGenerator.GetTopologicalTorsionGenerator(**commonKwargs, **fingKwargs)

        fingerprints = [fingGen.GetFingerprint(m) for m in mols]

    elif fing == 'maccs':
        fingerprints = [MACCSkeys.GenMACCSKeys(x) for x in mols]

    elif fing == 'pattern':
        fingerprints = [Chem.PatternFingerprint(x) for x in mols]

    elif fing == 'layered':
        fingerprints = [Chem.LayeredFingerprint(x) for x in mols]

    elif fing == 'avalon':
        fingerprints = [pyAvalonTools.GetAvalonFP(x) for x in mols]

    return fingerprints

def calculateDistances(fingerprints, distType='Tanimoto'):
    distType = distType.lower()

    distFunc = distDic[distType]
    dists = []
    for i in range(1, len(fingerprints)):
        dist = distFunc(fps[i], fps[:i], returnDistance=True)
        dists.extend([x for x in dist])

    return dists

def getMedoidsFromMatrix(cluster_indices, distance_matrix):
    medoids = []
    # Get indices of points in the current cluster
    for clIndex in cluster_indices:
        # Subset the distance matrix for the current cluster
        cluster_distances = distance_matrix[np.ix_(clIndex, clIndex)]

        # Sum of distances for each point in the cluster
        sum_distances = cluster_distances.sum(axis=1)

        # The medoid is the point with the smallest sum of distances
        medoid_index = clIndex[np.argmin(sum_distances)]
        medoids.append(medoid_index)

    return medoids


def getMedoidsFromCluster(fps, cluster_indices, distance='jaccard'):
    distance = 'jaccard' if distance.lower() == 'tanimoto' else distance
    medoids = []
    for clIndex in cluster_indices:
        # Get indices of points in this cluster
        cluster_points = fps[clIndex]
        # Compute pairwise distances within the cluster
        distances = pairwise_distances(cluster_points, metric=distance)

        # Sum of distances for each point (smallest = medoid)
        sum_distances = distances.sum(axis=1)
        medoid_idx = clIndex[np.argmin(sum_distances)]
        medoids.append(medoid_idx)

    return medoids

def buildClusters(fps, mols, clustType, paramsDic):
    clustType = clustType.lower()

    if clustType in CLUSTERINGS:
        if clustType not in ['birch', 'bitbirch']:
            dists = calculateDistances(fps, distType=paramsDic['distance'])
            distMatrix = squareform(dists)

        if clustType == 'butina':
            cluster_indices = Butina.ClusterData(dists, len(mols), paramsDic['cutoff'], isDistData=True)

        elif clustType == 'dbscan':
            db = DBSCAN(eps=paramsDic['cutoff'], min_samples=paramsDic['minSamples'], metric='precomputed')
            db.fit(distMatrix)
            cluster_indices = labelsToClusterIndex(db.labels_)

        elif clustType == 'hdbscan':
            hdb = HDBSCAN(min_cluster_size=paramsDic['minClusterSize'], metric='precomputed')
            hdb.fit(distMatrix)
            cluster_indices = labelsToClusterIndex(hdb.labels_)

        elif clustType == 'kmedoids':
            from sklearn_extra.cluster import KMedoids

            kmed = KMedoids(n_clusters=paramsDic['nClusters'], metric='precomputed')
            kmed.fit(distMatrix)
            cluster_indices = labelsToClusterIndex(kmed.labels_)

        if clustType == 'birch':
            bir = Birch(threshold=paramsDic['cutoff'], branching_factor=paramsDic['branchingFactor'], 
                        n_clusters=paramsDic['nClusters'])
            bir.fit(np.array(fps))
            cluster_indices = labelsToClusterIndex(bir.labels_)

        elif clustType == 'bitbirch':
            import bitbirch.bitbirch as bb
            fps = np.array(fps)

            bb.set_merge('diameter')
            brc = bb.BitBirch(threshold=paramsDic['cutoff'], branching_factor=paramsDic['branchingFactor'])
            brc.fit(fps)
            cluster_indices = brc.get_cluster_mol_ids()

        cluster_mols = [operator.itemgetter(*cluster)(mols) for cluster in cluster_indices]

        if clustType in ['butina', 'dbscan', 'hdbscan']:
            reps = getMedoidsFromMatrix(cluster_indices, distMatrix)
        elif clustType in ['birch', 'bitbirch']:
            reps = getMedoidsFromCluster(fps, cluster_indices, paramsDic['distance'])
        elif clustType == 'kmedoids':
            reps = kmed.medoid_indices_

    else:
        # this pickers do not create clusters, directly pick the representatives
        if clustType == 'leaderpick':
            lp = rdSimDivPickers.LeaderPicker()
            picks = lp.LazyBitVectorPick(fps, len(fps), paramsDic['cutoff'])

            cluster_mols = assignPointsToClusters(picks, fps, mols, paramsDic['distance'])

        elif clustType == 'maxminpick':
            lp = rdSimDivPickers.MaxMinPicker()
            picks = lp.LazyBitVectorPick(fps, len(fps), paramsDic['nClusters'])

            cluster_mols = assignPointsToClusters(picks, fps, mols, paramsDic['distance'])

        reps = [p for p in picks]

    rep_mols = [mols[i] for i in reps]
    cluster_mols = [[c] if isinstance(c, Chem.rdchem.Mol) else c for c in cluster_mols]
    return cluster_mols, rep_mols

###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['molFiles'])
    molFiles = paramsDic['molFiles']
    outputPath = paramsDic['outputPath']

    typeDic = {int: ['fingerSize', 'radius', 'minPath', 'maxPath', 'minDistance', 'maxDistance', 'branchingFactor',
                     'nClusters', 'minSamples', 'minClusterSize'],
               float: ['cutoff'],
               bool: ['useChiralty', 'useHs', 'use2D']}
    paramsDic = typeParamsDic(paramsDic, typeDic)

    mols_dict, mols = getMolFilesDic(molFiles)
    fps = makeFingerprints(mols, paramsDic['finger'], paramsDic)

    cluster_mols, rep_mols = buildClusters(fps, mols, paramsDic['cluster'], paramsDic)
    sizeClusters = [len(c) for c in cluster_mols]
    print(f'n clusters: {len(cluster_mols)}, mean size: {sum(sizeClusters)/len(sizeClusters)}')

    with open(outputPath, 'w') as f:
        for i, cluster in enumerate(cluster_mols):
            for mol in cluster:
                isRep = mol in rep_mols
                f.write(f'{mols_dict[mol]}\t{i}\t{isRep}\n')
