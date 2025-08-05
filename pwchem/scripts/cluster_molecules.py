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
    uniqueLabels = np.unique(labels)
    clusterGroups = {label: np.where(labels == label)[0].tolist() for label in uniqueLabels if label != -1}
    return list(clusterGroups.values())

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
        dists.extend(dist)

    return dists

def getMedoidsFromMatrix(clusterIndices, distanceMatrix):
    medoids = []
    # Get indices of points in the current cluster
    for clIndex in clusterIndices:
        # Subset the distance matrix for the current cluster
        clusterDistances = distanceMatrix[np.ix_(clIndex, clIndex)]

        # Sum of distances for each point in the cluster
        sumDistances = clusterDistances.sum(axis=1)

        # The medoid is the point with the smallest sum of distances
        medoidIndex = clIndex[np.argmin(sumDistances)]
        medoids.append(medoidIndex)

    return medoids


def getMedoidsFromCluster(fps, clusterIndices, distance='jaccard'):
    distance = 'jaccard' if distance.lower() == 'tanimoto' else distance
    medoids = []
    for clIndex in clusterIndices:
        # Get indices of points in this cluster
        clusterPoints = fps[clIndex]
        # Compute pairwise distances within the cluster
        distances = pairwise_distances(clusterPoints, metric=distance)

        # Sum of distances for each point (smallest = medoid)
        sumDistances = distances.sum(axis=1)
        medoidIdx = clIndex[np.argmin(sumDistances)]
        medoids.append(medoidIdx)

    return medoids

def optimizeParameter(
        bb, fps,  # Function to optimize (returns a value based on param)
        paramLow,  # Lowest possible parameter value
        paramHigh,  # Highest possible parameter value
        targetValue,  # Desired output of `func`
        step,  # Initial step size
        tolerance=1e-6,  # Stopping condition (how close to `targetValue`)
        maxIter=100,  # Prevent infinite recursion
        increasing=True,  # Whether `func` increases with parameter (True/False)
        verbose=False  # Print debug info
):
    """
    Recursively optimizes a parameter until `func(param)` approximates `targetValue`.
    Returns the best parameter value found.
    """
    bestParam = None
    bestDiff = float('inf')

    step = abs(step) if increasing else -abs(step)

    # Iterate through parameter values
    param = paramLow
    while (increasing and param <= paramHigh+1e-6) or (not increasing and param >= paramHigh-1e-6):
        brc = bb.BitBirch(threshold=param, branching_factor=paramsDic['branchingFactor'])
        brc.fit(fps)
        currentValue = len(brc.get_centroids())
        currentDiff = abs(currentValue - targetValue)

        if verbose:
            print(f"Trying param={param:.10f}, value={currentValue:.6f}, diff={currentDiff:.6f}")

        # Update best parameter if closer to target
        if currentDiff < bestDiff:
            bestDiff = currentDiff
            bestParam = param

        # Early exit if within tolerance or maxIter
        if currentDiff <= tolerance:
            return brc, bestParam, currentValue

        # Check if we passed the target (crossing point)
        if (increasing and currentValue > targetValue) or (not increasing and currentValue < targetValue):
            newLow = param - step
            newHigh = param
            newStep = step / 2  # Halve step

            if step <= 1e-6 or maxIter <= 1:
                return brc, bestParam, currentValue
            else:
                if verbose:
                    print(f"Crossed target! Refining between {newLow:.10f} and {newHigh:.10f} with step={newStep:.10f}")
                return optimizeParameter(
                    bb, fps, newLow, newHigh, targetValue, newStep,
                    tolerance, maxIter - 1, increasing, verbose
                )

        param += step

    # If no crossing found, return best param
    return brc, bestParam, currentValue


def buildClusters(fps, mols, clustType, paramsDic):
    clustType = clustType.lower()

    if clustType in CLUSTERINGS:
        if clustType not in ['birch', 'bitbirch']:
            dists = calculateDistances(fps, distType=paramsDic['distance'])
            distMatrix = squareform(dists)

        if clustType == 'butina':
            clusterIndices = Butina.ClusterData(dists, len(mols), paramsDic['cutoff'], isDistData=True)

        elif clustType == 'dbscan':
            db = DBSCAN(eps=paramsDic['cutoff'], min_samples=paramsDic['minSamples'], metric='precomputed')
            db.fit(distMatrix)
            clusterIndices = labelsToClusterIndex(db.labels_)

        elif clustType == 'hdbscan':
            hdb = HDBSCAN(min_cluster_size=paramsDic['minClusterSize'], metric='precomputed')
            hdb.fit(distMatrix)
            clusterIndices = labelsToClusterIndex(hdb.labels_)

        elif clustType == 'kmedoids':
            from sklearn_extra.cluster import KMedoids

            kmed = KMedoids(n_clusters=paramsDic['nClusters'], metric='precomputed')
            kmed.fit(distMatrix)
            clusterIndices = labelsToClusterIndex(kmed.labels_)

        if clustType == 'birch':
            bir = Birch(threshold=paramsDic['cutoff'], branching_factor=paramsDic['branchingFactor'], 
                        n_clusters=paramsDic['nClusters'])
            bir.fit(np.array(fps))
            clusterIndices = labelsToClusterIndex(bir.labels_)

        elif clustType == 'bitbirch':
            import bitbirch.bitbirch as bb
            fps = np.array(fps)

            bb.set_merge('diameter')
            highCutoff, lowCutoff = paramsDic['cutoff'], paramsDic['cutoffLow']
            objNClust, thStep = paramsDic['nClusters'], paramsDic['cutoffStep']

            brc, th, curClust = optimizeParameter(bb, fps, lowCutoff, highCutoff, objNClust, thStep, 1, 100, True, True)

            brc.reassign(fps, curClust)
            clusterIndices = brc.get_cluster_mol_ids()

        clusterMols = [operator.itemgetter(*cluster)(mols) for cluster in clusterIndices]

        if clustType in ['butina', 'dbscan', 'hdbscan']:
            reps = getMedoidsFromMatrix(clusterIndices, distMatrix)
        elif clustType in ['birch', 'bitbirch']:
            reps = getMedoidsFromCluster(fps, clusterIndices, paramsDic['distance'])
        elif clustType == 'kmedoids':
            reps = kmed.medoid_indices_

    else:
        # this pickers do not create clusters, directly pick the representatives
        if clustType == 'leaderpick':
            lp = rdSimDivPickers.LeaderPicker()
            reps = lp.LazyBitVectorPick(fps, len(fps), paramsDic['cutoff'])

            clusterMols = assignPointsToClusters(reps, fps, mols, paramsDic['distance'])

        elif clustType == 'maxminpick':
            lp = rdSimDivPickers.MaxMinPicker()
            reps = lp.LazyBitVectorPick(fps, len(fps), paramsDic['nClusters'])

            clusterMols = assignPointsToClusters(reps, fps, mols, paramsDic['distance'])

    repMols = [mols[i] for i in reps]
    clusterMols = [[c] if isinstance(c, Chem.rdchem.Mol) else c for c in clusterMols]
    return clusterMols, repMols

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
               float: ['cutoff', 'cutoffLow', 'cutoffStep'],
               bool: ['useChiralty', 'useHs', 'use2D']}
    paramsDic = typeParamsDic(paramsDic, typeDic)

    molsDict, mols = getMolFilesDic(molFiles)
    fps = makeFingerprints(mols, paramsDic['finger'], paramsDic)

    clusterMols, repMols = buildClusters(fps, mols, paramsDic['cluster'], paramsDic)
    sizeClusters = [len(c) for c in clusterMols]
    #print(f'n clusters: {len(clusterMols)}, mean size: {sum(sizeClusters)/len(sizeClusters)}')

    with open(outputPath, 'w') as f:
        for i, cluster in enumerate(clusterMols):
            for mol in cluster:
                isRep = mol in repMols
                f.write(f'{molsDict[mol]}\t{i}\t{isRep}\n')
