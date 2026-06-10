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


"""
This protocol is used to make a consensus over several input sets of pockets that might come from different
sources

"""

import os
import numpy as np
from Bio import pairwise2, AlignIO
from Bio.Align import substitution_matrices
from Bio.PDB import PDBParser, MMCIFParser
from scipy.optimize import linear_sum_assignment

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput, StructROI
from pwchem.utils import writePDBLine, splitPDBLine, flipDic, createPocketFile, getBaseName, invertDic
from pwchem.utils.utilsFasta import getMultipleAlignmentCline

import networkx as nx
from typing import List, Dict
from collections import Counter


class ResidueGroup:
    """A group of residues, each with a centroid coordinate."""

    def __init__(self, residueIds: List[int], coordinates: np.ndarray, inputId: int):
        """
        residueIds: list of residue identifiers
        coordinates: numpy array of shape (n_residues, 3)
        groupClass: input identifier (int)
        """
        self.residueIds = residueIds  # e.g., [1, 2, 3, 4]
        self.coords = coordinates  # shape (n_residues, 3)
        self.groupClass = inputId

    def __repr__(self):
        return f"ResidueGroup(class={self.groupClass}, residues={self.residueIds}, size={len(self.residueIds)})"

    def centroid(self) -> np.ndarray:
        """Geometric center of all residues in this group."""
        return np.mean(self.coords, axis=0) if len(self.coords) > 0 else None


def residuesSimilarity(group1: ResidueGroup, group2: ResidueGroup) -> float:
    """Set similarity based on residue IDs."""
    set1 = set(group1.residueIds)
    set2 = set(group2.residueIds)

    intersection = len(set1 & set2)
    smallOne = min(len(set1), len(set2))

    return intersection / smallOne


def spatialSimilarity(group1: ResidueGroup, group2: ResidueGroup, maxDistance: float = 20.0) -> float:
    """
    Spatial similarity based on the groups centroids distance
    """
    rmsd = np.linalg.norm(group1.centroid() - group2.centroid())
    return max(0, 1 - rmsd / maxDistance)


def combinedSimilarity(group1: ResidueGroup, group2: ResidueGroup, spatialWeight: float = 0.6) -> float:
    """
    Combined similarity metric.
    spatialWeight: 0 = only set similarity, 1 = only spatial similarity
    """
    jaccard = residuesSimilarity(group1, group2)
    spatial = spatialSimilarity(group1, group2)

    return (1 - spatialWeight) * jaccard + spatialWeight * spatial


def buildResidueSimilarityGraph(groups: List[ResidueGroup], minSimilarity: float = 0.2, 
                                spatialWeight: float = 0.6) -> nx.Graph:
    """
    Build weighted graph where nodes are residue groups
    and edges represent similarity between groups.
    """
    graphG = nx.Graph()

    # Add nodes with attributes
    for i, group in enumerate(groups):
        graphG.add_node(i,
                   residueIds=group.residueIds,
                   coordinates=group.coords,
                   centroid=group.centroid())

    # Add weighted edges
    n = len(groups)
    for i in range(n):
        for j in range(i + 1, n):
            similarity = combinedSimilarity(groups[i], groups[j], spatialWeight)

            if similarity >= minSimilarity:
                graphG.add_edge(i, j, weight=similarity)

    print(f"Built graph with {graphG.number_of_nodes()} nodes and {graphG.number_of_edges()} edges")
    return graphG


def detectClusters(graphG: nx.Graph, method: str = "louvain", resolution: float = 1.0) -> List[List[int]]:
    """
    Detect communities (Clusters) in the similarity graph.

    Returns:
        List of communities, each community is a list of group indices
    """
    if method == "louvain":
        import community as community_louvain
        partition = community_louvain.best_partition(
            graphG, weight='weight', resolution=resolution
        )

        # Group nodes by community
        communitiesDict = {}
        for node, commId in partition.items():
            communitiesDict.setdefault(commId, []).append(node)

        return list(communitiesDict.values())

    elif method == "connected_components":
        # Fallback: use connected components
        return [list(comp) for comp in nx.connected_components(graphG)]

    else:
        raise ValueError(f"Unsupported method: {method}")


def getClusterRep(communityGroups: List[ResidueGroup], method: str, classWeighted: bool, spatialWeight: float,
                  minFreq: float, minClasses: int, specificClass: int = None,
                  controlSize: bool = True, sizeThres: float = 0.2):
    """
    Extract a spatially dense representative set of residues for a community.

    Args:
        communityGroups: All ResidueGroups in this community
        method: 'centroid', 'intersection', or 'bigger'

    Returns:
        Set of residue IDs forming the dense core
    """

    if method == 'centroid':
        repResiduesList = getRepresentativeCentroid(communityGroups, classWeighted, spatialWeight, specificClass)
    elif method == 'intersection':
        repResiduesList = getRepresentativeIntersection(communityGroups, minFreq, minClasses, controlSize, sizeThres)
    elif method == 'contained':
        repResiduesList = getRepresentativeContained(communityGroups, specificClass)
    else:
        repResiduesList = getRepresentativeBigger(communityGroups, specificClass)

    return repResiduesList

def checkClassFilter(commGroups, minClasses, sameClass):
    if sameClass:
        doPass = len(commGroups) >= minClasses
    else:
        classDistribution = Counter(g.groupClass for g in commGroups)
        uniqueClasses = set(classDistribution.keys())
        doPass = len(uniqueClasses) >= minClasses
    return doPass
    

def residueCommunityPipeline(
        groups: List[ResidueGroup],
        minSimilarity: float = 0.2,
        spatialWeight: float = 0.6,
        communityMethod: str = "louvain",
        resolution: float = 1.0,
        minClasses: int = 2,
        sameClass: bool = False
):
    """
    Builds the clusters of pockets using the community detection in graph.

    Args:
        minClasses: Minimum number of distinct classes required in a cluster

    Returns:
        clusterGroups: clusters of ResidueGroup objects
        len(communities): number of original clusters prior to class filtering (for summary only)
    """
    graphG = buildResidueSimilarityGraph(groups, minSimilarity, spatialWeight)
    communities = detectClusters(graphG, communityMethod, resolution)

    clusterGroups = []
    for commId, commIndices in enumerate(communities):
        # Get groups in this community
        commGroups = [groups[idx] for idx in commIndices]
        if checkClassFilter(commGroups, minClasses, sameClass):
            clusterGroups.append(commGroups)
    
    return clusterGroups, len(communities)


def analyzeResults(results: List[Dict], nComms):
    """Print detailed analysis including class diversity."""

    s = "\n=== Consensus pockets analysis ==="

    if results:
        s += f"\nTotal clusters found: {nComms}"
        s += f"\nClusters passing class filter: {len(results)} ({len(results) / nComms:.1%})"

    s += f"\n\nFiltered clusters: {len(results)}"

    if not results:
        s += "\nNo clusters passed the class filter."
        return s

    # Summary statistics
    avgClasses = np.mean([r['n_classes'] for r in results])
    avgGroups = np.mean([r['n_groups'] for r in results])
    avgRepresentative = np.mean([r['n_representative'] for r in results])

    s += f"\nAverage input classes per clusters: {avgClasses:.1f}"
    s += f"\nAverage pockets per clusters: {avgGroups:.1f}"
    s += f"\nAverage representative size: {avgRepresentative:.1f} residues"

    # Detailed community info
    s += "\nDetailed clusters information:"
    for res in sorted(results, key=lambda x: (x['n_classes'], x['n_groups']), reverse=True):
        s += f"\n\nCluster {res['community_id']}:"
        s += f"\n  Contains {res['n_groups']} pockets from {res['n_classes']} input classes"
        s += f"\n  Classes: {sorted(res['uniqueClasses'])}"

        # Show class distribution
        s += "\n  Class distribution:"
        for cls, count in sorted(res['classDistribution'].items()):
            fraction = count / res['n_groups']
            s += f"\n    {cls}: {count} pockets ({fraction:.1%})"

        s += f"\n  Representative: {res['n_representative']} residues"
        s += f"\n  Coverage: {res['coverage']:.1%}"

        # Show first few representative residues
        repList = sorted(res['representative_residues'])
        if len(repList) <= 8:
            s += f"\n  Representative residues: {repList}"
        else:
            s += f"\n  Representative residues (first 8): {repList[:8]}..."
    return s


def getRepresentativeBigger(
        communityGroups: List[ResidueGroup],
        specificClass: int = None
):
    """
    Get representative as the group with the most residues.

    Args:
        tie_breaker: How to handle ties ("random", "first", "centroid")
    """
    if specificClass is not None:
        communityGroups = [g for g in communityGroups if g.groupClass == specificClass]

    if not communityGroups:
        return None, set()

    # Find groups with maximum size
    maxSize = max(len(g.residueIds) for g in communityGroups)
    candidateIndices = [
        i for i, g in enumerate(communityGroups)
        if len(g.residueIds) == maxSize
    ]

    # Handle ties
    if len(candidateIndices) == 1:
        selectedIdx = candidateIndices[0]
    else:
        # Among ties, choose the centroid
        avgSimilarities = []
        for idx in candidateIndices:
            # Calculate average similarity to all groups
            totalSim = 0
            for j, otherGroup in enumerate(communityGroups):
                if idx != j:
                    totalSim += combinedSimilarity(
                        communityGroups[idx], otherGroup
                    )
            avgSimilarities.append(totalSim / (len(communityGroups) - 1))

        selectedIdx = candidateIndices[np.argmax(avgSimilarities)]

    return [(selectedIdx, set(communityGroups[selectedIdx].residueIds))]

def getRepresentativeContained(
        communityGroups: List[ResidueGroup],
        specificClass: int = None
):
    """
    Get representative as the group that is contained in some other group

    """

    if not communityGroups:
        return None, set()

    containedGroups = []
    for i, g in enumerate(communityGroups):
        gRes = set(g.residueIds)
        for g2 in communityGroups:
            if g != g2:
                g2Res = set(g2.residueIds)
                if len(gRes - g2Res) == 0 and (not specificClass or g.groupClass == specificClass):
                    containedGroups.append((i, gRes))
                    break

    return containedGroups

def getRepresentativeIntersection(
        communityGroups: List[ResidueGroup],
        minFreq: float = 0.7,
        minClasses: int = 2,
        controlSize: bool = True,
        sizeThres: float = 0.2):
    """
    Intersection representative with class awareness.

    Args:
        minFreq: minimum frequency in a cluster for a residue to be considered intersection
        controlSize: the output representative size is controled to be as much as the average of the community groups
        sizeThres: if not None, the output group size should be close to the average of the community groups size.
                      This threshold establishes how flexible this limit is.
                      0.2 means it cannot be bigger than 20% of average
        minClasses: minimum number of classes a residue needs to appear in to be considered
    """
    if not communityGroups:
        return None, set()

    # Count occurrences and track classes
    residueInfo = {}  # residueId -> {'count': 0, 'classes': set()}
    groupSizes = []
    for group in communityGroups:
        groupClass = group.groupClass
        groupSizes.append(len(group.residueIds))
        for residueId in group.residueIds:
            if residueId not in residueInfo:
                residueInfo[residueId] = {'count': 0, 'classes': set()}
            residueInfo[residueId]['count'] += 1
            residueInfo[residueId]['classes'].add(groupClass)

    totalGroups = len(communityGroups)
    threshold = minFreq * totalGroups
    communitySize = sum(groupSizes) / len(groupSizes)
    sizeThreshold = communitySize * (1 + sizeThres)

    sortedItems = sorted(residueInfo.items(), key=lambda item: item[1]['count'], reverse=True)

    # Filter residues
    prevFreq = totalGroups
    representative = set()
    for residueId, info in sortedItems:
        if info['count'] >= threshold:
            toAdd = True
            if controlSize and not ((info['count'] < prevFreq and len(representative) < communitySize)
                                or (info['count'] == prevFreq and len(representative) < sizeThreshold)):
                toAdd = False

            if toAdd and len(info['classes']) >= minClasses:
                representative.add(residueId)
                prevFreq = info['count']

    return [(None, representative)]


def getRepresentativeCentroid(
        communityGroups: List[ResidueGroup],
        classWeighted: bool = True,
        spatialWeight: float = 0.6,
        specificClass: int = None):
    """
    Centroid representative with class weighting.
    Args:
        classWeighted: weight the importance of the ROIs based on how numerous they are (the less the better)
        spatialWeight: weight for the spatial similarity between ROIs (jaccard similarity weights the opposite proportion)
    """

    if not communityGroups or (specificClass is not None and specificClass not in [g.groupClass for g in communityGroups]):
        return None, set()

    if len(communityGroups) == 1:
        return 0, set(communityGroups[0].residueIds)

    n = len(communityGroups)
    similarities = np.zeros((n, n))
    classWeights = np.ones(n)

    # Calculate class diversity weights if requested
    if classWeighted:
        classCounts = Counter(g.groupClass for g in communityGroups)
        for i, group in enumerate(communityGroups):
            # Groups from rarer classes get higher weight
            classWeights[i] = 1.0 / classCounts[group.groupClass]
        # Normalize weights
        classWeights = classWeights / classWeights.sum()

    # Calculate weighted similarities
    for i in range(n):
        for j in range(i + 1, n):
            sim = combinedSimilarity(communityGroups[i], communityGroups[j], spatialWeight)
            # Weight by class importance
            weightedSim = sim * (classWeights[i] + classWeights[j]) / 2
            similarities[i, j] = weightedSim
            similarities[j, i] = weightedSim

    # Find medoid
    totalSimilarities = similarities.sum(axis=1)
    for i in sorted(range(len(totalSimilarities)), key=lambda i: totalSimilarities[i], reverse=True):
        if specificClass is None or communityGroups[i].groupClass == specificClass:
            medoidIdx = i

    return [(medoidIdx, set(communityGroups[medoidIdx].residueIds))]


CENTROID, INTERSEC, BIGGEST, CONTAIN = 0, 1, 2, 3
LOUV, CONNECT = 0, 1

def mapMsaResidues(alFile):
    '''Return the mapping of the residues resNumber -> AlignPos
    {seqIdx: {originalPos: alignPos}}
    '''
    alignment = AlignIO.read(alFile, "clustal")
    align2Ori = {}
    for seqIdx, record in enumerate(alignment):
        oriPos = 0
        align2Ori[seqIdx] = {}
        for i, residue in enumerate(record.seq):
            if residue != "-":
                align2Ori[seqIdx][oriPos] = i
                oriPos += 1
    return align2Ori


class ProtocolConsensusStructROIs(EMProtocol):
    """
    AI Generated:

        This protocol performs consensus analysis of structural regions of interest
        (structural ROIs or pockets) across multiple input datasets, potentially
        derived from different protein structures or prediction methods.

        It identifies spatially and sequence-consistent binding pockets by building
        a residue-based similarity graph and applying community detection algorithms
        to extract consensus pocket clusters.

        The protocol supports multiple strategies for defining representative
        pockets, including centroid-based selection, intersection-based consensus
        residue sets, and structure-size-based selection.

        Core Concepts
        -------------
        Structural ROI (Pocket):
            A set of protein residues forming a spatial region of interest,
            typically representing a binding pocket or functional site.

        Residue Group:
            An abstraction of a pocket represented by:
            - A list of residue identifiers
            - 3D coordinates of residues
            - Input dataset identifier (class)

        Consensus Pocket:
            A cluster of residue groups that represent a structurally and/or
            sequence-consistent binding region across multiple inputs.

        Residue Similarity:
            Quantifies overlap between residue sets (Jaccard-like overlap ratio).

        Spatial Similarity:
            Measures proximity between pocket centroids in 3D space.

        Combined Similarity:
            Weighted combination of residue overlap and spatial proximity.

        Workflow
        --------
        1. Input multiple SetOfStructROIs (pocket sets).
        2. Optionally map residues and chains across different proteins using
           sequence alignment (MSA + chain matching).
        3. Extract residue sets and 3D coordinates for each pocket.
        4. Build a similarity graph where nodes are pockets and edges represent
           combined residue + spatial similarity.
        5. Detect communities using:
           - Louvain algorithm (modularity optimization), or
           - Connected components (threshold-based clustering)
        6. Filter clusters based on:
           - Minimum number of input sources contributing
           - Optional class diversity constraints
        7. Generate consensus representatives:
           - Existing representative pocket from cluster, or
           - Synthetic pocket built from residue intersection
        8. Optionally generate independent consensus sets per input dataset.
        9. Output final consensus structural ROI sets and summary statistics.

        Clustering Criteria
        -------------------
        - minSimil:
            Minimum combined similarity required to connect two pockets in graph.

        - spatialW:
            Weight of spatial proximity in similarity calculation
            (0 = sequence overlap only, 1 = spatial only).

        - clustMethod:
            Community detection method:
            - Louvain (resolution-controlled modularity optimization)
            - Connected components (hard threshold clustering)

        - resolution:
            Controls granularity of Louvain clustering.

        Cluster Filtering
        -----------------
        - numOfOverlap:
            Minimum number of contributing input sets required per cluster.

        - sameClust:
            Whether multiple ROIs from the same input count separately.

        - minClasses:
            Minimum number of distinct input classes required in a cluster.

        Representative Selection Modes
        ------------------------------
        Centroid:
            Selects the pocket most similar to all others in the cluster,
            optionally weighted by class frequency and spatial similarity.

        Intersection:
            Constructs a new pocket composed of residues present in a minimum
            fraction of cluster members.

        Bigger:
            Selects the pocket with the largest number of residues.

        Contained:
            Selects pockets that are fully contained within others in the cluster.

        Internal Processing Logic
        -------------------------
        - Protein chains are optionally aligned across datasets using MSA.
        - Residue indices are mapped to a reference protein when necessary.
        - Pocket residues are converted into residue groups with coordinates.
        - A similarity graph is constructed using pairwise group comparison.
        - Community detection identifies structurally consistent regions.
        - Residue mapping ensures cross-protein comparability.
        - Representative pockets are selected or constructed per cluster.
        - Optional independent consensus outputs are generated per input set.
        - Summary statistics include class distribution and coverage metrics.

        Output
        ------
        - outputStructROIs:
            Set of consensus structural ROI objects representing clustered pockets.

        - outputStructROIs_<XXX> (optional):
            Independent consensus sets per input dataset.

        - summary.txt:
            Detailed report of cluster composition, class distribution,
            and representative residue coverage.

        Use Cases
        ---------
        - Identification of conserved binding pockets across protein variants
        - Comparison of predicted pockets from multiple methods
        - Detection of functionally relevant conserved regions
        - Reduction of redundant or overlapping pocket predictions
        - Cross-structure consensus analysis of binding sites
        """
    _label = 'Consensus structural ROIs'
    _possibleOutputs = PredictStructROIsOutput
    clustChoices = ['Louvain', 'Connected components']
    repChoices = ['Centroid', 'Intersection', 'Bigger', 'Contained']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        g1 = form.addGroup('Input')
        g1.addParam('inputStructROIsSets', params.MultiPointerParam,
                    pointerClass='SetOfStructROIs', allowsNull=False,
                    label="Input Sets of Structural ROIs: ",
                    help='Select the structural ROIs sets to make the consensus')

        g2 = form.addGroup('Clustering')
        g2.addParam('clustMethod', params.EnumParam, default=0,
                    label='Clustering method: ', choices=self.clustChoices,
                    help='Clustering method based on graph communities based on residues overlapping.'
                         '\nLouvain selects the best partition of communities.'
                         '\nConnected components select those clusters with enough similarity.')
        g2.addParam('minSimil', params.FloatParam, default=0.8, label='Minimum similarity: ',
                    help="Minimum similarity for two ROIs to be considered in the same cluster. "
                         "It is calculated as the weighted sum of residue proportion of overlap and spatial similarity")
        g2.addParam('spatialW', params.FloatParam, default=0.4, label='Spatial weight: ',
                    expertLevel=params.LEVEL_ADVANCED,
                    help="Weight of the spatial part of the similarity. It takes into account how close the residues "
                         "in two ROIs are, disregarding the overlap.")

        g2.addParam('resolution', params.FloatParam, default=0.9, label='Resolution for Louvain algorithm: ',
                    expertLevel=params.LEVEL_ADVANCED, condition='clustMethod==0',
                    help='Resolution used in the Louvain algorithm')


        g3 = form.addGroup('Representative')
        g3.addParam('repChoice', params.EnumParam, default=INTERSEC,
                    label='Representant choice: ', choices=self.repChoices,
                    help='How to choose the representative ROI from a cluster of overlapping ROIs. \n'
                         'Centroid: chooses the ROI with bigger similarity to the rest in the cluster.\n'
                         'Intersection: creates a new standard pocket with the residues shared by x proportion of '
                         'the cluster ROIS\nBigger: chooses the ROI with higher number of residues in contact.')
        g3.addParam('numOfOverlap', params.IntParam, default=2,
                    label='Minimun number of overlapping structural regions: ',
                    help="Min number of structural regions to be considered consensus StructROIs")
        g3.addParam('sameClust', params.BooleanParam, default=False,
                    label='Count ROIs from same input: ', expertLevel=params.LEVEL_ADVANCED,
                    help='Whether to count overlapping structural ROIs from the same input set when calculating the '
                         'cluster size')

        g3.addParam('minFreq', params.FloatParam, default=0.6, label='Minimum frequency for intersection: ',
                    expertLevel=params.LEVEL_ADVANCED, condition=f'repChoice=={INTERSEC}',
                    help='Minimum frequency a residue bust appear in a cluster to be considered intersection.')
        g3.addParam('controlSize', params.BooleanParam, default=True, label='Control intersection size: ',
                    expertLevel=params.LEVEL_ADVANCED, condition=f'repChoice=={INTERSEC}',
                    help='Control intersection size so it does not grow far from the cluster average.')
        g3.addParam('sizeThres', params.FloatParam, default=0.2, label='Control size threshold: ',
                    expertLevel=params.LEVEL_ADVANCED, condition=f'repChoice=={INTERSEC} and controlSize',
                    help='Hard threshold to control the size so it does not grow x over the cluster average.')
    
        g3.addParam('outIndv', params.BooleanParam, default=False, condition=f'repChoice not in [{INTERSEC}]',
                    label='Output for each input: ', expertLevel=params.LEVEL_ADVANCED,
                    help='Creates an output set related to each input set, with the elements from each input'
                         'present in the consensus clusters')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep(self.consensusStep)
        self._insertFunctionStep(self.createOutputStep)

    def consensusStep(self):
        self.doMap = False
        if self.checkDifferentInput():
            self.chainsMapDic = self.buildChainsMapDic()
            self.residuesMapDic = self.buildResiduesMapDic()
            self.resId2PositionDic, self.resPosition2IdDic = self.buildResPositionMapDic()
            self.doMap = True

        self.pocketDic = self.buildPocketDic()

        residueGroups, groupDic = [], {}
        for i, pockSet in enumerate(self.inputStructROIsSets):
            for newPock in pockSet.get():
                coords = self.getPocketCoords(newPock)
                residues = self.getPocketResidues(newPock)
                inSetId = self.pocketDic[newPock.getFileName()]
                group = ResidueGroup(residues, coords, inSetId)

                residueGroups.append(group)
                groupDic[group] = newPock.clone()

        pocketClusters, nComms = self.generatePocketClusters(residueGroups)
        self.consensusPockets, clustersConsensus = self.buildStructROIs(pocketClusters, groupDic)
        self.buildSummary(clustersConsensus, self.consensusPockets, nComms)

        if self.outIndv.get() and self.repChoice.get() != INTERSEC:
            self.indepConsensusSets = {}
            for inSetId in range(len(self.inputStructROIsSets)):
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId], _ = self.buildStructROIs(pocketClusters, groupDic, inSetId)

    def createOutputStep(self):
        self.consensusPockets = self.fillEmptyAttributes(self.consensusPockets)
        self.consensusPockets, idsDic = self.reorderIds(self.consensusPockets)

        outPockets = SetOfStructROIs(filename=self._getPath('ConsensusStructROIs_All.sqlite'))
        for i, outPock in enumerate(self.consensusPockets):
            newPock = outPock.clone()
            newPock.setObjId(i+1)
            if newPock.getVolume() is None:
                newPock.setVolume(newPock.getPocketVolume())
            outPockets.append(newPock)

        if outPockets.getSize() > 0:
            outPockets.buildPDBhetatmFile(suffix='_All')
            self._defineOutputs(**{self._possibleOutputs.outputStructROIs.name: outPockets})

        if self.outIndv.get():
            indepOutputs = self.createIndepOutputs()
            for setId in indepOutputs:
                #Index should be the same as in the input
                suffix = '_{:03d}'.format(setId+1)
                outName = 'outputStructROIs' + suffix
                outSet = indepOutputs[setId]
                if outSet.getSize() > 0:
                    outSet.buildPDBhetatmFile(suffix=suffix)
                    self._defineOutputs(**{outName: outSet})
                    self._defineSourceRelation(self.inputStructROIsSets[setId].get(), outSet)

    # --------------------------- INFO functions -----------------------------------

    def buildSummary(self, pocketClusters, reps, nComms):
        results = []
        for i, cluster in enumerate(pocketClusters):
            # Calculate community statistics
            allResidues = set()
            allCoords = []
            for group in cluster:
                allResidues.update(group.residueIds)
                allCoords.append(group.coords)

            classDistribution = Counter(g.groupClass for g in cluster)
            uniqueClasses = set(classDistribution.keys())
            repResidues = reps[i].getDecodedCResidues()

            # Store results including class info
            result = {
                'community_id': i, 'groups': cluster, 'n_groups': len(cluster),
                'classDistribution': classDistribution,
                'uniqueClasses': uniqueClasses, 'n_classes': len(uniqueClasses),
                'allResidues': allResidues, 'n_all_residues': len(allResidues),
                'representative_id': reps[i], 'representative_residues': repResidues,
                'n_representative': len(repResidues),
                'coverage': len(repResidues) / len(allResidues) if allResidues else 0,
            }
            results.append(result)

        with open(self.getSummaryFile(), 'w') as f:
            f.write(analyzeResults(results, nComms))

    def getSummaryFile(self):
        return self._getExtraPath('summary.txt')

    def _summary(self):
        summary = []
        if os.path.exists(self.getSummaryFile()):
            with open(self.getSummaryFile()) as f:
                summary += [f.read()]

        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        """ Try to find warnings on define params. """
        warnings = []
        if self.checkDifferentInput():
            names = self.getInputProtNames()
            warnings.append('The protein this structural ROIs are calculated might not be the same for all the inputs.'
                  '\nDetected protein names: {}'.format(' | '.join(set(names))))
        return warnings

    def checkDifferentInput(self):
        diff = False
        names = self.getInputProtNames()
        if len(set(names)) > 1:
            diff = True
        return diff

    def getInputProtNames(self):
        names = []
        for inSet in self.inputStructROIsSets:
            names.append(inSet.get().getProteinName())
        return names


    # --------------------------- UTILS functions -----------------------------------
    def buildStructROIs(self, clustersGroups, groupDic, specificClass=None):
        outPockets, newGroups = [], []
        asFile = self.inputStructROIsSets[0].get().getProteinFile()
        repMethod = self.getEnumText('repChoice').lower()

        i = 0
        for cluster in clustersGroups:
            repResiduesList = getClusterRep(cluster, repMethod, not self.sameClust.get(), self.spatialW.get(),
                                               self.minFreq.get(), self.numOfOverlap.get(), specificClass,
                                               self.controlSize.get(), self.sizeThres.get())

            for repId, repResidues in repResiduesList:
                i += 1
                if repId is not None:
                    outPocket = groupDic[cluster[repId]]
                    outPockets.append(outPocket)
                    if cluster not in newGroups:
                        newGroups.append(cluster)
                else:
                    if len(repResidues) > 0:
                        pocketFile = self._getExtraPath(f'pocketFile_{i}.cif')
                        coords = self.getResidueCoords(repResidues, asFile, avgResidue=False)
                        createPocketFile(coords, i, pocketFile)
                        outPocket = StructROI(pocketFile, asFile)
                        outPocket.setContactResidues(outPocket.encodeIds(repResidues))
                        outPockets.append(outPocket)
                        if cluster not in newGroups:
                            newGroups.append(cluster)
        return outPockets, newGroups

    def getPrefSetId(self):
        """Return the index of the set whose name matches fromSet"""
        targetName = self.fromSet.get()
        if targetName:
            for idx, setPointer in enumerate(self.inputStructROIsSets):
                roiSet = setPointer.get()
                if str(roiSet) == targetName:
                    return idx
        return None

    def isPocketInside(self, small, big):
        resSmall = set(small.getDecodedCResidues())
        resBig = set(big.getDecodedCResidues())
        return resSmall.issubset(resBig)

    def getInputChainSequences(self):
        '''Returns a dict: {inpIndex: {chainId: protSeq, ...}, ...}
        '''
        dic = {}
        for i, pPointer in enumerate(self.inputStructROIsSets):
            pSet = pPointer.get()
            dic[i] = pSet.getProteinSequencesDic()
        return dic

    def getSequencesSimilarity(self, seq1, seq2):
        blos = substitution_matrices.load('BLOSUM80')
        align = pairwise2.align.globaldx(seq1, seq2, blos)[0]
        return align.score

    def getPairWiseChainMapping(self, seqDic1, seqDic2):
        simMatrix = np.zeros((len(seqDic1), len(seqDic2)))
        for i, (_, newSeq) in enumerate(seqDic1.items()):
            for j, (_, refSeq) in enumerate(seqDic2.items()):
                simMatrix[i, j] = self.getSequencesSimilarity(newSeq, refSeq)

        rowInd, colInd = linear_sum_assignment(-simMatrix)
        matches = list(zip(rowInd, colInd))
        refKeys, newKeys = list(seqDic2.keys()), list(seqDic1.keys())
        return {newKeys[i]: refKeys[j] for i, j in matches}

    def buildChainsMapDic(self):
        '''Build a mapping of the chains to the reference protein (first setOfStructROIs protein).
        The chains will be mapped by pairwise alignment of the sequences to the reference protein
        It will be a dictionary where first level is the index input.
        Second level is a mapping where the key is the chainId of that input and the value is the reference chainId.
        Therefore, first item will be just the first input mapped to itself

        :return: {  inpIndex0: {refChainId0: refChainId0, ...},
                    inpIndex1: {chainId11: refChainId0, chainId10: refChainId1, ...},
                  }
        '''
        inpChainDic = self.getInputChainSequences()
        for i in inpChainDic:
            if i == 0:
                refSeqDic = inpChainDic[0]
                mapDic = {0: {k: k for k in refSeqDic}}
            else:
                seqDic = inpChainDic[i]
                mapDic[i] = self.getPairWiseChainMapping(seqDic, refSeqDic)

        return mapDic

    def getChainGroups(self):
        '''Returns a list of lists with the grouped chains:
        [[chain00, chain10, chain20], [chain01, chain11, chain21]]'''
        for idx, protDic in self.chainsMapDic.items():
            if idx == 0:
                chainGroups = [[chain] for chain in protDic]
            else:
                mDic = flipDic(protDic)
                for group in chainGroups:
                    if group[0] in mDic:
                        group.append(mDic[group[0]])
        return chainGroups

    def buildResiduesMapDic(self):
        '''Build a mapping of the residues to the reference protein (first setOfStructROIs protein).
        The residues of each chain will be mapped to the ones in the reference protein by a multiple alignment

        It will be a dictionary where first level is the index input.
        Second level is a mapping where the key is the chainId of that input and the value is the reference chainId.
        Therefore, first item will be just the first input mapped to itself

        :return: {  refChain0: {seqIdx: {originalPos: alignPos}, seqIdx2: {originalPos: alignPos}, ...}
                    refChain1: {seqIdx: {originalPos: alignPos}, seqIdx2: {originalPos: alignPos}, ...},
                  }
        '''
        chainGroups = self.getChainGroups()
        inpChainDic = self.getInputChainSequences()

        resMap = {}
        for group in chainGroups:
            refChain = group[0]
            alignFile = self.performMSA(inpChainDic, group, refChain)
            resMap[refChain] = mapMsaResidues(alignFile)
        return resMap

    def buildResPositionMapDic(self):
        '''Builds a dictionary that maps the ID of each residue in the inputs to their position index in the sequence.
        {seqIdx: {chainID: {resId: resPos}}}
        Also returns the inverse dic:
        {seqIdx: {chainID: {resPos: resId}}}
        '''
        dic, invDic = {}, {}
        for i, pPointer in enumerate(self.inputStructROIsSets):
            pSet = pPointer.get()
            dic[i] = pSet.getProteinSequencesResIdsDic()
            invDic[i] = {}
            for chain in dic[i]:
                invDic[i][chain] = invertDic(dic[i][chain])
        return dic, invDic

    def performMSA(self, inpChainDic, group, refChain):
        iFile, oFile = self._getTmpPath(f'chains_{refChain}.fa'), \
                       self._getExtraPath(f'chains_{refChain}_alignment.fa')
        with open(iFile, 'w') as f:
            for i, chain in enumerate(group):
                seq = inpChainDic[i][chain]
                f.write(f'>prot{i}_chain{chain}\n{seq}\n')

        cline = getMultipleAlignmentCline('MAFFT', iFile, oFile)
        self.runJob(cline, '')
        return oFile

    def buildPocketDic(self):
        dic = {}
        for i, pPointer in enumerate(self.inputStructROIsSets):
            for pocket in pPointer.get():
                dic[pocket.getFileName()] = i
        return dic

    def createIndepOutputs(self):
        outSets = {}
        for setId, pockSet in self.indepConsensusSets.items():
            suffix = '_{:03d}'.format(setId+1)
            newSet = SetOfStructROIs(filename=self._getExtraPath('ConsensusStructROIs{}.sqlite'.format(suffix)))
            if pockSet:
                for pock in pockSet:
                    newSet.append(pock.clone())
            outSets[setId] = newSet
        return outSets

    def getInpProteinFiles(self):
        inpProteinFiles = []
        for inpSet in self.inputStructROIsSets:
            inpProteinFiles.append(inpSet.get().getProteinFile())
        return inpProteinFiles

    def getPDBName(self):
        pdbFile = self.inputStructROIsSets[0].get().getFirstItem().getProteinFile().split('/')[-1]
        return pdbFile.split('_out')[0]

    def getPocketResidues(self, pock):
        res = pock.getDecodedCResidues()
        if self.doMap:
            inSetId = self.pocketDic[pock.getFileName()]
            res = self.performMapping(res, inSetId)
        return res

    def getResidueCoords(self, residues, asFile, avgResidue=True):
        coordDic = self.parseResidueCoords(asFile)
        atomCoordsList = [coordDic[resId] for resId in residues]
        if avgResidue:
            resCoords = [np.mean(atomCoords, axis=0) for atomCoords in atomCoordsList]
        else:
            resCoords = [atomCoord for atomCoords in atomCoordsList for atomCoord in atomCoords]
        return resCoords

    def getPocketCoords(self, pock):
        asFile = pock.getProteinFile()
        residues = pock.getDecodedCResidues()
        return self.getResidueCoords(residues, asFile)

    def generatePocketClusters(self, residueGroups):
        clustMethod = 'louvain' if self.clustMethod.get() == LOUV else 'connected_components'

        clusterGroups, nCommunities = residueCommunityPipeline(
            residueGroups,
            minSimilarity=self.minSimil.get(),
            spatialWeight=self.spatialW.get(),
            communityMethod=clustMethod,
            resolution=self.resolution.get(),
            minClasses=self.numOfOverlap.get(),
            sameClass=self.sameClust.get(),
        )

        return clusterGroups, nCommunities


    def getIndepClusters(self, clusters):
        indepClustersDic = {}
        for clust in clusters:
            if len(clust) >= self.numOfOverlap.get():
                curIndepCluster = {}
                for pock in clust:
                    curPockFile = pock.getFileName()
                    inSetId = self.pocketDic[curPockFile]
                    if inSetId in curIndepCluster:
                        curIndepCluster[inSetId] += [pock]
                    else:
                        curIndepCluster[inSetId] = [pock]

                for inSetId in curIndepCluster:
                    if inSetId in indepClustersDic:
                        indepClustersDic[inSetId] += [curIndepCluster[inSetId]]
                    else:
                        indepClustersDic[inSetId] = [curIndepCluster[inSetId]]
        return indepClustersDic

    def countPocketsInCluster(self, cluster):
        setIds = []
        for pock in cluster:
            setId = self.pocketDic[pock.getFileName()]
            if self.sameClust.get() or not setId in setIds:
                setIds.append(setId)
        return len(setIds)

    def cluster2representative(self, clusters, minSize=None, onlyRef=None):
        if minSize == None:
            minSize = self.numOfOverlap.get()

        representatives = []
        for i, clust in enumerate(clusters):
            # Check if cluster is big enough
            if self.countPocketsInCluster(clust) >= minSize:
                # Representative only from reference protein if they are different or preferent set
                filteredClust = clust
                if onlyRef is not None:
                    filteredClust = self.filterPocketsBySet(clust, onlyRef)
                elif self.getPrefSetId() is not None:
                    # If there is a Preference set, filter for it
                    filteredClust = self.filterPocketsBySet(clust, self.getPrefSetId())

                outPockets = None
                if self.repChoice.get() == MAXVOL and filteredClust:
                    outPockets = [self.getMaxVolumePocket(filteredClust)]
                elif self.repChoice.get() == MAXSURF and filteredClust:
                    outPockets = [self.getMaxSurfacePocket(filteredClust)]
                elif self.repChoice.get() == CONTAIN and filteredClust:
                    # The contained must be in the filtered set, but can be contained by non filtered ones
                    outPockets = self.getContainedPockets(filteredClust, clust)
                elif self.repChoice.get() == INTERSEC:
                    # As intersection creates a new pocket rather than choosing, it uses all pockets, no the filtered
                    outPockets = [self.getIntersectionPocket(clust, i)]

                if outPockets is not None:
                    representatives += outPockets

        return representatives

    def filterPocketsBySet(self, pockets, setId):
        nPocks = []
        for pock in pockets:
            curPockFile = pock.getFileName()
            inSetId = self.pocketDic[curPockFile]
            if inSetId == setId:
                nPocks.append(pock.clone())
        return nPocks

    def getMaxVolumePocket(self, cluster):
        '''Return the pocket with max volume in a cluster
        The volume is calculated from the convex hull of the contact atoms'''
        maxVol = -1
        for pocket in cluster:
            pocketVol = pocket.getSurfaceConvexVolume()
            if pocketVol > maxVol:
                maxVol = pocketVol
                outPocket = pocket.clone()
        return outPocket

    def getMaxSurfacePocket(self, cluster):
        '''Return the pocket with max surface area in a cluster.
        The surface is calculated from the convex hull of the contact atoms'''
        maxArea, maxVol = 0, 0
        for pocket in cluster:
            pocketArea = pocket.getSurfaceConvexArea()
            if pocketArea >= maxArea:
                maxArea = pocketArea
                outPocket = pocket.clone()
        return outPocket

    def parseResidueCoords(self, asFile):
        resCoordsDic = {}
        parser = PDBParser if asFile.endswith('.pdb') else MMCIFParser

        struct = parser().get_structure(getBaseName(asFile), asFile)[0]
        for chain in struct.get_chains():
            chainId = chain.get_id()
            for res in chain.get_residues():
                resId = res.id[1]
                residueId = f'{chainId}_{resId}'
                resCoordsDic[residueId] = [atom.coord.tolist() for atom in res.get_atoms()]

        return resCoordsDic

    def getContainedPockets(self, filteredCluster, cluster):
        reps = []
        for pock in filteredCluster:
            isInside = False
            for other in cluster:
                if pock.getFileName() != other.getFileName() and self.isPocketInside(pock, other):
                    isInside = True
                    break
            if isInside:
                reps.append(pock.clone())
        return reps

    def getIntersectionPocket(self, cluster, i):
        '''Return the pocket as set of intersection residues in a cluster.'''
        if len(cluster) > 1:
            inters = cluster[0].getDecodedCResidues()
            pdbFile = cluster[0].getProteinFile()
            for pocket in cluster:
                pRes = pocket.getDecodedCResidues()
                inters = set(inters).intersection(set(pRes))

            if len(inters) > 0:
                coords = []
                resCoordDic = self.parseResidueCoords(pdbFile)
                for res in inters:
                    coords += resCoordDic[res]

                pocketFile = self._getExtraPath(f'pocketFile_{i}.cif')
                createPocketFile(coords, i, pocketFile)
                outPocket = StructROI(pocketFile, pdbFile)
                outPocket.setObjId(i)
                outPocket.calculateContacts()
            else:
                outPocket = self.getMaxSurfacePocket(cluster)
        else:
            outPocket = self.getMaxSurfacePocket(cluster)

        return outPocket

    def getPocketsIntersection(self, pock1, pock2):
        res1, res2 = pock1.getDecodedCResidues(), pock2.getDecodedCResidues()
        if self.doMap:
            inSetId1, inSetId2 = self.pocketDic[pock1.getFileName()], self.pocketDic[pock2.getFileName()]
            res1, res2 = self.performMapping(res1, inSetId1), self.performMapping(res2, inSetId2)
        overlap = set(res1).intersection(set(res2))
        return overlap

    def performMapping(self, resList, setId):
        mapResList = []
        for resId in resList:
            chain, res = resId.split('_')
            resIdx = self.resId2PositionDic[setId][chain][int(res)]
            if chain in self.chainsMapDic[setId]:
                refChain = self.chainsMapDic[setId][chain]
                refPos = self.residuesMapDic[refChain][setId][resIdx]
                refRes = self.resPosition2IdDic[setId][refChain][refPos]
            else:
                refChain, refRes = setId, setId
            mapResList.append(f'{refChain}_{refRes}')
        return mapResList

    def calculateResiduesOverlap(self, pock1, pock2):
        res1, res2 = pock1.getDecodedCResidues(), pock2.getDecodedCResidues()
        overlap = self.getPocketsIntersection(pock1, pock2)
        return len(overlap) / min(len(res1), len(res2))

    def getAllPocketAttributes(self, pocketSets):
        '''Return a dic with {attrName: ScipionObj=None}'''
        attributes = {}
        for pockSet in pocketSets:
            item = pockSet.get().getFirstItem()
            attrKeys = item.getObjDict().keys()
            for attrK in attrKeys:
                if not attrK in attributes:
                    value = item.__getattribute__(attrK)
                    attributes[attrK] = value.clone()
                    attributes[attrK].set(None)
        return attributes

    def fillEmptyAttributes(self, inSet):
        '''Fill all items with empty attributes'''
        attributes = self.getAllPocketAttributes(self.inputStructROIsSets)
        for item in inSet:
            for attr in attributes:
                if not hasattr(item, attr):
                    item.__setattr__(attr, attributes[attr])
        return inSet

    def reorderIds(self, inSet):
        '''Return the set with the reordered ids and a mapper dictionary {newId: oldId}'''
        idsDic = {}
        for i, item in enumerate(inSet):
            idsDic[i+1] = item.getObjId()
            item.setObjId(i+1)
        return inSet, idsDic

    def getTemplateOutPDB(self):
        templatePocket = None
        for pock in self.consensusPockets:
            if pock.getPocketClass() != 'AutoLigand':
                templatePocket = pock
                break
        if templatePocket == None:
            templatePocket = self.consensusPockets[0]
        return self.parseATOMlines(templatePocket.getProteinFile())

    def parseHETATM(self, pdbFile, oldId, newId):
        outStr = ''
        with open(pdbFile) as f:
            for line in f:
                if line.startswith('HETATM') and int(splitPDBLine(line)[5]) == int(oldId):
                    outStr += self.pdbLineReplacement(line, str(oldId), str(newId))
        return outStr

    def parseATOMlines(self, pdbFile):
        outStr = ''
        with open(pdbFile) as f:
            for line in f:
                if line.startswith('ATOM'):
                    outStr += line
                elif line.startswith('HETATM') and splitPDBLine(line)[2] != 'APOL':
                    outStr += line
        return outStr

    def pdbLineReplacement(self, line, oldId, newId):
        oldLen, newLen = len(oldId), len(newId)
        oldStr = 'APOL STP C{}'.format((4-oldLen)*' ' + oldId)
        newStr = 'APOL STP C{}'.format((4-newLen)*' ' + newId)
        line = line.replace(oldStr, newStr)
        return line
