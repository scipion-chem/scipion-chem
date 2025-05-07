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

import numpy as np
from Bio import pairwise2, AlignIO
from Bio.Align import substitution_matrices
from scipy.optimize import linear_sum_assignment

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput, StructROI
from pwchem.utils import writePDBLine, splitPDBLine, flipDic
from pwchem.utils.utilsFasta import getMultipleAlignmentCline

MAXVOL, MAXSURF, INTERSEC = 0, 1, 2

def mapMsaResidues(alFile):
    '''Return the mapping of the residues resNumber -> AlignPos
    {seqIdx: {resNum: alignPos}}
    '''
    alignment = AlignIO.read(alFile, "clustal")
    align2Ori = {}
    for seqIdx, record in enumerate(alignment):
        oriPos = 0
        align2Ori[seqIdx] = {}
        for i, residue in enumerate(record.seq):
            if residue != "-":
                oriPos += 1
                align2Ori[seqIdx][oriPos] = i
    return align2Ori


class ProtocolConsensusStructROIs(EMProtocol):
    """
    Executes the consensus on the sets of pockets
    """
    _label = 'Consensus structural ROIs'
    _possibleOutputs = PredictStructROIsOutput
    repChoices = ['MaxVolume', 'MaxSurface', 'Intersection']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputStructROIsSets', params.MultiPointerParam,
                       pointerClass='SetOfStructROIs', allowsNull=False,
                       label="Input Sets of Structural ROIs: ",
                       help='Select the structural ROIs sets to make the consensus')
        form.addParam('outIndv', params.BooleanParam, default=False,
                      label='Output for each input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Creates an output set related to each input set, with the elements from each input'
                           'present in the consensus clusters')

        form.addParam('overlap', params.FloatParam, default=0.75, label='Proportion of residues for overlapping: ',
                      help="Min proportion of residues (from the smaller) of two structural regions to be considered "
                           "overlapping")
        form.addParam('repChoice', params.EnumParam, default=MAXSURF,
                      label='Representant choice: ', choices=self.repChoices,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='How to choose the representative ROI from a cluster of overlapping ROIs. MaxSurface and '
                           'MaxVolume chooses the existing pocket with maximum surface or volume, respectively. '
                           '\nIntersection creates a new standard pocket with the interecting residues from the '
                           'ROIs in the cluster (if any)')
        form.addParam('numOfOverlap', params.IntParam, default=2,
                      label='Minimun number of overlapping structural regions: ',
                      help="Min number of structural regions to be considered consensus StructROIs")
        form.addParam('sameClust', params.BooleanParam, default=False,
                      label='Count ROIs from same input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Whether to count overlapping structural ROIs from the same input set when calculating the '
                           'cluster size')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('consensusStep')
        self._insertFunctionStep('createOutputStep')

    def consensusStep(self):
        self.doMap = False
        if self.checkDifferentInput():
            self.chainsMapDic = self.buildChainsMapDic()
            self.residuesMapDic = self.buildResiduesMapDic()
            self.doMap = True

        self.pocketDic = self.buildPocketDic()
        pocketClusters = self.generatePocketClusters()
        # Getting the representative of the clusters from any of the inputs
        ref = 0 if self.doMap else None
        self.consensusPockets = self.cluster2representative(pocketClusters, onlyRef=ref)

        if self.outIndv.get():
            self.indepConsensusSets = {}
            for inSetId in range(len(self.inputStructROIsSets)):
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId] = self.cluster2representative(pocketClusters, onlyRef=inSetId)

    def createOutputStep(self):
        self.consensusPockets = self.fillEmptyAttributes(self.consensusPockets)
        self.consensusPockets, idsDic = self.reorderIds(self.consensusPockets)

        outPockets = SetOfStructROIs(filename=self._getPath('ConsensusStructROIs_All.sqlite'))
        for outPock in self.consensusPockets:
            newPock = outPock.clone()
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
    def _summary(self):
        summary = []
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
                    group.append(mDic[group[0]])
        return chainGroups

    def buildResiduesMapDic(self):
        '''Build a mapping of the residues to the reference protein (first setOfStructROIs protein).
        The residues of each chain will be mapped to the ones in the reference protein by a multiple alignment

        It will be a dictionary where first level is the index input.
        Second level is a mapping where the key is the chainId of that input and the value is the reference chainId.
        Therefore, first item will be just the first input mapped to itself

        :return: {  refChain0: {seqIdx: {resNum: alignPos}, seqIdx2: {resNum: alignPos}, ...}
                    refChain1: {seqIdx: {resNum: alignPos}, seqIdx2: {resNum: alignPos}, ...},
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
        for setId in self.indepConsensusSets:
            suffix = '_{:03d}'.format(setId+1)
            newSet = SetOfStructROIs(filename=self._getExtraPath('ConsensusStructROIs{}.sqlite'.format(suffix)))
            for pock in self.indepConsensusSets[setId]:
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

    def generatePocketClusters(self):
        '''Generate the pocket clusters based on the overlapping residues
        Return (clusters): [[pock1, pock2], [pock3], [pock4, pock5, pock6]]'''
        clusters = []
        #For each set of pockets
        for i, pockSet in enumerate(self.inputStructROIsSets):
            #For each of the pockets in the set
            for newPock in pockSet.get():
                newClusters, newClust = [], [newPock.clone()]
                #Check for each of the clusters
                for clust in clusters:
                    #Check for each of the pockets in the cluster
                    overClust = False
                    for cPocket in clust:
                        propOverlap = self.calculateResiduesOverlap(newPock, cPocket)
                        #If there is overlap with the new pocket from the set
                        if propOverlap > self.overlap.get():
                            overClust = True
                            break

                    #newClust: Init with only the newPocket, grow with each clust which overlaps with newPocket
                    if overClust:
                        newClust += clust
                    #If no overlap, the clust keeps equal
                    else:
                        newClusters.append(clust)
                #Add new cluster containing newPocket + overlapping previous clusters
                newClusters.append(newClust)
                clusters = newClusters.copy()
        return clusters

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
                # Representative only from reference (first) protein if they are different
                if onlyRef is not None:
                    clust = self.filterPocketsBySet(clust, onlyRef)

                if clust:
                    if self.repChoice.get() == MAXVOL:
                        outPocket = self.getMaxVolumePocket(clust)
                    elif self.repChoice.get() == MAXSURF:
                        outPocket = self.getMaxSurfacePocket(clust)
                    elif self.repChoice.get() == INTERSEC:
                        outPocket = self.getIntersectionPocket(clust, i)

                    representatives.append(outPocket)
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
        maxVol = 0
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

    def createPocketFile(self, coords, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(coords):
                f.write(writePDBLine(['HETATM', str(j), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
        return outFile

    def parsePDBResidueCoords(self, pdbFile):
        resCoordsDic = {}
        with open(pdbFile) as f:
            for line in f:
                if line.startswith('ATOM'):
                    line = splitPDBLine(line)
                    proteinChain, residueNumber = line[4:6]
                    residueId = '{}_{}'.format(proteinChain, residueNumber)
                    coord = list(map(float, line[6:9]))
                    if residueId in resCoordsDic:
                        resCoordsDic[residueId].append(coord)
                    else:
                        resCoordsDic[residueId] = [coord]
        return resCoordsDic


    def getIntersectionPocket(self, cluster, i):
        '''Return the pocket as set of intersection residues in a cluster.'''
        inters = cluster[0].getDecodedCResidues()
        pdbFile = cluster[0].getProteinFile()
        for pocket in cluster:
            pRes = pocket.getDecodedCResidues()
            inters = set(inters).intersection(set(pRes))

        if len(inters) > 0:
            coords = []
            resCoordDic = self.parsePDBResidueCoords(pdbFile)
            for res in inters:
                coords += resCoordDic[res]

            pocketFile = self.createPocketFile(coords, i)
            outPocket = StructROI(pocketFile, pdbFile)
            outPocket.calculateContacts()
        else:
            outPocket = self.getMaxSurfacePocket(cluster)

        return outPocket

    def getPocketsIntersection(self, pock1, pock2):
        res1, res2 = pock1.getDecodedCResidues(), pock2.getDecodedCResidues()
        inSetId1, inSetId2 = self.pocketDic[pock1.getFileName()], self.pocketDic[pock2.getFileName()]

        if self.doMap:
            res1, res2 = self.performMapping(res1, inSetId1), self.performMapping(res2, inSetId2)
        overlap = set(res1).intersection(set(res2))
        return overlap

    def performMapping(self, resList, setId):
        mapResList = []
        for resId in resList:
            chain, res = resId.split('_')
            refChain = self.chainsMapDic[setId][chain]
            refRes = self.residuesMapDic[refChain][setId][int(res)]
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
