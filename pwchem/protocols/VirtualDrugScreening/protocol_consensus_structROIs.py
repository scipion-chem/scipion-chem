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

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwem.objects.data import AtomStruct

from pwchem.objects import SetOfStructROIs, PredictStructROIsOutput
from pwchem.utils import *

from pwchem.constants import *

MAXVOL, MAXSURF = 0, 1

class ProtocolConsensusStructROIs(EMProtocol):
    """
    Executes the consensus on the sets of pockets
    """
    _label = 'Consensus structural ROIs'
    _possibleOutputs = PredictStructROIsOutput
    actionChoices = ['MaxVolume', 'MaxSurface']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputStructROIsSets', params.MultiPointerParam,
                       pointerClass='SetOfStructROIs', allowsNull=False,
                       label="Input Sets of Structural ROIs",
                       help='Select the structural ROIs sets to make the consensus')
        form.addParam('outIndv', params.BooleanParam, default=False,
                      label='Output for each input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Creates an output set related to each input set, with the elements from each input'
                           'present in the consensus clusters')

        form.addParam('overlap', params.FloatParam, default=0.75, label='Proportion of residues for overlapping',
                      help="Min proportion of residues (from the smaller) of two structural regions to be considered "
                           "overlapping")
        form.addParam('action', params.EnumParam, default=MAXSURF,
                      label='Action on overlapping', choices=self.actionChoices,
                      expertLevel=params.LEVEL_ADVANCED,
                      help='Action to take on overlapping structural regions, whether to merge them or keep just one '
                           'with some condition')
        form.addParam('numOfOverlap', params.IntParam, default=2,
                      label='Minimun number of overlapping structural regions',
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
        pocketDic = self.buildPocketDic()
        pocketClusters = self.generatePocketClusters()
        # Getting the representative of the clusters from any of the inputs
        self.consensusPockets = self.cluster2representative(pocketClusters, pocketDic)

        if self.outIndv.get():
            self.indepConsensusSets = {}
            # Separating clusters by input set
            indepClustersDic = self.getIndepClusters(pocketClusters, pocketDic)
            for inSetId in indepClustersDic:
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId] = self.cluster2representative(indepClustersDic[inSetId],
                                                                               pocketDic, minSize=1)


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
        names = []
        for inSet in self.inputStructROIsSets:
            names.append(inSet.get().getProteinName())
        if len(set(names)) > 1:
            warnings.append('The protein this structural ROIs are calculated might not be the same for all the inputs.'
                  '\nDetected protein names: {}'.format(' | '.join(set(names))))
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def buildPocketDic(self):
        dic = {}
        for i, pSet in enumerate(self.inputStructROIsSets):
            for pocket in pSet.get():
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

    def getIndepClusters(self, clusters, pocketDic):
        indepClustersDic = {}
        for clust in clusters:
            if len(clust) >= self.numOfOverlap.get():
                curIndepCluster = {}
                for pock in clust:
                    curPockFile = pock.getFileName()
                    inSetId = pocketDic[curPockFile]
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

    def countPocketsInCluster(self, cluster, pocketDic):
        setIds = []
        for pock in cluster:
            setId = pocketDic[pock.getFileName()]
            if self.sameClust.get() or not setId in setIds:
                setIds.append(setId)
        return len(setIds)

    def cluster2representative(self, clusters, pocketDic, minSize=None):
        if minSize == None:
            minSize = self.numOfOverlap.get()

        representatives = []
        for clust in clusters:
            if self.countPocketsInCluster(clust, pocketDic) >= minSize:
                if self.action.get() == MAXVOL:
                    outPocket = self.getMaxVolumePocket(clust)
                elif self.action.get() == MAXSURF:
                    outPocket = self.getMaxSurfacePocket(clust)
                representatives.append(outPocket)
        return representatives

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

    def calculateResiduesOverlap(self, pock1, pock2):
        res1, res2 = pock1.getDecodedCResidues(), pock2.getDecodedCResidues()
        overlap = set(res1).intersection(set(res2))
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

    def createOutPDB(self, idsDic):
        outStr = self.getTemplateOutPDB()
        for pocket in self.consensusPockets:
            outFile = pocket.getProteinFile()
            newId, oldId = pocket.getObjId(), idsDic[pocket.getObjId()]
            outStr += self.parseHETATM(outFile, oldId, newId)

        outPDBFile = self._getExtraPath(self.getPDBName()) + '_out.pdb'
        with open(outPDBFile, 'w') as f:
            f.write(outStr)
            f.write('\nTER\n')

        pmlFile = self._getExtraPath('{}.pml'.format(self.getPDBName()))
        with open(pmlFile, 'w') as f:
            f.write(PML_STR.format(outPDBFile.split('/')[-1]))
        return os.path.abspath(outPDBFile), os.path.abspath(pmlFile)

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
