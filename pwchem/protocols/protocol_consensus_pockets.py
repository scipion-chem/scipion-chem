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
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfPockets
from pwem.objects.data import AtomStruct
from ..constants import *
from pwchem.utils import *
import os

MAXVOL, MAXSURF = 0, 1

class ProtocolConsensusPockets(EMProtocol):
    """
    Executes the consensus on the sets of pockets
    """
    _label = 'Consensus pockets'
    actionChoices = ['MaxVolume', 'MaxSurface']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputPocketSets', params.MultiPointerParam,
                       pointerClass='SetOfPockets', allowsNull=False,
                       label="Input Sets of Pockets",
                       help='Select the pocket sets to make the consensus')
        form.addParam('overlap', params.FloatParam, default=0.75, label='Proportion of residues for overlapping',
                      help="Min proportion of residues (from the smaller) of two pockets to be considered overlapping")
        form.addParam('action', params.EnumParam, default=MAXSURF,
                      label='Action on overlapping', choices=self.actionChoices,
                      help='Action to take on overlapping pockets, whther to merge them or keep just one '
                           'with some condition')
        form.addParam('numOfOverlap', params.IntParam, default=2, label='Minimun number of overlapping pockets',
                      help="Min number of pockets to be considered consensus pockets")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('consensusStep')
        self._insertFunctionStep('createOutputStep')

    def consensusStep(self):
        pocketClusters = self.generatePocketClusters()
        self.consensusPockets = self.cluster2pocket(pocketClusters)
        self.indepConsensusSets = self.getIndepConsensus(pocketClusters)

    def createOutputStep(self):
        self.consensusPockets = self.fillEmptyAttributes(self.consensusPockets, self.getAllPocketAttributes())
        self.consensusPockets, idsDic = self.reorderIds(self.consensusPockets)

        outPockets = SetOfPockets(filename=self._getPath('consensusPocketsAll.sqlite'))
        outProtFile, outPmlFile = self.createOutPDB(idsDic)
        for outPock in self.consensusPockets:
            newPock = outPock.clone()
            newPock.setProteinFile(outProtFile)
            newPock.setPmlFile(outPmlFile)
            outPockets.append(newPock)
        self._defineOutputs(outputPocketsAll=outPockets)

        indepOutputs = self.createIndepOutputs()
        for outSet in indepOutputs:
            #Index should be the same as in the input
            i = self.getInpProteinFiles().index(outSet.getProteinFile())
            outName = 'outputPockets{}'.format(i+1)
            self._defineOutputs(**{outName:outSet})
            self._defineSourceRelation(self.inputPocketSets[i].get(), outSet)

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
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def createIndepOutputs(self):
        inputProteinFiles = self.getInpProteinFiles()
        outSets = []
        for proteinFile in self.indepConsensusSets:
            i = inputProteinFiles.index(proteinFile)
            newSet = SetOfPockets(filename=self._getPath('consensusPockets{}.sqlite'.format(i+1)))
            for pock in self.indepConsensusSets[proteinFile]:
                newSet.append(pock.clone())
            outSets.append(newSet)
        return outSets

    def getInpProteinFiles(self):
        inpProteinFiles = []
        for inpSet in self.inputPocketSets:
            inpProteinFiles.append(inpSet.get().getProteinFile())
        return inpProteinFiles

    def getPDBName(self):
        pdbFile = self.inputPocketSets[0].get().getFirstItem().getProteinFile().split('/')[-1]
        return pdbFile.split('_out')[0]

    def generatePocketClusters(self):
        '''Generate the pocket clusters based on the overlapping residues
        Return (clusters): [[pock1, pock2], [pock3], [pock4, pock5, pock6]]'''
        clusters = []
        #For each set of pockets
        for i, pockSet in enumerate(self.inputPocketSets):
            #For each of the pockets in the set
            for newPock in pockSet.get():
                if i == 0:
                    #First, each cluster is formed by just one pocket
                    clusters.append([newPock.clone()])
                else:
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

    def getIndepConsensus(self, clusters):
        outSets = {}
        for clust in clusters:
            if len(clust) >= self.numOfOverlap.get():
                for pock in clust:
                    curProtFile = pock.getProteinFile()
                    if curProtFile in outSets:
                        outSets[curProtFile] += [pock]
                    else:
                        outSets[curProtFile] = [pock]
        return outSets


    def cluster2pocket(self, clusters):
        pockets = []
        for clust in clusters:
            if len(clust) >= self.numOfOverlap.get():
                if self.action.get() == MAXVOL:
                    outPocket = self.getMaxVolumePocket(clust)
                elif self.action.get() == MAXSURF:
                    outPocket = self.getMaxVolumePocket(clust)
                pockets.append(outPocket)
        return pockets

    def getMaxVolumePocket(self, cluster):
        '''Return the pocket with max volume in a cluster'''
        maxVol = 0
        for pocket in cluster:
            pocketVol = pocket.getSurfaceConvexVolume()
            if pocketVol > maxVol:
                maxVol = pocketVol
                outPocket = pocket.clone()
        return outPocket

    def getMaxSurfacePocket(self, cluster):
        '''Return the pocket with max surface in a cluster.
        The surface is just interpolated to the number of contact atoms. In case of even, volume is used'''
        maxSurf, maxVol = 0, 0
        for pocket in cluster:
            pocketSurf = len(pocket.getDecodedCAtoms())
            pocketVol = pocket.getSurfaceConvexVolume()
            if pocketSurf > maxSurf:
                maxSurf, maxVol = pocketSurf, pocketVol
                outPocket = pocket.clone()
            elif pocketSurf == maxSurf:
                if pocketVol > maxVol:
                    maxSurf, maxVol = pocketSurf, pocketVol
                    outPocket = pocket.clone()

        return outPocket

    def calculateResiduesOverlap(self, pock1, pock2):
        res1, res2 = pock1.getDecodedCResidues(), pock2.getDecodedCResidues()
        overlap = set(res1).intersection(set(res2))
        return len(overlap) / min(len(res1), len(res2))

    def getAllPocketAttributes(self):
        attrs = set([])
        for pockSet in self.inputPocketSets:
            attrs = attrs.union(set(pockSet.get().getFirstItem().getObjDict().keys()))
        return attrs

    def fillEmptyAttributes(self, inSet, attributes):
        '''Fill all items with empty attributes'''
        for item in inSet:
            for attr in attributes:
                if not hasattr(item, attr):
                    setAttribute(item, attr, 'None')
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
                if line.startswith('HETATM') and int(line.split()[5]) == int(oldId):
                    outStr += self.pdbLineReplacement(line, str(oldId), str(newId))
        return outStr

    def parseATOMlines(self, pdbFile):
        outStr = ''
        with open(pdbFile) as f:
            for line in f:
                if line.startswith('ATOM'):
                    outStr += line
                elif line.startswith('HETATM') and line.split()[2] != 'APOL':
                    outStr += line
        return outStr

    def pdbLineReplacement(self, line, oldId, newId):
        oldLen, newLen = len(oldId), len(newId)
        oldStr = 'APOL STP C{}'.format((4-oldLen)*' ' + oldId)
        newStr = 'APOL STP C{}'.format((4-newLen)*' ' + newId)
        line = line.replace(oldStr, newStr)
        return line
