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
This protocol is used to make a consensus over several input sets of smallMolecules which have been docked by the
same or different software

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re

class ProtocolConsensusDocking(EMProtocol):
    """
    Executes the consensus on the sets of SmallMolecules. The poses of the different molecules are clustered
    with respect to their RMSD
    """
    _label = 'Consensus docking'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMoleculesSets', params.MultiPointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Sets of Docked Molecules",
                       help='Select the pocket sets to make the consensus')
        form.addParam('maxRMSD', params.FloatParam, default=1, label='Max RMSD for overlap: ',
                      help="Maximum RMSD for clustering different docked molecules")
        form.addParam('numOfOverlap', params.IntParam, default=2, label='Minimum number of overlapping dockings',
                      help="Min number of docked molecules to be considered consensus docking")
        form.addParam('action', params.StringParam, default='',
                      label='Criteria to choose cluster representative: ',
                      help='Criteria to follow on docking clusters to choose a representative')
        form.addParam('maxmin', params.BooleanParam, default=True,
                      label='Keep maximum values: ',
                      help='True to keep the maximum values. False to get the minimum')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('consensusStep')
        self._insertFunctionStep('createOutputStep')

    def consensusStep(self):
        molDic = self.buildMolDic()
        molClusters = self.generateDockingClusters()
        #Getting the representative of the clusters from any of the inputs
        self.consensusMols = self.cluster2representative(molClusters)

        self.indepConsensusSets = {}
        #Separating clusters by input set
        indepClustersDic = self.getIndepClusters(molClusters, molDic)
        for inSetId in indepClustersDic:
            # Getting independent representative for each input set
            self.indepConsensusSets[inSetId] = self.cluster2representative(indepClustersDic[inSetId], minSize=1)


    def createOutputStep(self):
        inputProteinFile = self.inputMoleculesSets[0].get().getProteinFile()

        self.relabelDic = {}
        self.consensusMols = self.fillEmptyAttributes(self.consensusMols)
        self.consensusMols, idsDic = self.reorderIds(self.consensusMols)

        outDocked = SetOfSmallMolecules(filename=self._getPath('consensusDocked_All.sqlite'))
        outDocked.setDocked(True)
        outDocked.setProteinFile(inputProteinFile)
        for outDock in self.consensusMols:
            newDock = outDock.clone()
            newDock = self.relabelPosId(newDock)
            outDocked.append(newDock)
        self._defineOutputs(outputSmallMolecules=outDocked)

        indepOutputs = self.createIndepOutputs()
        for setId in indepOutputs:
            #Index should be the same as in the input
            suffix = '_{:03d}'.format(setId+1)
            outName = 'outputSmallMolecules' + suffix
            outSet = indepOutputs[setId]
            outSet.setDocked(True)
            outSet.setProteinFile(inputProteinFile)

            self._defineOutputs(**{outName: outSet})
            self._defineSourceRelation(self.inputMoleculesSets[setId].get(), outSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        """ Try to find warnings on define params. """
        validations = []
        for pSet in self.inputMoleculesSets:
            pSet = pSet.get()
            if not pSet.isDocked():
                validations.append('Sets of input molecules must be docked first.\n'
                                'Set: {} has not been docked'.format(pSet))

        return validations

    # --------------------------- UTILS functions -----------------------------------
    def buildMolDic(self):
        dic = {}
        for i, pSet in enumerate(self.inputMoleculesSets):
            for mol in pSet.get():
                dic[mol.getPoseFile()] = i
        return dic

    def createIndepOutputs(self):
        outSets = {}
        for setId in self.indepConsensusSets:
            suffix = '_{:03d}'.format(setId+1)
            newSet = SetOfSmallMolecules(filename=self._getExtraPath('consensusSmallMolecules{}.sqlite'.format(suffix)))
            for mol in self.indepConsensusSets[setId]:
                newMol = mol.clone()
                newMol = self.relabelPosId(newMol)
                newSet.append(newMol)
            outSets[setId] = newSet
        return outSets

    def getInpProteinFiles(self):
        inpProteinFiles = []
        for inpSet in self.inputMoleculesSets:
            inpProteinFiles.append(inpSet.get().getProteinFile())
        return inpProteinFiles

    def getPDBName(self):
        pdbFile = self.inputMoleculesSets[0].get().getFirstItem().getProteinFile().split('/')[-1]
        return pdbFile.split('_out')[0]

    def generateDockingClusters(self):
        '''Generate the pocket clusters based on the overlapping residues
        Return (clusters): [[pock1, pock2], [pock3], [pock4, pock5, pock6]]'''
        clusters = []
        #For each set of pockets
        for i, molSet in enumerate(self.inputMoleculesSets):
            #For each of the pockets in the set
            for newMol in molSet.get():
                newClusters, newClust = [], [newMol.clone()]
                #Check for each of the clusters
                for clust in clusters:
                    #Check for each of the pockets in the cluster
                    append2Cluster = True
                    for cMol in clust:
                        if cMol.getMolBase() == newMol.getMolBase():
                            curRMSD = self.calculateMolsRMSD(newMol, cMol)

                            #If there is overlap with the new pocket from the set
                            if curRMSD > self.maxRMSD.get():
                                append2Cluster = False
                                break
                        else:
                            append2Cluster = False

                    #newClust: Init with only the newPocket, grow with each clust which overlaps with newPocket
                    if append2Cluster:
                        newClust += clust
                    #If no overlap, the clust keeps equal
                    else:
                        newClusters.append(clust)
                #Add new cluster containing newPocket + overlapping previous clusters
                newClusters.append(newClust)
                clusters = newClusters.copy()
        return clusters

    def getIndepClusters(self, clusters, molDic):
        indepClustersDic = {}
        for clust in clusters:
            if len(clust) >= self.numOfOverlap.get():
                curIndepCluster = {}
                for mol in clust:
                    curMolFile = mol.getPoseFile()
                    inSetId = molDic[curMolFile]
                    if inSetId in curIndepCluster:
                        curIndepCluster[inSetId] += [mol]
                    else:
                        curIndepCluster[inSetId] = [mol]

                for inSetId in curIndepCluster:
                    if inSetId in indepClustersDic:
                        indepClustersDic[inSetId] += [curIndepCluster[inSetId]]
                    else:
                        indepClustersDic[inSetId] = [curIndepCluster[inSetId]]
        return indepClustersDic

    def cluster2representative(self, clusters, minSize=None):
        if minSize == None:
            minSize = self.numOfOverlap.get()

        representatives = []
        for clust in clusters:
            if len(clust) >= minSize:
                repr = self.getRepresentativeMolecule(clust)
                representatives.append(repr)
        return representatives

    def getRepresentativeMolecule(self, cluster):
        '''Return the docked molecule with max score in a cluster'''
        maxScore = -100000 if self.maxmin.get() else 100000
        for mol in cluster:
            molScore = getattr(mol, self.action.get())
            if (molScore > maxScore and self.maxmin.get()) or \
                    (molScore < maxScore and not self.maxmin.get()):
                maxScore = molScore
                outMol = mol.clone()
        return outMol

    def getMinEnergyMolecule(self, cluster):
        '''Return the docked molecule with min energy in a cluster'''
        minEnergy = 10000
        for mol in cluster:
            molEnergy = mol.getEnergy()
            if molEnergy < minEnergy:
                minEnergy = molEnergy
                outMol = mol.clone()
        return outMol

    def calculateMolsRMSD(self, mol1, mol2):
        posDic1, posDic2 = mol1.getAtomsPosDic(), mol2.getAtomsPosDic()
        if self.checkSameKeys(posDic1, posDic2):
            rmsd=0
            for atomId in posDic1:
                rmsd += calculateDistance(posDic1[atomId], posDic2[atomId])
            return rmsd / len(posDic2)
        else:
            print('Atom ids of the molecules are different and cannot be compared')


    def checkSameKeys(self, d1, d2):
        return set(d1.keys()) == set(d2.keys())

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
        attributes = self.getAllPocketAttributes(self.inputMoleculesSets)
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

    def getOriginalReceptorFile(self):
        return self.inputMoleculesSets[0].get().getProteinFile()

    def relabelPosId(self, mol):
        molName = mol.getMolName()
        posFile = mol.getPoseFile()
        #If mol has not been relabel yet
        if not posFile in self.relabelDic:
            #Set the new posId (begining with one, for each molName)
            if molName in self.relabelDic:
                self.relabelDic[molName] += 1
            else:
                self.relabelDic[molName] = 1

            newPosId = self.relabelDic[molName]
            posExt = os.path.splitext(posFile)[1]

            #Moving new posFile to current protocol with new posId
            newPosFile = re.sub('_\d{}'.format(posExt), '_{}{}'.format(newPosId, posExt), posFile, 1)
            newPosFile = self._getPath(os.path.basename(newPosFile))
            shutil.copy(posFile, newPosFile)
            mol.setPoseFile(newPosFile)

            self.relabelDic[posFile] = newPosFile

        else:
            mol.setPoseFile(self.relabelDic[posFile])

        return mol






