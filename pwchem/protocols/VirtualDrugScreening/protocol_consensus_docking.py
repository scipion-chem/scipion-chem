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

import os, re, glob
from scipy.cluster.hierarchy import linkage, fcluster
from sklearn.metrics.pairwise import pairwise_distances

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *


class ProtocolConsensusDocking(EMProtocol):
    """
    Executes the consensus on the sets of SmallMolecules. The poses of the different molecules are clustered
    with respect to their RMSD
    """
    _label = 'Consensus docking'

    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMoleculesSets', params.MultiPointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Sets of Docked Molecules",
                       help='Select the pocket sets to make the consensus')

        group = form.addGroup('Clustering')
        group.addParam('doScipy', params.BooleanParam, default=True,
                      label='Use scipy for clustering: ',
                      help='Whether to use scipy hierarchical clustering methods or an '
                           'aggregative approach (usually slower).')
        group.addParam('linkage', params.EnumParam, default=0,
                       choices=['Single', 'Complete', 'Average', 'Ward', 'Centroid', 'Median'],
                       label='Hyerarchical clustering linkage: ', condition='doScipy',
                       help='Hyerarchical clustering linkage used on scipy clutering.\n'
                            'https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html')
        group.addParam('checkRep', params.BooleanParam, default=True,
                       label='Check only cluster representative: ', condition='not doScipy',
                       help='If True, the RMSD is only calculated on the cluster representative, speeding up '
                            'the process but making it less restrictive')
        form.addParam('outIndv', params.BooleanParam, default=False, condition='doScipy or not checkRep',
                      label='Output for each input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Creates an output set related to each input set, with the elements from each input'
                           'present in the consensus clusters')
        form.addParam('sameClust', params.BooleanParam, default=True,
                      label='Count poses from same input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Whether to count overlapping docked molecules from the same input set when calculating the '
                           'cluster size')

        group.addParam('maxRMSD', params.FloatParam, default=1, label='Max RMSD for overlap: ',
                      help="Maximum RMSD for clustering different docked molecules")
        group.addParam('numOfOverlap', params.IntParam, default=2, label='Minimum number of overlapping dockings: ',
                      help="Min number of docked molecules to be considered consensus docking")

        group = form.addGroup('Representative')
        group.addParam('repAttr', params.StringParam, default='',
                      label='Criteria to choose cluster representative: ',
                      help='Criteria to follow on docking clusters to choose a representative. '
                           'It will extract the representative as the pose with max/min (next argument) value '
                           'of this attribute')
        group.addParam('maxmin', params.BooleanParam, default=True,
                      label='Keep maximum values: ',
                      help='True to keep the maximum values. False to get the minimum')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        cStep = self._insertFunctionStep('convertInputStep', prerequisites=[])
        conStep = self._insertFunctionStep('consensusStep', prerequisites=[cStep])
        self._insertFunctionStep('createOutputStep', prerequisites=[conStep])

    def convertInputStep(self):
        allMols = self.getAllInputMols()
        outDir = self.getInputMolsDir()
        os.mkdir(outDir)

        maeMols, _ = getMAEMoleculeFiles(allMols)
        if len(maeMols) > 0:
            try:
                from pwchemSchrodinger.utils.utils import convertMAEMolSet
                convertMAEMolSet(maeMols, outDir, self.numberOfThreads.get(), updateSet=False, subset=False)
            except ImportError:
                print('Conversion of MAE input files could not be performed because schrodinger plugin is not installed')

    def consensusStep(self):
        molDic = self.buildMolDic()
        minSize = self.numOfOverlap.get()

        doInd = True
        if self.doScipy:
            molClusters = self.generateDockingClustersScipy()
            self.consensusMols = self.cluster2representative(molClusters, molDic, minSize)
        else:
            if not self.checkRep:
                molClusters = self.generateDockingClusters()
                self.consensusMols = self.cluster2representative(molClusters, molDic, minSize)
            else:
                self.consensusMols = self.buildRepresentativeClusters(molDic)
                doInd = False

        if self.outIndv.get() and doInd:
            self.indepConsensusSets = {}
            #Separating clusters by input set
            indepClustersDic = self.getIndepClusters(molClusters, molDic)
            for inSetId in indepClustersDic:
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId] = self.cluster2representative(indepClustersDic[inSetId], molDic, minSize=1)


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
            #newDock = self.relabelPosId(newDock)
            outDocked.append(newDock)
        self._defineOutputs(outputSmallMolecules=outDocked)

        if self.outIndv.get() and hasattr(self, 'indepConsensusSets'):
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

        if not self.repAttr.get() or not self.repAttr.get().strip():
            validations.append('You must specify an attribute to choose a ligand representative of the generated '
                               'cluster. Typically, the energy or the score, to choose the best out of the cluster')

        return validations

    # --------------------------- UTILS functions -----------------------------------
    def getInputMolsDir(self):
        return os.path.abspath(self._getExtraPath('inputMolecules'))

    def buildMolDic(self):
        dic = {}
        for i, pSet in enumerate(self.inputMoleculesSets):
            for mol in pSet.get():
                dic[mol.getUniqueName()] = i
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
        '''Generate the docking clusters based on the RMSD of the ligands
        Return (clusters): [[dock1, dock2], [dock3], [dock4, dock5, dock6]]'''
        clusters = []
        #For each set of molecules
        for i, molSet in enumerate(self.inputMoleculesSets):
            #For each of the molecules in the set
            for newMol in molSet.get():
                newClusters, newClust = [], [newMol.clone()]
                #Check for each of the clusters
                for clust in clusters:
                    #Check for each of the molecules in the cluster
                    append2Cluster = False
                    for cMol in clust:
                        curRMSD = self.calculateMolsRMSD(newMol, cMol)

                        #If the RMSD threshold is passed
                        if curRMSD <= self.maxRMSD.get():
                            append2Cluster = True
                            break

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

    def buildRepresentativeClusters(self, molDic):
        clusters = []
        # For each set of molecules
        for i, molSet in enumerate(self.inputMoleculesSets):
            # For each of the molecules in the set
            for newMol in molSet.get():
                newClusters, newClust = [], [newMol.clone()]
                # Check for each of the clusters
                for clust in clusters:
                    # Check for each of the molecules in the cluster
                    repMol = clust[0] #Only checking first element (representative)
                    repRMSD = self.calculateMolsRMSD(newMol, repMol)
                    # If the RMSD threshold is passed
                    if repRMSD <= self.maxRMSD.get():
                        repScore, newScore = getattr(repMol, self.repAttr.get()), getattr(newMol, self.repAttr.get())
                        if (newScore > repScore and self.maxmin.get()) or \
                                (newScore < repScore and not self.maxmin.get()):
                            # Becomes new representative (first position)
                            newClust += clust
                        else:
                            newClust = clust + newClust
                    else:
                        # newClust: Init with only the newPocket, grow with each clust which overlaps with newPocket
                        newClusters.append(clust)

                # Add new cluster containing newPocket + overlapping previous clusters
                newClusters.append(newClust)
                clusters = newClusters.copy()

        #Getting representatives
        reps = []
        for cl in clusters:
            if self.countMolsInCluster(cl, molDic) >= self.numOfOverlap.get():
                reps.append(cl[0])
        return reps

    def generateDockingClustersScipy(self):
        '''Generate the docking clusters based on the RMSD of the ligands
        Return (clusters): [[dock1, dock2], [dock3], [dock4, dock5, dock6]]'''
        allMols = self.getAllInputMols()
        posArrays, molsDic = {}, {}
        for mol in allMols:
            posDic = mol.getAtomsPosDic()
            keys = ''.join(sorted(list(posDic.keys())))
            if keys in posArrays:
                posArrays[keys].append([x for coord in sorted(posDic.items()) for x in coord[1]])   # flattening coords
                molsDic[keys].append(mol)
            else:
                posArrays[keys] = [[x for coord in sorted(posDic.items()) for x in coord[1]]]  # flattening coords
                molsDic[keys] = [mol]

        finalClusters = []
        for keys in posArrays:
            if len(posArrays[keys]) > 1: #Needed more than 1 mol to cluster
                # Matrix distance and clustering perform over different molecules separatedly
                rmsds = pairwise_distances(posArrays[keys], metric=self.calculateFlattenRMSD,
                                           n_jobs=self.numberOfThreads.get())   # Parallel distance matrix calculation

                linked = linkage(rmsds, self.getEnumText('linkage').lower())
                clusters = fcluster(linked, self.maxRMSD.get(), 'distance')

                clusterMol = {}
                for i, cl in enumerate(clusters):
                    if cl in clusterMol:
                        clusterMol[cl] += [molsDic[keys][i]]
                    else:
                        clusterMol[cl] = [molsDic[keys][i]]
                finalClusters += list(clusterMol.values())
            else:
                finalClusters += [molsDic[keys]]

        return finalClusters

    def calculateFlattenRMSD(self, u, v):
        '''From a flat vector of 3D coords, unflattens it into a matrix and calculates the RMSD of those coordinates'''
        return calculateRMSD(self.unflattenCoords(u), self.unflattenCoords(v))

    def unflattenCoords(self, v):
        return [v[i:i+3] for i in range(0,len(v),3)]


    def getIndepClusters(self, clusters, molDic):
        indepClustersDic = {}
        for clust in clusters:
            if self.countMolsInCluster(clust, molDic) >= self.numOfOverlap.get():
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

    def countMolsInCluster(self, cluster, molDic):
        setIds = []
        for mol in cluster:
            setId = molDic[mol.getUniqueName()]
            if self.sameClust.get() or not setId in setIds:
                setIds.append(setId)
        return len(setIds)

    def cluster2representative(self, clusters, molDic, minSize):
        representatives = []
        for clust in clusters:
            if self.countMolsInCluster(clust, molDic) >= minSize:
                repr = self.getRepresentativeMolecule(clust)
                representatives.append(repr)
        return representatives

    def getRepresentativeMolecule(self, cluster):
        '''Return the docked molecule with max score in a cluster'''
        maxScore = -10e15 if self.maxmin.get() else 10e15
        for mol in cluster:
            molScore = getattr(mol, self.repAttr.get())
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
            rmsd = calculateRMSDKeys(posDic1, posDic2)
            return rmsd
        else:
            return 100


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

    def getConvMolNames(self):
        molDir = self.getInputMolsDir()
        convMolsName  = {}
        for molFile in glob.glob(os.path.join(molDir, '*')):
            convMolsName[getBaseName(molFile)] = molFile
        return convMolsName

    def getAllInputMols(self):
        mols = []
        convMolNames = self.getConvMolNames()
        for molSet in self.inputMoleculesSets:
            for mol in molSet.get():
                newMol = mol.clone()
                if mol.getUniqueName() in convMolNames:
                    newMol.setPoseFile(os.path.relpath(convMolNames[newMol.getUniqueName()]))
                mols.append(newMol)
        return mols






