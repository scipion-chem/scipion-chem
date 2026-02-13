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

MIN, MAX = 0, 1

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
        group.addParam('linkage', params.EnumParam, default=0,
                       choices=['Single', 'Complete', 'Average', 'Ward', 'Centroid', 'Median'],
                       label='Hierarchical clustering linkage: ',
                       help='Hierarchical clustering linkage used on scipy clustering.\n'
                            'https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html')
        form.addParam('outIndv', params.BooleanParam, default=False,
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
        group.addParam('repMode', params.EnumParam, default=0, choices=['Centroid', 'By attribute'],
                       label='Representative choice mode: ',
                       help='How to choose the representative from the clusters')
        group.addParam('repAttr', params.StringParam, default='', condition='repMode==1',
                      label='Attribute to choose cluster representative: ',
                      help='Criteria attribute  to follow on docking clusters to choose a representative. '
                           'It will extract the representative as the pose with max/min (next argument) value '
                           'of this attribute')
        group.addParam('maxmin', params.EnumParam, default=MAX, choices=['Minimum', 'Maximum'],
                       label='Keep values: ', condition='repMode==1',
                       help='Whether to keep the minimum or maximum values of the criteria attribute to select the '
                            'representative')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        cStep = self._insertFunctionStep(self.convertInputStep, prerequisites=[])
        conStep = self._insertFunctionStep(self.consensusStep, prerequisites=[cStep])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[conStep])

    def convertInputStep(self):
        allMols = self.getAllInputMols()
        outDir = self.getInputMolsDir()
        os.mkdir(outDir)

        tmpDirPdb, tmpDirMae = (os.path.abspath(self._getTmpPath('pdbFiles')),
                                os.path.abspath(self._getTmpPath('maeFiles')))
        os.mkdir(tmpDirPdb), os.mkdir(tmpDirMae)
        for mol in allMols:
            molFile = os.path.abspath(mol.getPoseFile())
            _, ext = os.path.splitext(molFile)
            if ext in ('.pdbqt', '.pdb'):
                iDir = tmpDirPdb
            elif ext in ('.mae', '.maegz'):
                iDir = tmpDirMae
            else:
                iDir = outDir
            shutil.copy(molFile, os.path.join(iDir, mol.getUniqueName() + ext))
        self.convertPDBFiles(tmpDirPdb, outDir)
        self.convertMAEFiles(tmpDirMae, outDir)

    def consensusStep(self):
        molSetDic = self.buildMolSetDic()
        minSize = self.numOfOverlap.get()

        molClusters = self.generateDockingClusters()
        self.consensusMols = self.cluster2representative(molClusters, molSetDic, minSize)

        if self.outIndv.get():
            self.indepConsensusSets = {}
            #Separating clusters by input set
            indepClustersDic = self.getIndepClusters(molClusters, molSetDic)
            for inSetId in indepClustersDic:
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId] = self.cluster2representative(indepClustersDic[inSetId],
                                                                               molSetDic, minSize=1)

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

        if self.repMode.get() == 1 and (not self.repAttr.get() or not self.repAttr.get().strip()):
            validations.append('You must specify an attribute to choose a ligand representative of the generated '
                               'cluster. Typically, the energy or the score, to choose the best out of the cluster')

        return validations

    # --------------------------- UTILS functions -----------------------------------
    def convertPDBFiles(self, pdbDir, outDir):
        pdbFiles = list(os.listdir(pdbDir))
        if len(pdbFiles) > 0:
            args = f' --multiFiles -iD "{pdbDir}" --pattern "*" -of sdf --outputDir "{outDir}"'
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

    def convertMAEFiles(self, maeDir, outDir):
        try:
            from pwchemSchrodinger.utils.utils import convertMAE2Mol2File
        except ImportError:
            print('Conversion of MAE input files could not be performed because schrodinger plugin is not installed')
            return

        maeFiles = list(os.listdir(maeDir))
        if len(maeFiles) > 0:
            for maeFile in maeFiles:
                convertMAE2Mol2File(maeFile, outDir)

    def getInputMolsDir(self):
        return os.path.abspath(self._getExtraPath('inputMolecules'))

    def buildMolSetDic(self):
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

    def sumDistancesVectorized(self, clusters, flatDistances):
        """
        Vectorized version for better performance with larger datasets.
        """
        n = len(clusters)

        # Create full distance matrix
        distMatrix = np.zeros((n, n))
        triuIndices = np.triu_indices(n, k=1)
        distMatrix[triuIndices] = flatDistances
        distMatrix = distMatrix + distMatrix.T

        sums = []
        for clusterId in set(clusters):
            mask = np.array(clusters) == clusterId
            clusterIndices = np.where(mask)[0]

            for i in clusterIndices:
                # Sum distances to all other elements in same cluster
                sumDist = np.sum(distMatrix[i, clusterIndices]) - distMatrix[i, i]
                sums.append(sumDist)

        return sums

    def generateDockingClusters(self):
        '''Generate the docking clusters based on the RMSD of the ligands
        Return the clusters of molecules, where each molecule item is a tuple with the mol object and its associated
        of RMSDs to evaluate its centrality in the cluster.

        Return (clusters): [[(dock1, sRMSD1), (dock2, sRMSD2)], [(dock3, sRMSD3)],...]'''
        allMols = self.getAllInputMols(conv=True)
        _, molDic = buildMolDic(allMols)

        finalClusters = []
        rmsdDic = runParallelRdkitRMSD(allMols, nJobs=self.numberOfThreads.get())
        for molName, rmsds in rmsdDic.items():
            if rmsds: #Needed more than 1 mol to cluster
                linked = linkage(rmsds, self.getEnumText('linkage').lower())
                clusters = fcluster(linked, self.maxRMSD.get(), 'distance')
                rmsdSums = self.sumDistancesVectorized(clusters, rmsds)

                clusterMol = {}
                for i, cl in enumerate(clusters):
                    if cl in clusterMol:
                        clusterMol[cl] += [(molDic[molName][i], rmsdSums[i])]
                    else:
                        clusterMol[cl] = [(molDic[molName][i], rmsdSums[i])]
                finalClusters += list(clusterMol.values())
            else:
                finalClusters += [(molDic[molName][0], 0)]

        return finalClusters


    def getIndepClusters(self, clusters, molSetDic):
        indepClustersDic = {}
        for clust in clusters:
            if self.countMolsInCluster(clust, molSetDic) >= self.numOfOverlap.get():
                curIndepCluster = {}
                for mol, sRMSD in clust:
                    curMolFile = mol.getPoseFile()
                    inSetId = molSetDic[curMolFile]
                    if inSetId not in curIndepCluster:
                        curIndepCluster[inSetId] = []
                    curIndepCluster[inSetId] += [(mol, sRMSD)]

                for inSetId, setMols in curIndepCluster.items():
                    if inSetId not in indepClustersDic:
                        indepClustersDic[inSetId] = []
                    indepClustersDic[inSetId] += [setMols]
        return indepClustersDic

    def countMolsInCluster(self, cluster, molSetDic):
        setIds = []
        for mol, _ in cluster:
            setId = molSetDic[mol.getUniqueName()]
            if self.sameClust.get() or not setId in setIds:
                setIds.append(setId)
        return len(setIds)

    def cluster2representative(self, clusters, molSetDic, minSize):
        representatives = []
        for clust in clusters:
            if self.countMolsInCluster(clust, molSetDic) >= minSize:
                repr = self.getRepresentativeMolecule(clust)
                representatives.append(repr)
        return representatives

    def getRepresentativeMolecule(self, cluster):
        '''Return the docked molecule with max score in a cluster'''
        if self.repMode.get() == 0:
            medoid = min(cluster, key=lambda x: x[1])[0]
            outMol = medoid.clone()
        else:
            maxScore = -10e15 if self.maxmin.get() == MAX else 10e15
            for mol, _ in cluster:
                molScore = getattr(mol, self.repAttr.get())
                if (molScore > maxScore and self.maxmin.get() == MAX) or \
                        (molScore < maxScore and self.maxmin.get() == MIN):
                    maxScore = molScore
                    outMol = mol.clone()
        return outMol


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

    def getAllInputMols(self, conv=False):
        if conv:
            convMolsDic = self.getConvMolsDic()

        mols = []
        for molSet in self.inputMoleculesSets:
            for mol in molSet.get():
                newMol = mol.clone()
                if conv:
                    newMol.setPoseFile(convMolsDic[newMol.getUniqueName()])
                mols.append(newMol)
        return mols

    def getConvMolsDic(self):
        convMolsDic = {}
        inDir = self.getInputMolsDir()
        for molSet in self.inputMoleculesSets:
            for mol in molSet.get():
                molUName = mol.getUniqueName()
                molFile = glob.glob(os.path.join(inDir, f"{molUName}.*"))[0]
                convMolsDic[molUName] = molFile
        return convMolsDic






