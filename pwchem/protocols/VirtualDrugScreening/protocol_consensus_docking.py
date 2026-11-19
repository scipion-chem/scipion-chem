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
import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pyworkflow.object import Float
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *

MIN, MAX = 0, 1
REF_AS, REF_SET = 0, 1

# Cluster provenance attributes written on every representative (see clusterStats)
STAT_ATTRS = ['consensusMethods', 'consensusPoses', 'consensusScore', 'clusterRMSD']

# Outlier rejection of the receptor superposition (fitAtoms), as in PyMOL align.
FIT_CYCLES, FIT_CUTOFF = 5, 2.0

class ProtocolConsensusDocking(EMProtocol):
    """
    AI Generated:

    This protocol is used to perform consensus docking analysis across multiple
    sets of docked small molecules. It clusters docking poses based on RMSD
    similarity and identifies consensus binding modes supported by multiple
    inputs (either from the same docking run or different docking software).

    The protocol can optionally generate independent consensus sets per input
    molecule set and supports multiple representative selection strategies,
    including medoid-based and attribute-based selection (e.g., docking score).

    Core Concepts
    -------------
    Docked Molecules:
        Small molecules with multiple docking poses associated with a receptor.

    Pose Clustering:
        Docking poses are clustered using hierarchical clustering based on
        pairwise RMSD distances.

    Consensus Cluster:
        A cluster of poses representing a consistent binding mode shared across
        multiple input sets.

    Representative Pose:
        The selected pose that best represents a cluster. It can be chosen as:
        - Medoid (minimum total RMSD to other members)
        - Best/worst value of a user-defined attribute (e.g., docking score)

    Workflow
    --------
    1. Input multiple SetOfSmallMolecules (docked pose sets).
    2. Convert input structures into standardized file formats if needed
       (PDB, MAE, etc.).
    3. Compute pairwise RMSD distances between all docking poses.
    4. Perform hierarchical clustering using a selected linkage method.
    5. Cut clusters based on a maximum RMSD threshold.
    6. Filter clusters by minimum number of supporting poses (consensus size).
    7. Select representative molecule for each cluster.
    8. Optionally generate:
       - Global consensus set
       - Independent consensus sets per input molecule set
    9. Output clustered and representative docking poses.

    Clustering Criteria
    -------------------
    - RMSD threshold (maxRMSD):
        Defines the maximum allowed structural deviation within a cluster.

    - Minimum overlap (numOfOverlap):
        Minimum number of poses required for a cluster to be considered valid.

    - Same-set counting (sameClust):
        Controls whether multiple poses from the same input set count separately
        when evaluating cluster support.

    Representative Selection Modes
    ------------------------------
    Medoid:
        Selects the molecule with the smallest sum of RMSD distances to all
        other cluster members.

    Attribute-based:
        Selects molecule with minimum or maximum value of a specified attribute
        (e.g., docking score, binding energy).

    Internal Processing Logic
    -------------------------
    - Input molecules are grouped and assigned set identifiers.
    - Files are converted to compatible formats when necessary.
    - Pairwise RMSD matrices are computed in parallel for efficiency.
    - Hierarchical clustering is applied using SciPy linkage methods.
    - Flat clusters are extracted using distance cutoff (maxRMSD).
    - Cluster centrality is computed using summed RMSD distances.
    - Clusters are filtered based on consensus constraints.
    - Representative poses are cloned to avoid side effects.

    Output
    ------
    - outputSmallMolecules:
        Set containing consensus representative docking poses.

    - outputSmallMolecules_<XXX> (optional):
        Independent consensus sets per input set when enabled.

    - outputSmallMolecules_All:
        Full consensus docking set including all clustered poses.

    Use Cases
    ---------
    - Identification of stable binding modes across docking runs
    - Cross-docking consensus analysis
    - Reduction of redundant docking poses
    - Prioritization of reproducible ligand–receptor interactions
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

        group = form.addGroup('Receptor alignment')
        group.addParam('doAlign', params.BooleanParam, default=False,
                       label='Align receptors: ',
                       help='Superpose every pose onto a common reference before clustering.\n'
                            'The RMSD used for the clustering is calculated on the coordinates as they are, '
                            'with no fitting, so it is only meaningful when every pose has been generated in '
                            'the same coordinate frame. That is the case for docking into a single prepared '
                            'receptor, but not for co-folding methods, which predict a receptor of their own '
                            'for each model.\nEach pose is moved with the transformation that superposes its '
                            'own receptor onto the reference one, matching CA atoms by chain and residue '
                            'number. The receptor of a pose is its own protein file if it has one '
                            '(co-folding, where every model has its own structure), else the protein file of '
                            'its input set (docking, where a whole set shares the receptor).')
        group.addParam('refOrigin', params.EnumParam, default=REF_AS, condition='doAlign',
                       choices=['AtomStruct', 'From input set'], label='Reference: ', display=params.EnumParam.DISPLAY_HLIST,
                       help='Structure every pose is superposed onto.')
        group.addParam('refStruct', params.PointerParam, pointerClass='AtomStruct',
                       condition='doAlign and refOrigin==%s' % REF_AS, allowsNull=True,
                       label='Reference structure: ',
                       help='Receptor (or complex) used as reference for the superposition.')
        group.addParam('refSet', params.StringParam, default='',
                       condition='doAlign and refOrigin==%s' % REF_SET,
                       label='Reference input set: ',
                       help='The receptor of the first molecule of this input set is used as reference. '
                            'Use the wizard to pick one of the input sets.')
        group.addParam('maxFitRMSD', params.FloatParam, default=5.0, condition='doAlign',
                       expertLevel=params.LEVEL_ADVANCED, label='Warn above fit RMSD: ',
                       help='A superposition whose CA RMSD is above this value is reported in the log. A bad '
                            'fit is indistinguishable from a method that disagrees once the coordinates have '
                            'been written, so it is worth checking.')

        group = form.addGroup('Clustering')
        group.addParam('linkage', params.EnumParam, default=0,
                       choices=['Single', 'Complete', 'Average', 'Ward', 'Centroid', 'Median'],
                       label='Hierarchical clustering linkage: ',
                       help='Hierarchical clustering linkage used on scipy clustering.\n'
                            'https://docs.scipy.org/doc/scipy/reference/generated/scipy.cluster.hierarchy.linkage.html')
        group.addParam('maxRMSD', params.FloatParam, default=1, label='Max RMSD for overlap: ',
                      help="Maximum RMSD for clustering different docked molecules")
        group.addParam('numOfOverlap', params.IntParam, default=2, label='Minimum number of overlapping dockings: ',
                      help="Min number of docked molecules to be considered consensus docking")
        group.addParam('sameClust', params.BooleanParam, default=True,
                      label='Count poses from same input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Whether to count overlapping docked molecules from the same input set when calculating the '
                           'cluster size.\nWhen True, the size of a cluster is its number of poses, so a single input '
                           'set can reach the minimum on its own. When False, it is the number of input sets present '
                           'in the cluster, so the minimum is the number of inputs that must agree.')
        group.addParam('minPosesPerSet', params.IntParam, default=1, condition='not sameClust',
                       expertLevel=params.LEVEL_ADVANCED, label='Min poses per input to count: ',
                       help='How many poses an input set must place in a cluster before that set counts towards '
                            'the cluster size. 1 means a single pose is enough.\nRaise it when each input produces '
                            'several correlated poses (the diffusion samples of a co-folding method, for example), '
                            'so that an input counts with a pose it reproduces and not with its luckiest sample.')

        group = form.addGroup('Representative')
        group.addParam('repMode', params.EnumParam, default=0, choices=['Medoid', 'By attribute'],
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
        group.addParam('oneRepPerMol', params.BooleanParam, default=False,
                       label='One representative per molecule: ',
                       help='Keep only the best supported cluster of each molecule, instead of one representative '
                            'per cluster reaching the minimum size. Clusters are compared by number of input sets, '
                            'then number of poses, then tightness (mean RMSD inside the cluster).\nMolecules with no '
                            'cluster reaching the minimum produce no output.')
        group.addParam('outIndv', params.BooleanParam, default=False,
                      label='Output for each input: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Creates an output set related to each input set, with the elements from each input'
                           'present in the consensus clusters')

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

        if self.doAlign.get():
            self.alignMolecules()

    def consensusStep(self):
        molSetDic = self.buildMolSetDic()
        self.posesPerSetMol = self.buildPosesPerSetMol()
        minSize = self.numOfOverlap.get()

        molClusters = self.generateDockingClusters()
        self.consensusMols = self.cluster2representative(molClusters, molSetDic, minSize)

        if self.oneRepPerMol.get():
            self.consensusMols = self.selectBestPerMol(self.consensusMols)

        if self.outIndv.get():
            self.indepConsensusSets = {}
            #Separating clusters by input set
            indepClustersDic = self.getIndepClusters(molClusters, molSetDic)
            for inSetId in indepClustersDic:
                # Getting independent representative for each input set
                self.indepConsensusSets[inSetId] = self.cluster2representative(indepClustersDic[inSetId],
                                                                               molSetDic, minSize=1,
                                                                               addStats=False)

    def createOutputStep(self):
        inputProteinFile = self.inputMoleculesSets[0].get().getProteinFile()

        self.relabelDic = {}
        self.consensusMols = self.fillEmptyAttributes(self.consensusMols)
        self.consensusMols, idsDic = self.reorderIds(self.consensusMols)

        # Each representative keeps the pose and the receptor of its own input protocol
        origDic, recDic = self.buildOrigPoseDic(), self.buildReceptorDic()

        outDocked = SetOfSmallMolecules(filename=self._getPath('consensusDocked_All.sqlite'))
        outDocked.setDocked(True)
        outDocked.setProteinFile(inputProteinFile)
        for outDock in self.consensusMols:
            newDock = outDock.clone()
            #newDock = self.relabelPosId(newDock)
            newDock = self.restoreOriginalFiles(newDock, origDic, recDic)
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

        if not self.sameClust.get() and self.minPosesPerSet.get() < 1:
            validations.append('The minimum number of poses per input to count must be at least 1.')

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

    def restoreOriginalFiles(self, mol, origDic, recDic):
        '''Point a representative back at the files of the protocol that generated it: its pose as it
        was written and its own receptor, which are already in the same frame.'''
        posFile = origDic.get(mol.getUniqueName())
        if posFile:
            mol.setPoseFile(self.getRelPath(posFile))

        recFile = recDic.get(mol.getUniqueName())
        if recFile:
            mol.setProteinFile(self.getRelPath(recFile))
        return mol

    def buildOrigPoseDic(self):
        '''{uniqueName: pose file as the input protocol wrote it}'''
        return {mol.getUniqueName(): mol.getPoseFile()
                for molSet in self.inputMoleculesSets for mol in molSet.get()}

    def getRelPath(self, filePath):
        '''Path relative to the project, so that the output does not depend on where it is stored'''
        projPath = os.path.abspath(self.getProject().getPath())
        absPath = os.path.abspath(filePath)
        return os.path.relpath(absPath, projPath) if absPath.startswith(projPath) else filePath

    ################## Receptor alignment ##################
    def alignMolecules(self):
        '''Move every converted pose with the transformation superposing its own receptor onto the
        reference one. Done on the converted files, so the clustering RMSD, which does not fit
        anything, compares poses living in the same coordinate frame'''
        recDic, molsDir = self.buildReceptorDic(), self.getInputMolsDir()
        self.transDic, self.fitDic, self.nFitDic, self.nCAs = {}, {}, {}, {}

        for molFile in glob.glob(os.path.join(molsDir, '*')):
            uName = getBaseName(molFile)
            if uName not in recDic:
                raise Exception(f'No input molecule found for the converted file {molFile}')

            rot, tran = self.getTransform(recDic[uName])
            self.transformMolFile(molFile, rot, tran)

        self.reportFits()

    def reportFits(self):
        '''Summary of the receptor superpositions. A fit is reported as bad when its RMSD is high or when
        most of the paired CAs had to be discarded as outliers (the RMSD alone only covers the kept ones)'''
        if not self.fitDic:
            return
        rmsds = np.array(list(self.fitDic.values()))
        print(f'Receptor alignment: {len(rmsds)} receptor(s) fitted, CA RMSD mean {rmsds.mean():.2f} A, '
              f'max {rmsds.max():.2f} A')

        bad = {rec: (rms, *self.nFitDic[rec]) for rec, rms in self.fitDic.items()
               if rms > self.maxFitRMSD.get() or self.nFitDic[rec][0] < 0.5 * self.nFitDic[rec][1]}
        if bad:
            print(f'WARNING: {len(bad)} receptor(s) badly superposed (CA RMSD above {self.maxFitRMSD.get()} A '
                  f'or less than half of the paired CAs within {FIT_CUTOFF} A). Worst ones:')
            for rec, (rms, nUsed, nPairs) in sorted(bad.items(), key=lambda x: -x[1][0])[:3]:
                print(f'  {os.path.basename(rec)}: RMSD {rms:.2f} A, {nUsed}/{nPairs} CAs used, '
                      f'CAs in file {self.nCAs.get(rec)} vs reference {self.nCAs.get(self.getReferenceFile())}')


    def getTransform(self, recFile):
        '''(rotation, translation) taking recFile onto the reference. Cached by receptor file, so a
        receptor shared by a whole input set is only fitted once'''
        recFile = os.path.abspath(recFile)
        if recFile not in self.transDic:
            refChains, movChains = self.getReferenceChains(), self.getChains(recFile)
            refCAs, movCAs = self.pairCAs(refChains, movChains)
            if len(refCAs) < 3:
                raise Exception(f'Only {len(refCAs)} CA atoms could be paired between {recFile} and the '
                                f'reference {self.getReferenceFile()}. Check that both structures are the '
                                f'same protein.')

            rotran, rms, nUsed = self.fitAtoms(refCAs, movCAs)
            self.transDic[recFile], self.fitDic[recFile] = rotran, rms
            self.nFitDic[recFile] = (nUsed, len(refCAs))
            self.nCAs[recFile] = sum(len(c) for c in self.getChains(recFile).values())
        return self.transDic[recFile]

    def fitAtoms(self, refAtoms, movAtoms):
        '''Superpose movAtoms onto refAtoms, rejecting the pairs further apart than FIT_CUTOFF and
        fitting again, as PyMOL align does. Returns (rotran, RMSD, atoms used in the last cycle)'''
        sup = Superimposer()
        sup.set_atoms(refAtoms, movAtoms)

        for _ in range(FIT_CYCLES):
            rot, tran = sup.rotran
            dists = [np.linalg.norm(np.dot(mov.coord, rot) + tran - ref.coord)
                     for ref, mov in zip(refAtoms, movAtoms)]
            keep = [i for i, dist in enumerate(dists) if dist <= FIT_CUTOFF]
            if len(keep) < 3 or len(keep) == len(refAtoms):
                break

            refAtoms = [refAtoms[i] for i in keep]
            movAtoms = [movAtoms[i] for i in keep]
            sup = Superimposer()
            sup.set_atoms(refAtoms, movAtoms)

        return sup.rotran, sup.rms, len(refAtoms)

    def pairCAs(self, refChains, movChains):
        '''Pair the CA atoms of two structures chain by chain, in the order they appear in the file.
        Residue numbering is not used, so a model numbered from 1 can be superposed on a reference
        with its own numbering, as long as both contain the same residues in the same order'''
        chainPairs = [(cId, cId) for cId in refChains if cId in movChains]
        if not chainPairs:
            chainPairs = list(zip(refChains.keys(), movChains.keys()))

        refCAs, movCAs = [], []
        for refId, movId in chainPairs:
            refAtoms, movAtoms = refChains[refId], movChains[movId]

            refCAs += [atom for atom, _ in zip(refAtoms, movAtoms)]
            movCAs += [atom for _, atom in zip(refAtoms, movAtoms)]
        return refCAs, movCAs

    def getRefSetId(self):
        """Return the index of the input set whose name matches refSet"""
        targetName = self.refSet.get()
        if targetName:
            for idx, setPointer in enumerate(self.inputMoleculesSets):
                if str(setPointer.get()) == targetName:
                    return idx
        return None

    def getReferenceFile(self):
        if self.refOrigin.get() == REF_AS:
            return os.path.abspath(self.refStruct.get().getFileName())

        molSet = self.inputMoleculesSets[self.getRefSetId()].get()
        refFile = molSet.getFirstItem().getProteinFile() or molSet.getProteinFile()
        return os.path.abspath(refFile)

    def getReferenceChains(self):
        if not hasattr(self, 'refChains'):
            self.refChains = self.getChains(self.getReferenceFile())
            self.nCAs[self.getReferenceFile()] = sum(len(c) for c in self.refChains.values())
        return self.refChains

    def getChains(self, structFile):
        '''{chainId: [CA atoms]} of the first model of a structure, in file order'''
        parser = MMCIFParser(QUIET=True) if structFile.lower().endswith(('.cif', '.mmcif')) \
            else PDBParser(QUIET=True)
        struct = parser.get_structure(getBaseName(structFile), structFile)

        chains = {}
        for model in struct:
            for chain in model:
                atoms = [residue['CA'] for residue in chain if 'CA' in residue]
                if atoms:
                    chains[chain.id] = atoms
            break
        return chains

    def transformMolFile(self, molFile, rot, tran):
        '''Apply a rotation and a translation to the coordinates of a (multi)molecule SDF file'''
        if not molFile.lower().endswith(('.sdf', '.sd')):
            raise Exception(f'Receptor alignment is only implemented for SDF poses, got {molFile}')

        with open(molFile) as f:
            lines = f.readlines()

        outLines, i = [], 0
        while i < len(lines):
            countsIdx = i + 3
            if countsIdx >= len(lines):
                outLines += lines[i:]
                break

            if 'V3000' in lines[countsIdx]:
                raise Exception(f'V3000 SDF files are not supported by the receptor alignment: {molFile}')

            nAtoms = int(lines[countsIdx][0:3])
            outLines += lines[i:countsIdx + 1]

            for atomLine in lines[countsIdx + 1: countsIdx + 1 + nAtoms]:
                coords = [float(atomLine[0:10]), float(atomLine[10:20]), float(atomLine[20:30])]
                newCoords = np.dot(coords, rot) + tran
                outLines.append('%10.4f%10.4f%10.4f%s' % (*newCoords, atomLine[30:]))

            i = countsIdx + 1 + nAtoms
            while i < len(lines):
                outLines.append(lines[i])
                i += 1
                if outLines[-1].startswith('$$$$'):
                    break

        tmpFile = molFile + '.tmp'
        with open(tmpFile, 'w') as f:
            f.writelines(outLines)
        os.replace(tmpFile, molFile)

    def buildMolSetDic(self):
        dic = {}
        for i, pSet in enumerate(self.inputMoleculesSets):
            for mol in pSet.get():
                dic[mol.getUniqueName()] = i
        return dic

    def buildPosesPerSetMol(self):
        '''{(molName, inputSetId): number of poses of that molecule in that input set}.
        Denominator of the consensusScore'''
        dic = {}
        for i, pSet in enumerate(self.inputMoleculesSets):
            for mol in pSet.get():
                key = (mol.getMolName(), i)
                dic[key] = dic.get(key, 0) + 1
        return dic

    def buildReceptorDic(self):
        '''{uniqueName: receptor file}. The receptor of a pose is its own protein file if it has one
        (every co-folding model predicts its own structure), else the one of its input set'''
        recDic = {}
        for molSet in self.inputMoleculesSets:
            molSet = molSet.get()
            setRec = molSet.getProteinFile()
            for mol in molSet:
                recDic[mol.getUniqueName()] = mol.getProteinFile() or setRec
        return recDic

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
                finalClusters.append([(molDic[molName][0], 0)])

        return finalClusters


    def getIndepClusters(self, clusters, molSetDic):
        indepClustersDic = {}
        for clust in clusters:
            if self.countMolsInCluster(clust, molSetDic) >= self.numOfOverlap.get():
                curIndepCluster = {}
                for mol, sRMSD in clust:
                    inSetId = molSetDic[mol.getUniqueName()]
                    if inSetId not in curIndepCluster:
                        curIndepCluster[inSetId] = []
                    curIndepCluster[inSetId] += [(mol, sRMSD)]

                for inSetId, setMols in curIndepCluster.items():
                    if inSetId not in indepClustersDic:
                        indepClustersDic[inSetId] = []
                    indepClustersDic[inSetId] += [setMols]
        return indepClustersDic

    def getSetCounts(self, cluster, molSetDic):
        '''Return {inputSetId: number of poses of that input set in the cluster}'''
        counts = {}
        for mol, _ in cluster:
            setId = molSetDic[mol.getUniqueName()]
            counts[setId] = counts.get(setId, 0) + 1
        return counts

    def countMolsInCluster(self, cluster, molSetDic):
        '''Size of a cluster: either its number of poses (sameClust) or the number of input sets
        placing at least minPosesPerSet poses in it'''
        counts = self.getSetCounts(cluster, molSetDic)
        if self.sameClust.get():
            return sum(counts.values())

        minPoses = self.minPosesPerSet.get()
        return len([setId for setId, nPoses in counts.items() if nPoses >= minPoses])

    def clusterStats(self, cluster, molSetDic):
        '''Support of a cluster, to be stored on its representative:
        - consensusMethods: number of input sets present in it
        - consensusPoses: number of poses in it
        - consensusScore: sum over the input sets of the fraction of their poses (for that molecule)
                          present in the cluster. Ranges from 0 to the number of input sets
        - clusterRMSD: mean RMSD between the poses of the cluster. The sums of distances are already
                       calculated for the medoid, and each pairwise distance appears in two of them
        '''
        counts = self.getSetCounts(cluster, molSetDic)
        molName, nPoses = cluster[0][0].getMolName(), len(cluster)

        score = 0
        for setId, nSetPoses in counts.items():
            score += nSetPoses / self.posesPerSetMol[(molName, setId)]

        meanRMSD = sum([sRMSD for _, sRMSD in cluster]) / (nPoses * (nPoses - 1)) if nPoses > 1 else 0
        return {'consensusMethods': len(counts), 'consensusPoses': nPoses,
                'consensusScore': score, 'clusterRMSD': meanRMSD}

    def cluster2representative(self, clusters, molSetDic, minSize, addStats=True):
        representatives = []
        for clust in clusters:
            if self.countMolsInCluster(clust, molSetDic) >= minSize:
                repr = self.getRepresentativeMolecule(clust)
                if addStats:
                    for attr, value in self.clusterStats(clust, molSetDic).items():
                        setattr(repr, attr, Float(value))
                representatives.append(repr)
        return representatives

    def selectBestPerMol(self, representatives):
        '''Keep a single representative per molName: the one coming from the best supported cluster.
        Clusters are compared by number of input sets, then number of poses, then tightness'''
        bestDic = {}
        for mol in representatives:
            molName = mol.getMolName()
            rank = (mol.consensusMethods.get(), mol.consensusPoses.get(), -mol.clusterRMSD.get())
            if molName not in bestDic or rank > bestDic[molName][0]:
                bestDic[molName] = (rank, mol)
        return [mol for _, mol in bestDic.values()]

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






