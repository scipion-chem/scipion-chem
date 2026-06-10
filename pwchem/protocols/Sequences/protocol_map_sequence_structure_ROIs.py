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
This protocol maps a set of sequence ROIs to structure ROIs for an AtomStruct

"""
import json
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.MMCIFParser import MMCIFParser

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from pwchem.utils.utilsFasta import pairwiseAlign
from pwchem import Plugin
from pwchem.constants import MGL_DIC
from pwchem.protocols import ProtDefineStructROIs


class ProtMapSequenceROI(ProtDefineStructROIs):
    """
    AI Generated:

    This protocol maps SequenceROIs (sequence-defined regions of interest)
    onto 3D structural coordinates in an AtomStruct, generating structural
    ROIs (StructROIs).

    It converts linear sequence regions into spatial clusters of atoms or
    surface points in a protein structure, optionally clustering nearby
    coordinates into coherent structural pockets.

    Overview
    --------
    The protocol bridges sequence-based annotations with structural biology
    by projecting ROI definitions from sequence space into 3D space.

    It performs:
    - Sequence ↔ structure alignment
    - Residue index mapping
    - Extraction of atomic coordinates for ROI residues
    - Optional projection onto molecular surface
    - Optional clustering of spatial coordinates
    - Generation of StructROI objects

    Input
    -----
    inputSequenceROIs:
        SetOfSequenceROIs defining regions on a reference sequence.

    inputAtomStruct:
        Structural model (PDB/mmCIF/compatible formats).

    chain_name:
        JSON-like string specifying model and chain selection:
        {
            "chain": "A",
            "model": 0
        }

    Options
    -------
    doCluster:
        Whether to cluster spatial coordinates into discrete structural ROIs.

    maxIntraDistance / maxDepth:
        Parameters controlling spatial clustering and surface projection.

    Workflow
    --------
    1. Extract reference sequence and structure-derived sequence
    2. Perform pairwise alignment (sequence ↔ structure sequence)
    3. Build residue index mapping (sequence → structure residues)
    4. For each ROI:
       - Map sequence indices to structure residues
       - Extract all atomic coordinates for those residues
    5. Optional:
       - Project coordinates to molecular surface (MSMS)
       - Filter by surface proximity threshold
       - Cluster coordinates into spatial groups
    6. Generate pocket files (CIF format)
    7. Create StructROI objects and compute structural contacts

    Output
    ------
    outputStructROIs:
        SetOfStructROIs containing:
        - Structural ROI coordinates (clusters or per-ROI mapping)
        - Pocket files in CIF format
        - Contact analysis metadata

    Key Features
    ------------
    - Sequence-to-structure ROI mapping
    - Support for multi-residue ROIs
    - Atom-level coordinate extraction
    - Optional surface projection using MSMS
    - Distance-based clustering of spatial regions
    - Compatible with multiple structure formats (PDB/mmCIF/Schrodinger)

    Spatial Processing
    ------------------
    1. Residue mapping:
       Sequence indices → structure residue IDs via alignment

    2. Coordinate extraction:
       All atom coordinates from mapped residues

    3. Surface mapping (optional):
       - Compute molecular surface using MSMS
       - Project ROI coordinates to nearest surface points
       - Filter by maximum depth threshold

    4. Clustering (optional):
       - Merge ROI coordinates
       - Group spatially close points into clusters
       - Each cluster becomes a structural ROI

    Outputs per ROI
    ---------------
    Each ROI is converted into:
    - A CIF pocket file
    - A StructROI object
    - Number of spatial points
    - Contact information with structure

    Notes
    -----
    - Requires accurate sequence alignment for correct mapping
    - Missing alignment positions are skipped
    - Chain and model selection must be correctly specified
    - Surface-based projection improves biological interpretability
    """
    _label = 'Map sequence ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input sequence')
        group.addParam('inputSequenceROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       allowsNull=False, label="Input sequence ROIs: ",
                       help='Select the AtomStruct object where the pockets will be defined')

        group = form.addGroup('Input structure')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      allowsNull=False, label="Input AtomStruct: ",
                      help='Select the AtomStruct object where the pockets will be defined')

        group.addParam('chain_name', params.StringParam,
                      allowsNull=False, label='Chain of interest',
                      help='Specify the chain of the residue of interest')

        group.addParam('preview', params.LabelParam,
                       label='Preview alignment: ', condition='chain_name!=""',
                       help='Preview the alignment of the specified Sequence and AtomStruct chain')

        group = form.addGroup('Distances')
        group.addParam('doCluster', params.BooleanParam, default=True,
                       label='Cluster output coordinates: ',
                       help='Whether to cluster the ROI coordinates and extract those clusters as the final '
                            'structural ROIs or define a structural ROI for each sequence ROI input.')
        group = self._defineClusterParams(group, condition='doCluster')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('alignSequencesStep')
        self._insertFunctionStep('definePocketsStep')
        self._insertFunctionStep('defineOutputStep')

        self._mapWarning = 'Mapping of ROIs not posible, check the alignment in ', self._getPath("pairWise.aln")

    def alignSequencesStep(self):
        seq, seqName = self.getInputSequence()
        seqAS, seqNameAS = self.getInputASSequence()

        # Alignment
        out_file = os.path.abspath(self._getPath("pairWise.fasta"))
        pairwiseAlign(seq, seqAS, out_file, seqName1=seqName, seqName2=seqNameAS)
        pairwiseAlign(seq, seqAS, out_file.replace('.fasta', '.aln'), seqName1=seqName, seqName2=seqNameAS)

    def definePocketsStep(self):
        cifFile = cifFromASFile(self.getInputPath(), self._getCifFile(), atomStruct=self.inputAtomStruct.get())

        structModel = MMCIFParser().get_structure(self._getInputName(), cifFile)[0] # 0: modelID?

        self.structSurface = Plugin.runMSMS(structModel)
        mapDic = self.mapResidues(structModel)

        pocketCoords = self.getROICoords(mapDic, structModel)
        if pocketCoords:
            if self.surfaceCoords:
                pocketCoords = self.mapSurfaceCoords(pocketCoords)

            if self.doCluster:
                allCoords = self.mergeCoords(pocketCoords)
                coordsClusters = clusterSurfaceCoords(allCoords, self.maxIntraDistance.get())
                coordsClusters = dict(enumerate(coordsClusters))
            else:
                coordsClusters = pocketCoords

            self.saveCoordClusters(coordsClusters)
        else:
            print(self._mapWarning)

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        coordsClusters = self.readCoordClusters()
        if coordsClusters:
            for i, clust in coordsClusters.items():
                if clust:
                    pocketFile = self._getExtraPath(f'pocketFile_{i}.cif')
                    createPocketFile(clust, i, pocketFile)
                    pocket = StructROI(pocketFile, self._getCifFile())
                    pocket.setNumberOfPoints(len(clust))
                    if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                        pocket._maeFile = String(os.path.relpath(inpStruct.getFileName()))
                    pocket.calculateContacts()
                    outPockets.append(pocket)
                else:
                    print(self._mapWarning)

        if len(outPockets) > 0:
            outPockets.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outPockets)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getInputPath(self):
        return self.inputAtomStruct.get().getFileName()

    def getInputFileName(self):
        return self.getInputPath().split('/')[-1]

    def _getPdbFile(self):
        return os.path.abspath(self._getExtraPath(self._getInputName() + '.pdb'))

    def _getCifFile(self):
        return os.path.abspath(self._getExtraPath(self._getInputName() + '.cif'))

    def _getInputName(self):
        return getBaseName(self.getInputPath())

    def getInputSequence(self):
        inputObj = getattr(self, 'inputSequenceROIs').get()
        seq = inputObj.getSequence()
        seqName = inputObj.getSequenceObj().getId()
        if not seqName:
            seqName = inputObj.getSeqName()
        return seq, seqName

    def getInputASSequence(self):
        from pwem.convert.atom_struct import AtomicStructHandler
        asFile = self.getInputPath()
        seqName = os.path.basename(asFile)
        handler = AtomicStructHandler(asFile)
        chainName = getattr(self, 'chain_name').get()

        # Parse chainName for model and chain selection
        struct = json.loads(chainName)  # From wizard dictionary
        chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

        seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))
        return seq, seqName

    def mapResidues(self, structModel):
        '''Returns a dictionary which maps the idxs of the residues of  the sequence and the sequence from a structure
        idxSeq: idxStrSeq}'''
        seq, seqAS = self.parseSequences()
        chain = structModel[json.loads(self.chain_name.get())['chain']]
        resIdxs = self.getChainResidueIdxs(chain)

        mapDic = {}
        i, j = 0, 0
        for k in range(len(seq)):
            if seq[k] == '-':
                i = i + 1
            if seqAS[k] == '-':
                j = j + 1
            if seq[k] != '-' and seqAS[k] != '-':
                mapDic[k-i+1] = resIdxs[k-j]

        return mapDic

    def mergeCoords(self, coordDic):
        coords = []
        for roiId, roiCoords in coordDic.items():
            for coord in roiCoords:
                if coord not in coords:
                    coords.append(coord)
        return coords

    def getROICoords(self, mapDic, structModel):
        resIdxs = {}
        for roi in self.inputSequenceROIs.get():
            roiId = roi.getObjId()
            resIdxs[roiId] = []
            for roiIdx in range(roi.getROIIdx(), roi.getROIIdx2() + 1):
                if roiIdx in mapDic:
                    #Only added a residue if there was correspondance in the alignment
                    resIdxs[roiId].append(mapDic[roiIdx])

        coords = {}
        chainId = json.loads(self.chain_name.get())['chain']
        for roiId in resIdxs:
            coords[roiId] = []
            for resId in resIdxs[roiId]:
                residue = structModel[chainId][resId]
                atoms = residue.get_atoms()
                for a in atoms:
                    coords[roiId].append(list(a.get_coord()))
        return coords

    def mapSurfaceCoords(self, oCoords):
        sCoords = {}
        for roiId in oCoords:
            sCoords[roiId] = []
            for coord in oCoords[roiId]:
                closerSCoords = self.closerSurfaceCoords(coord)
                for cCoord in closerSCoords:
                    cCoord = list(cCoord)
                    if cCoord not in sCoords[roiId]:
                      sCoords[roiId].append(cCoord)

        return sCoords

    def closerSurfaceCoords(self, coord):
        distances = distance.cdist([coord], self.structSurface)
        closestIndexes = distances < self.maxDepth.get()
        return list(self.structSurface[closestIndexes[0]])

    def parseSequences(self):
        seq, seqAS = '', ''
        first = True
        with open(self._getPath("pairWise.fasta")) as f:
            f.readline()
            for line in f:
                if not line.startswith('>'):
                    if first:
                        seq += line.strip()
                    else:
                        seqAS += line.strip()
                else:
                    first = False
        return seq, seqAS

    def getChainResidueIdxs(self, chain):
        resIdxs = []
        for res in chain:
            resIdxs.append(res.get_id()[1])
        return resIdxs

    def getCoordsFile(self):
        return self._getExtraPath('coordClusters.txt')

    def saveCoordClusters(self, coordClusters):
        with open(self.getCoordsFile(), "w") as f:
            json.dump(coordClusters, f)

    def readCoordClusters(self):
        with open(self.getCoordsFile(), "r") as f:
            coordClusters = json.load(f)
        return coordClusters








