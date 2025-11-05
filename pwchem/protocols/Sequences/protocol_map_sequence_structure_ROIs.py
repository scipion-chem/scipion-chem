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
    Maps a set of SequenceROIs to their respective structure ROIs in an AtomStruct
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
        seq, seq_name = self.getInputSequence()
        seqAS, seq_nameAS = self.getInputASSequence()

        # Alignment
        out_file = os.path.abspath(self._getPath("pairWise.fasta"))
        pairwiseAlign(seq, seqAS, out_file, seqName1=seq_name, seqName2=seq_nameAS)
        pairwiseAlign(seq, seqAS, out_file.replace('.fasta', '.aln'), seqName1=seq_name, seqName2=seq_nameAS)

    def definePocketsStep(self):
        cifFile = cifFromASFile(self.getInputPath(), self._getCifFile(), AS=self.inputAtomStruct.get())

        structModel = MMCIFParser().get_structure(self._getInputName(), cifFile)[0] # 0: modelID?
        self.structSurface = get_surface(structModel,
                                         MSMS=Plugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))

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
        seq_name = inputObj.getSequenceObj().getId()
        if not seq_name:
            seq_name = inputObj.getSeqName()
        return seq, seq_name

    def getInputASSequence(self):
        from pwem.convert.atom_struct import AtomicStructHandler
        ASFile = self.getInputPath()
        seq_name = os.path.basename(ASFile)
        handler = AtomicStructHandler(ASFile)
        chainName = getattr(self, 'chain_name').get()

        # Parse chainName for model and chain selection
        struct = json.loads(chainName)  # From wizard dictionary
        chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

        seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))
        return seq, seq_name

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








