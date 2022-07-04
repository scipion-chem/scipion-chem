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
import os, json
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, min_dist, residue_depth
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb, alignClustalSequences

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from pwchem.utils.utilsFasta import pairwiseAlign
from pwchem import Plugin
from pwchem.constants import MGL_DIC


class ProtMapSequenceROI(EMProtocol):
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
        group.addParam('maxDepth', params.FloatParam, default='3.0',
                       label='Maximum atom depth (A): ',
                       help='Maximum atom distance to the surface to be considered and mapped')
        group.addParam('maxIntraDistance', params.FloatParam, default='2.0',
                       label='Maximum distance between pocket points (A): ',
                       help='Maximum distance between two pocket atoms to considered them same pocket')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('alignSequencesStep')
        self._insertFunctionStep('definePocketsStep')
        self._insertFunctionStep('defineOutputStep')

    def alignSequencesStep(self):
        seq, seq_name = self.getInputSequence()
        seqAS, seq_nameAS = self.getInputASSequence()

        # Alignment
        out_file = os.path.abspath(self._getPath("pairWise.fasta"))
        pairwiseAlign(seq, seqAS, out_file, seqName1=seq_name, seqName2=seq_nameAS)
        pairwiseAlign(seq, seqAS, out_file.replace('.fasta', '.aln'), seqName1=seq_name, seqName2=seq_nameAS)

    def definePocketsStep(self):
        parser = PDBParser()
        structModel = parser.get_structure(self.getASName(), self.getASFileName())[0] # 0: modelID?
        self.structSurface = get_surface(structModel,
                                         MSMS=Plugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))

        mapDic = self.mapResidues(structModel)

        originCoords = self.getROICoords(mapDic, structModel)
        if originCoords:
            surfaceCoords = self.mapSurfaceCoords(originCoords)
            self.coordsClusters = self.clusterSurfaceCoords(surfaceCoords)
        else:
            print('Mapping of ROIs not posible, check the alignment in ', self._getPath("pairWise.aln"))

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        if self.coordsClusters:
            outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
            for i, clust in enumerate(self.coordsClusters):
                pocketFile = self.createPocketFile(clust, i)
                pocket = StructROI(pocketFile, self.getASFileName())
                pocket.setNumberOfPoints(len(clust))
                if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                    pocket._maeFile = String(os.path.abspath(inpStruct.getFileName()))
                pocket.calculateContacts()
                outPockets.append(pocket)

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
    def getASFileName(self):
        inpStruct = self.inputAtomStruct.get()
        inpFile = inpStruct.getFileName()
        basename = os.path.basename(inpFile).split('.')[0]
        if inpFile.endswith('.cif'):
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            cifToPdb(inpFile, inpPDBFile)

        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            inpStruct.convert2PDB(outPDB=inpPDBFile)

        elif inpFile.endswith('.pdbqt'):
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(inpFile), inpPDBFile)
            runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        else:
            inpPDBFile = inpFile

        return inpPDBFile

    def getASName(self):
      return os.path.splitext(os.path.basename(self.getASFileName()))[0]

    def getInputSequence(self):
        inputObj = getattr(self, 'inputSequenceROIs').get()
        seq = inputObj.getSequence()
        seq_name = inputObj.getSequenceObj().getId()
        if not seq_name:
            seq_name = inputObj.getSeqName()
        return seq, seq_name

    def getInputASSequence(self):
        from pwem.convert.atom_struct import AtomicStructHandler
        inputObj = getattr(self, 'inputAtomStruct').get()
        ASFile = self.getASFileName()
        seq_name = os.path.basename(ASFile)
        handler = AtomicStructHandler(ASFile)
        chainName = getattr(self, 'chain_name').get()

        # Parse chainName for model and chain selection
        struct = json.loads(chainName)  # From wizard dictionary
        chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

        seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))
        return seq, seq_name

    def mapResidues(self, structModel):
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

    def getROICoords(self, mapDic, structModel):
        resIdxs = []
        for roi in self.inputSequenceROIs.get():
            for roiIdx in range(roi.getROIIdx(), roi.getROIIdx2() + 1):
                if roiIdx in mapDic:
                    #Only added a residue if there was correspondance in the alignment
                    resIdxs.append(mapDic[roiIdx])

        coords = []
        chainId = json.loads(self.chain_name.get())['chain']
        for resId in resIdxs:
            residue = structModel[chainId][resId]
            atoms = residue.get_atoms()
            for a in atoms:
                coords.append(list(a.get_coord()))
        return coords

    def mapSurfaceCoords(self, oCoords):
        sCoords = []
        for coord in oCoords:
            closerSCoords = self.closerSurfaceCoords(coord)
            for cCoord in closerSCoords:
                cCoord = list(cCoord)
                if not cCoord in sCoords:
                  sCoords.append(cCoord)

        return sCoords

    def closerSurfaceCoords(self, coord):
        distances = distance.cdist([coord], self.structSurface)
        closestIndexes = distances < self.maxDepth.get()
        return list(self.structSurface[closestIndexes[0]])

    def clusterSurfaceCoords(self, surfCoords):
        clusters = []
        for coord in surfCoords:
            newClusters = []
            newClust = [coord]
            for clust in clusters:
                merge = False
                for cCoord in clust:
                    dist = calculateDistance(coord, cCoord)
                    if dist < self.maxIntraDistance.get():
                        merge = True
                        break

                if merge:
                    newClust += clust
                else:
                    newClusters.append(clust)

            newClusters.append(newClust)
            clusters = newClusters.copy()
        return clusters

    def createPocketFile(self, clust, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(clust):
                f.write(writePDBLine(['HETATM', str(j), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
        return outFile

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







