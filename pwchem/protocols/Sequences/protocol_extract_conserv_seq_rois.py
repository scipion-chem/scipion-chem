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
This protocol is used define sequence ROIs based on the conservation in a multiple sequence alignment

"""
import os, math, json
import numpy as np
from scipy.stats import entropy

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toPdb, toCIF, AtomicStructHandler, addScipionAttribute

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *
from pwchem.utils.utilsFasta import pairwiseAlign
from pwchem import Plugin

SHANNON, SIMPSON, KABAT, PROP = 'Shannon Entropy', 'Simpson Diversity Index', 'Wu-kabat Variability coefficient', \
                                'Maximum Proportion'

def simpson(counts):
  sumi = 0
  N = sum(counts)
  for c in counts:
    sumi += c * (c - 1)
  return 1 - (sumi / (N * (N - 1)))


def kabat(counts):
  return (sum(counts) * len(counts)) / max(counts)

class ProtExtractSeqsROI(EMProtocol):
    """
    Extract a set of ROIs from a set of sequences based on the conservation in a alignment
    """
    _label = 'Extract sequences ROIs'
    _ATTRNAME = 'Conservation'
    _OUTNAME = 'outputAtomStruct'
    _possibleOutputs = {_OUTNAME: AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSequences', params.PointerParam, pointerClass='SetOfSequences',
                      allowsNull=False, label="Input aligned sequences: ",
                      help='Select the set of sequences object where the ROI will be defined')
        group.addParam('useCons', params.BooleanParam, default=True,
                       label='Use consensus sequence for output: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Use consensus sequence for the sequence ROIs output or a specific input sequence')
        group.addParam('outSeq', params.StringParam, default='', condition='not useCons',
                       label='Output sequence: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Input sequence to use for the output sequence ROIs')

        group.addParam('inputAS', params.PointerParam, pointerClass='AtomStruct',
                       allowsNull=True, label="Input protein structure: ",
                       help='Input protein structure. If included, the conservation values will be drawed over'
                            'the surface of  the protein in the results analysis.')
        group.addParam('chain_name', params.StringParam,
                       allowsNull=True, label='Chain of interest:', condition='inputAS',
                       help='Specify the chain of interest')

        group = form.addGroup('Variability measure')
        group.addParam('method', params.EnumParam, default=3,
                       choices=[SHANNON, SIMPSON, KABAT, PROP],
                       label='Method to measure variability: ',
                       help='Method to measure the conservation / variability in each position of the sequences.\n'
                            'http://imed.med.ucm.es/PVS/pvs-help.html#vmth')

        group = form.addGroup('Cluster regions')
        group.addParam('direction', params.EnumParam, default=0,
                       choices=['Low variability', 'High variability'],
                       label='Look for regions with: ', display=params.EnumParam.DISPLAY_HLIST,
                       help='Whether to look for regions with high or low variability.'
                            'high variability will look for high variability values of the metrics and viceversa')

        group.addParam('highLabel', params.LabelParam, condition='direction==0',
                       label='Keep below threshold',
                       help='High variability will look for high conservation values of the metrics')
        group.addParam('lowLabel', params.LabelParam, condition='direction==1',
                       label='Keep over threshold',
                       help='Low variability will look for low conservation values of the metrics')

        group.addParam('thres', params.FloatParam,
                       label='Main threshold for conserved region: ',
                       help='Main threshold that checks the values of the conservation and defines that a position is '
                            'or is not conserved depending on its conservation values.\n'
                            'Conservation values are normalized from 0 (high conservation) to 1 (low conservation)')
        group.addParam('flexThres', params.FloatParam,
                       label='Flexible threshold for conserved region: ', allowsNull=True,
                       help='This threshold allows an already started conserved region to keep growing if the '
                            'conservation values passes it even if it does not pass the main threshold.\n'
                            'Conservation values are normalized from 0 (high conservation) to 1 (low conservation)')
        group.addParam('numFlex', params.IntParam, default=1,
                       label='Number of consecutive flexible: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Number of consecutive times that the main threshold can be violated and the flexible '
                            'threshold stands.')

        group.addParam('minSize', params.IntParam, default=1,
                       label='Minimum region size: ',
                       help='Minimum size for a region to be considered. Gaps do not count.')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculateConservationStep')
        if self.inputAS.get():
            self._insertFunctionStep('mapConservationStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateConservationStep(self):
        outFile = self.calcConservation()

    def mapConservationStep(self):
        oriSeqsFile = self._getExtraPath('originalAlignment.fasta')
        addSeqsFile = self._getExtraPath('structureAddAlignment.fasta')
        self.inputSequences.get().exportToFile(oriSeqsFile)
        seqAS, seq_nameAS, seqFile = self.getInputASSequence()

        cline = 'mafft --add {} --mapout {} > {}'.format(seqFile, oriSeqsFile, addSeqsFile)
        self.runJob(cline, '')
        if os.path.getsize(addSeqsFile) == 0:
            auxAlignFile = self._getExtraPath('mafftAlignment.fasta')
            print('Alignment of the structure seuqnece failed, most likely because mafft detected that input sequences '
                  'were not aligned.\nAligning sequences with mafft in {}'.format(auxAlignFile))
            cline = 'mafft --auto {} > {}'.format(oriSeqsFile, auxAlignFile)
            self.runJob(cline, '')

            cline = 'mafft --add {} --mapout {} > {}'.format(seqFile, auxAlignFile, addSeqsFile)
            self.runJob(cline, '')

    def defineOutputStep(self):
        if self.useCons:
            outStr = str(self.calcConsensus())
            outSeq = self.inputSequences.get().getFirstItem()
            outSeq.setSequence(outStr)
            outSeq.setSeqName('ConsensusSequence')
        else:
            chosenIdx = int(self.outSeq.get().split(')')[0]) - 1
            for i, seq in enumerate(self.inputSequences.get()):
              if i == chosenIdx:
                  outSeq, outStr = seq, seq.getSequence()
                  break

        roiIdxs = self.getROIs()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for idxs in roiIdxs:
            roi = outStr[idxs[0]-1:idxs[1]-1]
            if len(roi.replace('-', '')) >= self.minSize.get():
                roiSeq = Sequence(sequence=roi, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs))
                seqROI = SequenceROI(sequence=outSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
                outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)

        if self.inputAS.get():
            outStructFileName = self._getPath('outputStructure.cif')
            # Write conservation in a section of the output cif file
            ASH = AtomicStructHandler()
            consScoresDic = self.mapConservation()
            inpAS = toCIF(self.inputAS.get().getFileName(), self._getTmpPath('inputStruct.cif'))
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, consScoresDic, self._ATTRNAME)
            ASH._writeLowLevel(outStructFileName, cifDic)

            AS = AtomStruct(filename=outStructFileName)
            self._defineOutputs(outputAtomStruct=AS)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        warns = []
        if not hasattr(self.inputSequences.get(), 'aligned') or not getattr(self.inputSequences.get(), 'aligned'):
            warns.append('Input sequences must be aligned to perform the conservation analysis.')
        return warns

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    #Main functions
    def calcConsensus(self):
        inFasta = self._getPath('alignedSequences.fasta')
        inSeqs = self.fillInputSequences()
        for seq in inSeqs:
            seq.appendToFile(inFasta)
        alignment = AlignIO.read(inFasta, 'fasta')
        summary_align = AlignInfo.SummaryInfo(alignment)
        return summary_align.gap_consensus(0.5)

    def calcConservation(self):
        seqsArr = self.getSequencesArray()
        if self.getEnumText('method') == SHANNON:
            values = self.calcShannon(seqsArr)
        elif self.getEnumText('method') == SIMPSON:
            values = self.calcSimpson(seqsArr)
        elif self.getEnumText('method') == KABAT:
            values = self.calcKabat(seqsArr)
        elif self.getEnumText('method') == PROP:
            values = self.calcProp(seqsArr)

        outFile = self.getConservationFile()
        with open(outFile, 'w') as f:
            for value in values:
                f.write('{}\t'.format(value))
        return outFile

    def getROIs(self):
        with open(self.getConservationFile()) as f:
            consValues = f.readline().split()

        flexThres = self.flexThres.get()
        if not flexThres:
            flexThres = self.thres.get()

        rois = []
        inRoi, fails = False, 0
        for i, v in enumerate(consValues):
            v = float(v)
            if (v > self.thres.get() and self.direction.get() == 1) or \
                    (v < self.thres.get() and self.direction.get() == 0):
                fails = 0
                if not inRoi:
                    inRoi = True
                    iniRoi = i + 1

            elif ((v > flexThres and self.direction.get() == 1) or \
                    (v < flexThres and self.direction.get() == 0)) and inRoi:
                if fails >= self.numFlex.get():
                    rois.append([iniRoi, i + 1])
                    inRoi = False
                else:
                    fails += 1

            elif inRoi:
                rois.append([iniRoi, i + 1])
                inRoi = False
        return rois

    def mapConservation(self):
        '''Return a dictionary with {spec: value}
        "spec" should be a chimera specifier. In this case:  chainId:residueIdx'''
        # Map positions of the original sequence with the actual index in the structure (might not start on 1)
        # {OriPos: Idx}
        structDic = {}
        parser = PDBParser()
        chain_id, modelId = json.loads(self.chain_name.get())["chain"], json.loads(self.chain_name.get())["model"]
        structModel = parser.get_structure(self.getASName(), self.getASFileName())[int(modelId)]
        for i, residue in enumerate(structModel[json.loads(self.chain_name.get())['chain']]):
          structDic[i + 1] = residue.get_id()[1]

        #Maps original positions of the sequence with respect to the new aligned
        # {OriPos: AlignPos}
        alignDic = {}
        with open(self._getExtraPath('structSequence.fasta.map')) as fIn:
            for i in range(2):
                fIn.readline()
            for line in fIn:
                oriPos, alignPos = int(line.split(',')[1].strip()), line.split(',')[2].strip()
                if alignPos != '-':
                    alignDic[oriPos] = int(alignPos)

        # Map positions of the alignment with conservation values
        # {AlignPos: ConsValue}
        consDic = {}
        with open(self.getConservationFile()) as fcons:
            line = fcons.readline()
            for i, value in enumerate(line.split()):
                consDic[i+1] = value

        print(consDic)

        #Final mapping: {"chain:Idx": ConsValue}
        mapDic = {}
        for oriPos in structDic:
            spec = '{}:{}'.format(chain_id, structDic[oriPos])
            if oriPos in alignDic:
                mapDic[spec] = float(consDic[alignDic[oriPos]])
            else:
                print('Original position {} of the structure sequence could not be mapped to alignment'.format(oriPos))
        return mapDic


    ##################

    def getInputSequence(self):
        inputObj = getattr(self, 'inputSequences')[0].get()
        seq = inputObj.getSequence()
        seq_name = inputObj.getSequenceObj().getId()
        if not seq_name:
          seq_name = inputObj.getSeqName()
        return seq, seq_name

    def getInputASSequence(self):
        inputObj = getattr(self, 'inputAS').get()
        seq_name = os.path.basename(inputObj.getFileName())
        handler = AtomicStructHandler(inputObj.getFileName())
        chainName = getattr(self, 'chain_name').get()

        # Parse chainName for model and chain selection
        struct = json.loads(chainName)  # From wizard dictionary
        chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

        seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))
        seqFile = self._getExtraPath('structSequence.fasta')
        with open(seqFile, "w") as f:
          f.write(('>{}\n{}\n'.format(seq_name, seq)))
        return seq, seq_name, seqFile

    def mapResidues(self, structModel):
        '''Returns a dictionary which maps the idxs of the residues of  the sequence and the sequence from a structure
        {idxSeq: idxStrSeq}'''
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

    def getASFileName(self):
        inpStruct = self.inputAS.get()
        inpFile = inpStruct.getFileName()
        if inpFile.endswith('.cif'):
            inpPDBFile = self._getExtraPath(os.path.basename(inpFile).replace('.cif', '.pdb'))
            cifToPdb(inpFile, inpPDBFile)

        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpPDBFile = self._getExtraPath(os.path.basename(inpFile).replace(inpStruct.getExtension(), '.pdb'))
            inpStruct.convert2PDB(outPDB=inpPDBFile)

        else:
            inpPDBFile = inpFile

        return inpPDBFile

    def getASName(self):
      return os.path.splitext(os.path.basename(self.getASFileName()))[0]

    def getMaxLenSeq(self, seqSet):
      maxi = 0
      for seq in seqSet:
          leni = len(seq.getSequence())
          if leni > maxi:
              maxi = leni
      return maxi

    def fillToLen(self, listi, leni, val='-'):
        if len(listi) < leni:
            if type(listi) == list:
                listi += [val] * (leni-len(listi))
            elif type(listi) == str:
                listi += val * (leni - len(listi))
        return listi

    def fillInputSequences(self):
        seqSet = []
        inSet = self.inputSequences.get()
        maxLen = self.getMaxLenSeq(inSet)
        for seq in inSet:
            seqSet.append(Sequence(sequence=self.fillToLen(seq.getSequence(), maxLen)))
        return seqSet

    def getElementNumber(self, listi):
        counts = []
        for elem in set(listi):
          counts.append(list(listi).count(elem))
        return counts

    def getSequencesArray(self):
        seqMat = []
        inSeqs = self.fillInputSequences()
        for seq in inSeqs:
            seqMat.append(np.array(list(seq.getSequence())))
        return np.array(seqMat)

    def getConservationFile(self):
        return self._getPath('conservationValues.tsv')

    def calcShannon(self, seqsArr, normalized=True):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            ent = entropy(counts, base=2)
            if normalized:
                ent = ent / math.log2(21)
            entrs.append(ent)
        return entrs

    def calcSimpson(self, seqsArr):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            entrs.append(simpson(counts))
        return entrs

    def calcKabat(self, seqsArr, normalized=True):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            ent = kabat(counts)
            if normalized:
                ent = (ent-1) / 440
            entrs.append(ent)
        return entrs

    def calcProp(self, seqsArr):
        '''Proportion of the more common residue (inverse: 1- prop)'''
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:,j])
            ent = 1 - max(counts) / sum(counts)
            entrs.append(ent)
        return entrs
