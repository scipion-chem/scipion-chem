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
This protocol is used to calculate conservation over a set of sequences and store it on the object SequenceChem

"""
import os, math, json
import numpy as np
from scipy.stats import entropy

from Bio import AlignIO
from Bio.Align import AlignInfo
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.objects import AtomStruct
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toCIF, AtomicStructHandler, addScipionAttribute, cifToPdb

from pwchem import Plugin as pwchemPlugin
from pwchem.objects import SetOfSequenceROIs, Sequence, SequenceChem
from pwchem.constants import BIOCONDA_DIC

SHANNON, SIMPSON, KABAT, PROP = 'Shannon Entropy', 'Simpson Diversity Index', 'Wu-kabat Variability coefficient', \
                                'Maximum Proportion Conservation'

def simpson(counts):
  sumi = 0
  N = sum(counts)
  for c in counts:
    sumi += c * (c - 1)
  return 1 - (sumi / (N * (N - 1)))


def kabat(counts):
  return (sum(counts) * len(counts)) / max(counts)

class ProtSeqCalculateConservation(EMProtocol):
    """
    Calculates conservation over a set of sequence and adds it into the object
    """
    _label = 'Calculate conservation'
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

        group.addParam('mapCons', params.EnumParam, default=0, label='Map sequence conservation to: ',
                       choices=['None', 'Structure', 'Sequence ROIs'],
                       help='Map the conservation from the sequence alignment into a structure or sequence ROIs '
                            'whose sequence must align to the sequence alignment too')
        group.addParam('useCons', params.BooleanParam, default=True, condition='mapCons==0',
                       label='Use consensus sequence for output: ',
                       help='Use consensus sequence for the sequence ROIs output or a specific input sequence')
        group.addParam('outSeq', params.StringParam, default='', condition='mapCons==0 and not useCons',
                       label='Output sequence: ',
                       help='Input sequence to use for the output sequence ROIs')

        group.addParam('inputAS', params.PointerParam, pointerClass='AtomStruct', condition='mapCons==1',
                       allowsNull=True, label="Input protein structure: ",
                       help='Input protein structure. If included, the conservation values will be drawed over'
                            'the surface of  the protein in the results analysis.')
        group.addParam('chain_name', params.StringParam, condition='mapCons==1', label='Chain of interest:',
                       help='Specify the chain of interest')
        
        group.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs', condition='mapCons==2',
                       allowsNull=True, label="Input sequence ROIs: ",
                       help='Map the conservation over the sequence ROIs')

        group = form.addGroup('Variability measure')
        group.addParam('method', params.EnumParam, default=3, label='Method to measure variability: ',
                       choices=[SHANNON, SIMPSON, KABAT, PROP],
                       help='Method to measure the conservation / variability in each position of the sequences.\n'
                            'http://imed.med.ucm.es/PVS/pvs-help.html#vmth')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('calculateConservationStep')
        if self.mapCons.get() != 0:
            self._insertFunctionStep('mapConservationStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateConservationStep(self):
        outFile = self.calcConservation()

    def mapConservationStep(self):
        oriSeqsFile = self._getExtraPath('originalAlignment.fasta')
        addSeqsFile = self._getExtraPath('structureAddAlignment.fasta')
        self.inputSequences.get().exportToFile(oriSeqsFile)
        
        _, _, seqFile = self.getMapSequence()
        cline0 = f'{pwchemPlugin.getEnvActivationCommand(BIOCONDA_DIC)} && '

        cline = cline0 + f'mafft --add {seqFile} --mapout {oriSeqsFile} > {addSeqsFile}'
        self.runJob(cline, '')
        if os.path.getsize(addSeqsFile) == 0:
            auxAlignFile = self._getExtraPath('mafftAlignment.fasta')
            print(f'Alignment of the structure sequence failed, most likely because mafft detected that input sequences'
                  f' were not aligned.\nAligning sequences with mafft in {auxAlignFile}')
            cline = cline0 + f'mafft --auto {oriSeqsFile} > {auxAlignFile}'
            self.runJob(cline, '')

            cline = cline0 + f'mafft --add {seqFile} --mapout {auxAlignFile} > {addSeqsFile}'
            self.runJob(cline, '')

    def defineOutputStep(self):
        outSeqObj = self.getOutputSequence()

        consDic = self.getConsDic()
        outSeq = SequenceChem(attributesFile=self._getExtraPath('sequenceAttributes.txt'))
        outSeq.copy(outSeqObj)
        outSeq.addAttributes({self.getEnumText('method'): list(consDic.values())})

        self._defineOutputs(outputSequence=outSeq)

        if self.mapCons.get() == 1:
            outStructFileName = self._getPath('outputStructure.cif')
            # Write conservation in a section of the output cif file
            ASH = AtomicStructHandler()
            consScoresDic = self.mapStructConservation()
            inpAS = toCIF(self.inputAS.get().getFileName(), self._getTmpPath('inputStruct.cif'))
            cifDic = ASH.readLowLevel(inpAS)
            cifDic = addScipionAttribute(cifDic, consScoresDic, self._ATTRNAME)
            ASH._writeLowLevel(outStructFileName, cifDic)

            AS = AtomStruct(filename=outStructFileName)
            self._defineOutputs(outputAtomStruct=AS)

        elif self.mapCons.get() == 2:
            alignDic = self.getAlignMapping()

            newROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
            for roi in self.inputROIs.get():
                idx0, idx1 = roi.getROIIdxs()
                consVal = self.getMeanConservation(consDic, alignDic[idx0], alignDic[idx1-1])
                setattr(roi, self.getEnumText('method'), params.Float(consVal))
                newROIs.append(roi)

            self._defineOutputs(outputROIs=newROIs)

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
            seq.appendToFile(inFasta, doClean=False)
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

    def getAlignMapping(self):
        '''Returns a dictionary of the mapping of the original indexes to the ones aligned to the MSA'''
        alignDic = {}
        with open(self._getExtraPath('mapSequence.fasta.map')) as fIn:
            for _ in range(2):
                fIn.readline()
            for line in fIn:
                oriPos, alignPos = line.split(',')[1].strip(), line.split(',')[2].strip()
                if alignPos != '-':
                    alignDic[int(oriPos)] = int(alignPos)
        return alignDic

    def mapStructConservation(self):
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
        alignDic = self.getAlignMapping()

        # Map positions of the alignment with conservation values
        # {AlignPos: ConsValue}
        consDic = self.getConsDic()

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

    def getMeanConservation(self, consDic, idx0, idx1):
        c = 0
        for i in range(idx0, idx1 + 1):
            c += float(consDic[i])
        return c / (idx1-idx0+1)

    def getConsDic(self):
        '''Return a dictionary with the conservation values over the sequence'''
        consDic = {}
        with open(self.getConservationFile()) as fcons:
            line = fcons.readline()
            for i, value in enumerate(line.split()):
                consDic[i + 1] = value

        return consDic

    def getOutputSequence(self):
        if self.mapCons.get() in [1, 2]:
            # Map to structure or ROI
            outSeq, _, _ = self.getMapSequence()

        else:
            if not self.useCons:
                outSeqStr = self.outSeq.get()
                for outSeq in self.inputSequences.get():
                    if outSeq.__str__() == outSeqStr:
                        break
            else:
                outStr = str(self.calcConsensus())
                outSeq = self.inputSequences.get().getFirstItem()
                outSeq.setSequence(outStr)
                outSeq.setSeqName('ConsensusSequence')

        return outSeq

    def getMapSequence(self):
        if self.mapCons.get() == 1:
            seqObj, seqName, seqFile = self.getInputASSequence()
        elif self.mapCons.get() == 2:
            seqObj, seqName, seqFile = self.getInputROIsSequence()
        return seqObj, seqName, seqFile

    def getMapSequenceFile(self):
        return self._getExtraPath('mapSequence.fasta')

    def getInputASSequence(self):
        inputObj = getattr(self, 'inputAS').get()
        seq_name = os.path.basename(inputObj.getFileName())
        handler = AtomicStructHandler(inputObj.getFileName())
        chainName = getattr(self, 'chain_name').get()

        # Parse chainName for model and chain selection
        struct = json.loads(chainName)  # From wizard dictionary
        chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

        seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))
        seqFile = self.getMapSequenceFile()
        with open(seqFile, "w") as f:
          f.write(('>{}\n{}\n'.format(seq_name, seq)))

        seqObj = self.inputSequences.get().getFirstItem()
        seqObj.setSequence(seq)
        seqObj.setSeqName('StructureSequence')
        return seqObj, seq_name, seqFile

    def getInputROIsSequence(self):
        seqObj = self.inputROIs.get().getSequenceObj()
        seqName = seqObj.getSeqName()
        seqObj.exportToFile(self.getMapSequenceFile())
        return seqObj, seqName, self.getMapSequenceFile()

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
            counts = self.getElementNumber(seqsArr[:, j])
            ent = entropy(counts, base=2)
            if normalized:
                ent = ent / math.log2(21)
            entrs.append(1 - ent)
        return entrs

    def calcSimpson(self, seqsArr):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:, j])
            entrs.append(1 - simpson(counts))
        return entrs

    def calcKabat(self, seqsArr, normalized=True):
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:, j])
            ent = kabat(counts)
            if normalized:
                ent = (ent-1) / 440
            entrs.append(1 - ent)
        return entrs

    def calcProp(self, seqsArr):
        '''Proportion of the more common residue'''
        entrs = []
        for j in range(seqsArr.shape[1]):
            counts = self.getElementNumber(seqsArr[:, j])
            ent = max(counts) / sum(counts)
            entrs.append(ent)
        return entrs
