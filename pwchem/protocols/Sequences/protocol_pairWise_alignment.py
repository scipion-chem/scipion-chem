# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, json

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import StringParam, PointerParam, BooleanParam
from pwchem.objects import SetOfSequencesChem, Sequence
from pwchem.utils.utilsFasta import pairwiseAlign, parseAlnFile


class ProtChemPairWiseAlignment(EMProtocol):
    """Pairwise alignment using Clustal Omega"""
    _label = 'pairwise alignment'

    def _defineParams(self, form):
        form.addSection(label='Input')

        group = form.addGroup('Sequence 1')
        group.addParam('condAtomStruct1', BooleanParam,
                       label='From atomic structure: ', default=True,
                       help="Whether the first input for the alignment is an atomic structure or a sequence")
        group.addParam('inputAtomStruct1', PointerParam, pointerClass='AtomStruct', condition='condAtomStruct1',
                       label='Input first atomic structure: ',
                       help="Atomic structure whose sequence will be align")
        group.addParam('chain_name1', StringParam, condition='condAtomStruct1',
                       default='', label='Structure chain:',
                       help="Chain of the structure to align")
        group.addParam('inputSequence1', PointerParam, pointerClass='Sequence',
                       condition='not condAtomStruct1',
                       label='Input first sequence:', help="First sequence for alignment")


        group = form.addGroup('Sequence 2')
        group.addParam('condAtomStruct2', BooleanParam,
                       label='From atomic structure: ', default=True,
                       help="Whether the first input for the alignment is an atomic structure or a sequence")
        group.addParam('inputAtomStruct2', PointerParam, pointerClass='AtomStruct', condition='condAtomStruct2',
                       label='Input second atomic structure: ',
                       help="Atomic structure whose sequence will be align")
        group.addParam('chain_name2', StringParam, condition='condAtomStruct2',
                       default='', label='Structure chain:',
                       help="Chain of the structure to align")
        group.addParam('inputSequence2', PointerParam, pointerClass='Sequence',
                       condition='not condAtomStruct2',
                       label='Input second sequence:', help="Second sequence for alignment")


    def _insertAllSteps(self):
        self._insertFunctionStep('alignmentFasta')

    def alignmentFasta(self):
        # Fasta file with sequences
        seq1, seq_name1 = self.getInputSequence(1)
        seq2, seq_name2 = self.getInputSequence(2)

        out_file = os.path.abspath(self._getPath("clustalo_pairwise.aln"))
        pairwiseAlign(seq1, seq2, out_file, seqName1=seq_name1, seqName2=seq_name2)
        outSeqDic = parseAlnFile(out_file)

        outSeqs = SetOfSequencesChem.create(outputPath=self._getPath())
        for seqId in outSeqDic:
            outSeqs.append(Sequence(name=seqId, sequence=outSeqDic[seqId], id=seqId))

        outSeqs.setAlignmentFileName(os.path.relpath(out_file))
        outSeqs.setAligned(True)
        self._defineOutputs(outputSequences=outSeqs)

    def _validate(self):
        errors = []
        return errors

    def getInputSequence(self, idx):
        condAS = getattr(self, 'condAtomStruct{}'.format(idx))
        if condAS:
            from pwem.convert.atom_struct import AtomicStructHandler
            inputObj = getattr(self, 'inputAtomStruct{}'.format(idx)).get()
            seq_name = os.path.basename(inputObj.getFileName())
            handler = AtomicStructHandler(inputObj.getFileName())
            chainName = getattr(self, 'chain_name{}'.format(idx))

            # Parse chainName for model and chain selection
            struct = json.loads(chainName.get())  # From wizard dictionary
            chain_id, modelId = struct["chain"].upper().strip(), int(struct["model"])

            seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chain_id))

        else:
            inputObj = getattr(self, 'inputSequence{}'.format(idx)).get()
            seq = inputObj.getSequence()
            seq_name = inputObj.getId()
            if not seq_name:
                seq_name = inputObj.getSeqName()
        return seq, seq_name

    def parseSequences(self, alignFile):
        pass