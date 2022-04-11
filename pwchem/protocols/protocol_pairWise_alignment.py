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

from pyworkflow.object import Pointer
from pwem.objects import Sequence, AtomStruct
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import StringParam, PointerParam, BooleanParam
from pwem.convert.sequence import alignClustalSequences

from pwchem.objects import SequencesAlignment
from pwchem.utils.utilsFasta import pairwiseAlign


class ProtChemPairWiseAlignment(EMProtocol):
    """Pairwise alignment using Clustal Omega"""
    _label = 'pairwise alignment'

    def _defineParams(self, form):
        form.addSection(label='Input')

        group = form.addGroup('Sequence 1')
        group.addParam('condAtomStruct1', BooleanParam,
                       label='Is it a PDB file?: ', default=True, help="PDB for alignment")
        group.addParam('inputAtomStruct1', PointerParam, pointerClass='AtomStruct', condition='condAtomStruct1',
                       label='Input first PBD file:',
                       help="Atomic structure where a sequence is extracted")
        group.addParam('chain_name1', StringParam, condition='condAtomStruct1',
                       default='', label='PDB chains:',
                       help="Chains in the PDB")
        group.addParam('inputSequence1', PointerParam, pointerClass='Sequence',
                       condition='not condAtomStruct1',
                       label='Input first sequence:', help="Sequence for alignment")


        group = form.addGroup('Sequence 2')
        group.addParam('condAtomStruct2', BooleanParam,
                       label='Is it a PDB file?:', default=True, help="PDB for sequence alignment")
        group.addParam('inputAtomStruct2', PointerParam, pointerClass='AtomStruct', condition='condAtomStruct2',
                       label='Input second PBD file:',
                       help="Atomic Structure where a sequence is extracted")
        group.addParam('chain_name2', StringParam, condition='condAtomStruct2',
                       default='', label='List of chains:',
                       help="Chains in the PDB input file")
        group.addParam('inputSequence2', PointerParam, pointerClass='Sequence',
                       condition='not condAtomStruct2',
                       label='Input second sequence:', help="Sequence for alignment")


    def _insertAllSteps(self):
        self._insertFunctionStep('alignmentFasta')

    def alignmentFasta(self):
        # Fasta file with sequences
        seq1, seq_name1 = self.getInputSequence(1)
        seq2, seq_name2 = self.getInputSequence(2)

        out_file = os.path.abspath(self._getPath("clustalo_pairWise.aln"))
        pairwiseAlign(seq1, seq2, out_file, seqName1=seq_name1, seqName2=seq_name2)

        out_fileClustal=SequencesAlignment(alignmentFileName=out_file)
        self._defineOutputs(outputVariants=out_fileClustal)

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
