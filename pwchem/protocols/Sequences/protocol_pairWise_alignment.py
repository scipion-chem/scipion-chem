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
    """
    This protocol performs a pairwise sequence alignment between two biological sequences.

    The sequences can be provided either directly as amino acid sequences or
    extracted from atomic structures (e.g., PDB models). In the case of structures,
    a specific chain and model can be selected for alignment.

    The alignment is computed using a pairwise alignment method (Clustal Omega-based
    implementation), and the result is stored as a standard alignment file.

    The protocol then parses the alignment output and converts it into a set of
    sequence objects, allowing further analysis within the Scipion framework.

    Workflow
    --------
    1. Input selection of two sequences:
       - Either atomic structures (with chain/model selection)
       - Or raw biological sequences

    2. Sequence extraction:
       - If structure is provided, the amino acid sequence is extracted from the selected chain
       - Otherwise, the input sequence is used directly

    3. Pairwise alignment:
       - The two sequences are aligned using an external alignment tool
       - The resulting alignment is saved to disk

    4. Output processing:
       - The alignment file is parsed
       - Aligned sequences are reconstructed
       - A SetOfSequencesChem object is created with the aligned sequences

    5. Output generation:
       - A collection of aligned sequences is produced
       - The alignment file is attached for reproducibility

    Key Concepts
    ------------
    Atomic Structure Input:
        A protein structure (e.g., PDB) from which a specific chain is selected
        and converted into an amino acid sequence.

    Sequence Input:
        A raw amino acid sequence provided directly by the user.

    Pairwise Alignment:
        A method that aligns two sequences to identify conserved regions,
        mismatches, and gaps.

    Alignment Output:
        A text-based representation of aligned sequences with gaps inserted
        to maximize similarity.

    Output
    ------
    outputSequences:
        A SetOfSequencesChem object containing the two aligned sequences.

        Each sequence is stored with:
        - Identifier
        - Aligned amino acid sequence
        - Original metadata (name, source)

    The alignment file is also stored for downstream analysis.

    Use Cases
    ---------
    - Structural vs sequence comparison
    - Homology analysis
    - Mutation and conservation analysis
    - Protein sequence validation
    """
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
            seq_name = inputObj.getSeqName()
        return seq, seq_name

    def parseSequences(self, alignFile):
        pass