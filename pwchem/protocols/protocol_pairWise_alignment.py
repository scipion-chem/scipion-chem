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

# import lxml.etree as ET
import os
import sys
import urllib.request

import Bio
from Bio.Align.Applications import ClustalwCommandline, ClustalOmegaCommandline
from Bio import AlignIO

from pwem.convert import AtomicStructHandler
from pwem.objects import Sequence, AtomStruct, PdbFile
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import StringParam, EnumParam, FileParam, PointerParam, BooleanParam, TextParam, LabelParam
from pwchem.objects import DatabaseID, SetOfDatabaseID, SequenceFasta, SequenceVariants, SequencesAlignment
from pwem.convert.sequence import (alignClustalSequences, alignBioPairwise2Sequences, alignMuscleSequences)

from pyworkflow.object import (Object, Float, Integer, String,
                               OrderedDict, CsvList, Boolean, Set, Pointer,
                               Scalar)


class ProtChemPdbFastaAlignment(EMProtocol):
    """Download a fasta from Uniprot and insert variants"""
    _label = 'pdb fasta alignment'

    def _defineParams(self, form):
        form.addSection(label='Input')

        group = form.addGroup('Sequence 1')
        group.addParam('CondAtomStruct1', BooleanParam,
                       label='It is a PDB file? : ', default=True, help="PDB for sequence alignment")

        group.addParam('inputAtomStruct', PointerParam, pointerClass='AtomStruct', condition='CondAtomStruct1',
                       label='Input PBD file:', allowsNull=False, help="Atomic Structure where a sequence is extracted")

        group.addParam('chain_name', StringParam, condition='CondAtomStruct1',
                       default='', label='List of chains: ',
                       help="Chains in the PDB input file")

        group.addParam('inputOtherFileSequence1', PointerParam, pointerClass='Sequence', condition='not CondAtomStruct1',
                      label='Input Sequence:', allowsNull=False, help="Original Sequence")




        group = form.addGroup('Sequence 2')
        group.addParam('CondAtomStruct2', BooleanParam,
                      label='It is a PDB file? : ', default=True, help="PDB for sequence alignment")

        group.addParam('inputAtomStruct2', PointerParam, pointerClass='AtomStruct', condition='CondAtomStruct2',
                      label='Input PBD file:', allowsNull=False, help="Atomic Structure where a sequence is extracted")

        group.addParam('chain_name2', StringParam, condition='CondAtomStruct2',
                      default='', label='List of chains: ',
                      help="Chains in the PDB input file")

        group.addParam('inputOtherFileSequence2', PointerParam, pointerClass='Sequence', condition='not CondAtomStruct2',
                      label='Input Sequence file:', allowsNull=False, help="Sequence where a sequence is extracted")





    def _insertAllSteps(self):
        self._insertFunctionStep('alignmentFasta')

    def alignmentFasta(self):

        # Sequemce 1

        conditionTypeFile = self.CondAtomStruct1.get()

        if conditionTypeFile == True:
            inputOtherFile_class = self.inputAtomStruct.get()

        if conditionTypeFile == False:
            inputOtherFile_class = self.inputOtherFileSequence1.get()

        if type(inputOtherFile_class) == Pointer:
            inputOtherFile_class = inputOtherFile_class.get()

        if issubclass(type(inputOtherFile_class), Sequence):
            seq1 = inputOtherFile_class.getSequence()

        elif issubclass(type(inputOtherFile_class), AtomStruct):
            from pwem.convert.atom_struct import AtomicStructHandler
            handler = AtomicStructHandler(inputOtherFile_class.getFileName())
            chainName = self.chain_name.get()

            # Parse chainName for model and chain selection
            split_chainName = chainName.split(':')
            print('split_chainName', split_chainName)
            modelID = split_chainName[1]
            chainID = split_chainName[2]
            value_modelID = int(modelID.split(',')[0])
            prevalue_chainID = str(chainID.split(',')[0])
            split_prevalue_chainID = prevalue_chainID.split('"')
            value_chainID = split_prevalue_chainID[1]

            seq1 = str(handler.getSequenceFromChain(modelID=value_modelID, chainID=value_chainID))


        # Sequence 2

        conditionTypeFile = self.CondAtomStruct2.get()

        if conditionTypeFile == True:
            inputOtherFile_class = self.inputAtomStruct2.get()

        if conditionTypeFile == False:
            inputOtherFile_class = self.inputOtherFileSequence2.get()

        if type(inputOtherFile_class) == Pointer:
            inputOtherFile_class = inputOtherFile_class.get()

        if issubclass(type(inputOtherFile_class), Sequence):
            seq2 = inputOtherFile_class.getSequence()

        elif issubclass(type(inputOtherFile_class), AtomStruct):
            from pwem.convert.atom_struct import AtomicStructHandler
            handler = AtomicStructHandler(inputOtherFile_class.getFileName())
            chainName = self.chain_name2.get()

            # Parse chainName for model and chain selection
            split_chainName = chainName.split(':')
            print('split_chainName', split_chainName)
            modelID = split_chainName[1]
            chainID = split_chainName[2]
            value_modelID = int(modelID.split(',')[0])
            prevalue_chainID = str(chainID.split(',')[0])
            split_prevalue_chainID = prevalue_chainID.split('"')
            value_chainID = split_prevalue_chainID[1]

            seq2 = str(handler.getSequenceFromChain(modelID=value_modelID, chainID=value_chainID))

        # Fasta file with sequences
        pairWaise_fasta = 'pairWaise' + '.fasta'
        fn_pairWaise = open(pairWaise_fasta, "w")

        fn_pairWaise.write('>Sequence_1' + "\n")
        fn_pairWaise.write(seq1 + "\n")
        fn_pairWaise.write('>Sequence_2' + "\n")
        fn_pairWaise.write(seq2)

        fn_pairFile = self._getPath('pairWaise.fasta')
        print('fn_pairFile: ', fn_pairFile)
        fn_pairWaise = open(pairWaise_fasta, "r")
        for linePair in fn_pairWaise:
               print(linePair)

        # Alignment
        in_file = pairWaise_fasta
        out_file = self._getPath("clustalo_pairWise.aln")
        cline = alignClustalSequences(in_file, out_file)
        args = ' --outfmt=clu'
        self.runJob(cline, args)

        outputName = os.path.basename(out_file)
        out_fileClustal=SequencesAlignment(alignmentFileName=outputName)

        self._defineOutputs(outputVariants=out_fileClustal)


    def _validate(self):
        errors = []
        return errors