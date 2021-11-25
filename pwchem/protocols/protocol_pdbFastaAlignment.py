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
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

from pwem.convert import AtomicStructHandler
from pwem.objects import Sequence, AtomStruct, PdbFile
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import StringParam, EnumParam, FileParam, PointerParam, BooleanParam
from pwchem.objects import DatabaseID, SetOfDatabaseID, SequenceFasta, SequenceVariants


class ProtChemPdbFastaAlignment(EMProtocol):
    """Download a fasta from Uniprot and insert variants"""
    _label = 'pdb fasta alignment'

    def _defineParams(self, form):
        form.addSection(label='Input')


        form.addParam('inputOtherFile', PointerParam, pointerClass='Sequence, AtomStruct',
                      label='Input file:', allowsNull=False, help="File where a sequence is extracted")


        form.addParam('inputOtherFile_class', BooleanParam,
                      label='Chain Selected: ', default=True,
                      help="Chain from PDB selected")


        form.addParam('inputchain', StringParam, condition='inputOtherFile_class',
                      label='Chain:', allowsNull=False,
                      help="chain selected")

        form.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      label='Input Sequence:', allowsNull=False, help="Original Sequence")

    def _insertAllSteps(self):
        self._insertFunctionStep('alignmentFasta')

    def alignmentFasta(self):

        inputOtherFile_class = self.inputOtherFile.get()
        # print('inputOtherFile_class: ', inputOtherFile_class)

        inputSequence = self.inputSequence.get()

        if issubclass(type(inputOtherFile_class), AtomStruct):
            inputSequence.pairWiseAlignment(self.inputOtherFile.get(), chainID=self.inputchain.get())

        if not self.fromTypeID:
            inputSequence.pairWiseAlignment(self.inputOtherFile.get())

        aln_file = open(os.getcwd() + '/' + 'clustalw_pairWise.aln', "r")

        fn_clustal = self._getPath('clustalw_pairWise.aln')
        clustal_file = open(fn_clustal, "w")

        for line in aln_file:
            clustal_file.write(line)

        aln_file.close()
        clustal_file.close()

        self._defineOutputs(outputAlignment=inputSequence)

    def _validate(self):
        errors = []
        return errors