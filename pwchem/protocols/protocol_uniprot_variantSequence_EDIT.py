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

from pwem.objects import Sequence
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import *
from pwchem.objects import DatabaseID, SetOfDatabaseID, SequenceFasta, SequenceVariants

class InsertVariants(EMProtocol):
    """Download a fasta from Uniprot and insert variants"""
    _label = 'insert variants in a sequence'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSequence', PointerParam, pointerClass='Sequence',
                      label='Input Sequence:', allowsNull=False, help="Original Sequence")

        form.addParam('selectVariant', StringParam,
                      label='Select a variant: ',
                      help="Variants to be insert into a sequence")

        form.addParam('addVariant', LabelParam,
                      label='Add variant to list of variants: ')

        form.addParam('toMutateList', TextParam, width=70,
                      default='', label='List of variants: ',
                      help="Variants to be applied into the sequence")

    def _insertAllSteps(self):
        self._insertFunctionStep('insertVariants')

    def insertVariants(self):
        # llamar al objeto SequenceVariants

        fnVars = self.inputSequence.get().getVariantsFileName()
        print('fnVars: ' , fnVars)
        fnSeq = self.inputSequence.get().getSequence()
        posSelected_input = self.toMutateList.get()
        print('posSelected_input: ', posSelected_input)

        middle_variants = posSelected_input.rstrip().split('\n')
        print('middle_variants: ', middle_variants)

        fileVariants  = middle_variants
        print('fileVariants: ', fileVariants)

        aminoacid_original = []
        aminoacid_exchanged = []
        position_inSequence = []
        for mutation in fileVariants:
            mutation = mutation.split()[0]
            aminoacid_original.append(mutation[0])
            aminoacid_exchanged.append(mutation[-1])
            position_inSequence.append(mutation[1:-1])

        # print('aminoacid_original: ', aminoacid_original)
        # print('aminoacid_exchanged: ', aminoacid_exchanged)
        # print('position_inSequence: ', position_inSequence)
        letters_fnSeq = list(fnSeq)
        # print('letters_fnSeq: ', letters_fnSeq)
        for ele_position_inSequence in range(0, len(position_inSequence)):
            letters_fnSeq[int(position_inSequence[ele_position_inSequence]) - 1] = aminoacid_exchanged[int(ele_position_inSequence)]

        mutatedSequence = ''.join(letters_fnSeq)
        print(mutatedSequence)

        seqMutant = Sequence()
        seqMutant.setSequence(mutatedSequence)

        lineage = SequenceVariants()
        lineage.getMutationsInLineage(fnVars)

        self._defineOutputs(outputSequence=seqMutant)

    def _validate(self):
        errors=[]
        return errors


