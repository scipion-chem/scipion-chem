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

from pwem.objects import Sequence
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import *

class ProtChemInsertVariants(EMProtocol):
    """Generates a variant / lineage or insert mutations in a given sequence"""
    _label = 'Generates a modified sequence'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSequenceVariants', PointerParam, pointerClass='SequenceVariants',
                      label='Input Sequence:', allowsNull=False,
                      help="Input sequence containing the information about the variants "
                           "and their corresponding mutations")

        form.addParam('fromVariant', BooleanParam, label='Generate predefined variant: ',
                      default='True',
                      help='Generates a variant (set of mutations)')
        form.addParam('selectMutation', StringParam,
                      label='Select a mutation: ', condition='not fromVariant',
                      help="Mutation to be inserted into the sequence")

        form.addParam('addMutation', LabelParam, condition='not fromVariant',
                      label='Add the mutation to the list: ',
                      help='Add the mutation to the list')

        form.addParam('toMutateList', TextParam, width=70, condition='not fromVariant',
                      default='', label='List of mutations: ',
                      help="Mutations selected to generate the new sequence")

        form.addParam('selectVariant', StringParam,
                      label='Select a variant: ', condition='fromVariant',
                      help="Variant to be generated")

    def _insertAllSteps(self):
        if not self.fromVariant:
            self._insertFunctionStep('insertMutationsStep')
        else:
            self._insertFunctionStep('generateVariantStep')

    def insertMutationsStep(self):
        fnSeq = self.inputSequenceVariants.get().getSequence()
        posSelected_input = self.toMutateList.get()
        listVariants = posSelected_input.rstrip().split('\n')
        aminoacid_original = []
        aminoacid_exchanged = []
        position_inSequence = []
        mutList = []
        for mutation in listVariants:
            mutation = mutation.split()[0]
            mutList.append(mutation)
            aminoacid_original.append(mutation[0])
            aminoacid_exchanged.append(mutation[-1])
            position_inSequence.append(mutation[1:-1])
        print('List of mutations: ', mutList)

        letters_fnSeq = list(fnSeq)
        for ele_position_inSequence in range(0, len(position_inSequence)):
            letters_fnSeq[int(position_inSequence[ele_position_inSequence]) - 1] = \
                aminoacid_exchanged[int(ele_position_inSequence)]
        mutatedSequence = ''.join(letters_fnSeq)
        id_sequence = self.inputSequenceVariants.get().getId()
        seqMutant = Sequence(sequence=mutatedSequence, name=(id_sequence + '_MutantsSequence'))
        seqMutant.setSequence(mutatedSequence)

        self._defineOutputs(outputSequence=seqMutant)

    def generateVariantStep(self):
        selectedVariant = self.selectVariant.get()
        inputSequenceVar = self.inputSequenceVariants.get()
        varDic = inputSequenceVar.getMutationsInLineage()
        print('List of mutations: ', varDic[selectedVariant])
        id_sequence = self.inputSequenceVariants.get().getId()

        variantSequence = inputSequenceVar.generateVariantLineage(selectedVariant)
        variantSequence = Sequence(sequence=variantSequence, name=(id_sequence + '_VariantSequence'),
                                   id=id_sequence + '_variant')
        self._defineOutputs(outputVariantSequence=variantSequence)

    def _validate(self):
        errors = []
        return errors


