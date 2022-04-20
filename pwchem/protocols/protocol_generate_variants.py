# **************************************************************************
# *
# * Authors:  Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
# *           Eugenia Ulzurrun (mariaeugenia.ulzurrun@cib.csic.es)
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

from pwem.objects import Sequence, SetOfSequences
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import *


class ProtChemGenerateVariants(EMProtocol):
    """Generates a set of variants from a sequence and a list of mutations / variants"""
    _label = 'Generates a set of sequence variants'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSequenceVariants', PointerParam, pointerClass='SequenceVariants',
                      label='Input Sequence Variants:', allowsNull=False,
                      help="Sequence containing the information about the variants and mutations")
        form.addParam('fromVariant', BooleanParam, label='Generate predefined variant: ',
                      default='True',
                      help='Generates a variant (set of mutations)')

        form.addParam('selectVariant', StringParam, condition='fromVariant',
                      label='Select a predefined variant:',
                      help="Variant to be generated")

        form.addParam('selectMutation', StringParam,
                      label='Select some mutations: ', condition='not fromVariant',
                      help="Mutations to be inserted into the sequence.\nYou can do multiple selection")

        form.addParam('addVariant', LabelParam,
                      label='Add variant or mutation:',
                      help='Add variant or mutation to the list')

        form.addParam('toMutateList', TextParam, width=70,
                      default='', label='List of variants or mutations:',
                      help='Generates a predefined list of variants or/and sequences with single mutations')


    def _insertAllSteps(self):
        self._insertFunctionStep('generateVariantStep')

    def generateVariantStep(self):
        posSelected_input = self.toMutateList.get()
        listVariants = posSelected_input.rstrip().split('\n')

        sequenceSet = SetOfSequences().create(outputPath=self._getPath())
        inputSequenceVar = self.inputSequenceVariants.get()
        for variantName in listVariants:
            mutsInfo = variantName.split(':')[1].strip()
            if 'Var:' in variantName:
                variantSequence = inputSequenceVar.generateVariantLineage(mutsInfo)
            elif 'Muts:' in variantName:
                variantSequence = inputSequenceVar.performSubstitutions(mutsInfo.split(','))

            seqVarObj = Sequence(sequence=variantSequence, name=variantName, id=variantName)
            sequenceSet.append(seqVarObj)

        sequenceSet.exportToFile(self._getPath('SetOfSequences.fasta'))
        self._defineOutputs(outputSetOfSequences=sequenceSet)


    def _validate(self):
        errors=[]
        return errors


