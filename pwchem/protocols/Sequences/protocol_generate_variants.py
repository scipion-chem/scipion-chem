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
    """
    This protocol is used to generate a set of sequence variants from a SequenceVariants object.

    The protocol allows the creation of new sequences derived from an input reference sequence
    by applying either predefined variant definitions or explicit mutation lists. Each generated
    variant is stored as a Sequence object within a SetOfSequences output.

    Overview
    --------
    The protocol supports two modes of sequence generation:
    - Variant-based generation: applies predefined mutation sets (lineages)
    - Mutation-based generation: applies user-specified residue substitutions

    Input
    -----
    inputSequenceVariants:
        SequenceVariants object containing the reference sequence and mutation definitions.

    fromVariant:
        Defines the generation mode:
        - Variant: generate sequences from predefined variant definitions
        - Mutations: generate sequences from explicit mutation lists

    selectVariant:
        Name of a predefined variant to generate (used when fromVariant = Variant)

    selectMutation:
        List of mutations to apply to the sequence (used when fromVariant = Mutations)

    toMutateList:
        Internal list of selected variants/mutations to be processed

    Workflow
    --------
    1. Parse user-defined variant or mutation list
    2. Identify generation mode for each entry
    3. Generate sequence:
       - Original sequence (no modification)
       - Variant lineage reconstruction
       - Substitution-based mutant generation
    4. Create Sequence objects for each result
    5. Store all sequences in a SetOfSequences container
    6. Export resulting sequences to FASTA format

    Output
    ------
    outputSequences:
        SetOfSequences containing:
        - Reference sequence (if selected)
        - Variant-derived sequences
        - Mutation-derived sequences

    Key Features
    ------------
    - Supports both variant lineage reconstruction and explicit mutation application
    - Produces standardized Sequence objects
    - Exports results in FASTA format for interoperability
    - Integrates with SequenceVariants mutation models

    Notes
    -----
    - Variants are assumed to be defined in the input SequenceVariants object
    - Mutations are applied as residue substitutions on the reference sequence
    - Sequence identifiers encode origin and mutation information for traceability
    """
    _label = 'Generates a set of sequence variants'
    _originOptions = ['Variant', 'Mutations']

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSequenceVariants', PointerParam, pointerClass='SequenceVariants',
                       label='Input Sequence Variants:', allowsNull=False,
                       help="Sequence containing the information about the variants and mutations")

        group = form.addGroup('Add variant or mutation')
        group.addParam('fromVariant', EnumParam, label='Generate sequence from: ',
                       default=0, choices=self._originOptions, display=EnumParam.DISPLAY_HLIST,
                       help='Generates a variant (set of mutations)')

        group.addParam('selectVariant', StringParam, condition='fromVariant==0',
                       label='Select a predefined variant:',
                       help="Variant to be generated")

        group.addParam('selectMutation', StringParam,
                       label='Select some mutations: ', condition='fromVariant==1',
                       help="Mutations to be inserted into the sequence.\nYou can do multiple selection")

        group.addParam('addVariant', LabelParam,
                       label='Add variant or mutation:',
                       help='Add variant or mutation to the list')

        group = form.addGroup('Summary')
        group.addParam('toMutateList', TextParam, width=70,
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
            if 'Original' == mutsInfo.strip():
                variantSequence = inputSequenceVar.getSequence()
                seqId = inputSequenceVar.getId()

            #Variant origin
            elif '{}:'.format(self._originOptions[0]) in variantName:
                variantSequence = inputSequenceVar.generateVariantLineage(mutsInfo)
                seqId = '{}_{}'.format(inputSequenceVar.getId(), mutsInfo)

            #Mutations origin
            elif '{}:'.format(self._originOptions[1]) in variantName:
                variantSequence = inputSequenceVar.performSubstitutions(mutsInfo.split(','))
                seqId = '{}_mutant{}'.format(inputSequenceVar.getId(), variantName.split(')')[0])

            seqVarObj = Sequence(sequence=variantSequence, name=seqId, id=seqId)
            sequenceSet.append(seqVarObj)

        sequenceSet.exportToFile(self._getPath('SetOfSequences.fasta'))
        self._defineOutputs(outputSequences=sequenceSet)


    def _validate(self):
        errors=[]
        return errors

    # ADD WIZARD UTILS
    def getOriginLabel(self):
      if self.fromVariant.get() == 0:
        inputLabel, sumLabel = 'selectVariant', 'Variant'
      elif self.fromVariant.get() == 1:
        inputLabel, sumLabel = 'selectMutation', 'Mutations'
      return inputLabel, sumLabel

    def buildSumLine(self, type=None):
      inputLabel, sumLabel = self.getOriginLabel()
      roiInfo = getattr(self, inputLabel).get()
      return f'{sumLabel}: {roiInfo}'



