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
# import pwem.objects.data as data, EMFile

from pwem.objects import Sequence, EMFile, SetOfSequences
from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import *
from pyworkflow.protocol import params
from pwchem.objects import DatabaseID, SetOfDatabaseID, SequenceVariants

class ProtChemGenerateSequences(EMProtocol):
    """Generates a set of sequence from a sequence and a list of mutations / variants"""
    _label = 'Generates a set of sequences'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputSequenceVariants', PointerParam, pointerClass='SequenceVariants',
                      label='Input Sequence:', allowsNull=False,
                      help="Sequence containing the information about the variants and mutations")

        form.addParam('selectVariant', StringParam,
                      label='Select a predefined variant or mutation:',
                      help="Variant to be generated or mutation to be inserted")

        form.addParam('addVariant', LabelParam,
                      label='Add variant or mutation:',
                      help='Add variant or mutation')

        form.addParam('toMutateList', TextParam, width=70,
                      default='', label='List of variants or mutations:',
                      help='Generates a predefined variants or a customized sequences by a single mutation')


    def _insertAllSteps(self):
        self._insertFunctionStep('generateVariantStep')

    def generateVariantStep(self):
        posSelected_input = self.toMutateList.get()
        listVariants = posSelected_input.rstrip().split('\n')
        print('List of variants: ', listVariants)

        sequenceSet = SetOfSequences().create(outputPath=self._getPath())

        for variant in range(len(listVariants)):
            variantName = listVariants[variant]
            inputSequenceVar = self.inputSequenceVariants.get()
            variantSequence = inputSequenceVar.generateVariantLineage(variantName)

            seqVarObj = Sequence(sequence=variantSequence, name=variantName, id=variantName)
            sequenceSet.append(seqVarObj)

        sequenceSet.exportToFile(self._getPath('SetOfSequences.fasta'))

        self._defineOutputs(outputSetOfSequences=sequenceSet)


    def _validate(self):
        errors=[]
        return errors


