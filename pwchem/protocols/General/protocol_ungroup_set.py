# **************************************************************************
# *
# * Authors:     Your Name (your.email@example.com)
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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params


class ProtUngroupSet(EMProtocol):
    """Ungroup a Set into individual objects.

    Takes a SetOf[Something] (e.g. SetOfAtomStructs) and outputs each element
    as a separate individual object (e.g. AtomStruct).
    """
    _label = 'ungroup set'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass='EMSet',
                      label='Input set: ',
                      help='Set of objects to ungroup into individual items.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('defineOutputStep')

    # --------------------------- STEPS subfunctions --------------------------
    def defineOutputStep(self):
        inputSet = self.inputSet.get()

        for i, item in enumerate(inputSet):
            newItem = item.clone()

            # Give the cloned item a fresh, unique output path if it has one
            if hasattr(newItem, 'getSeqName') and newItem.getSeqName():
                newItem.setSeqName(item.getSeqName())
            elif hasattr(newItem, 'getFileName') and newItem.getFileName():
                newItem.setFileName(item.getFileName())


            outputName = f'output{i + 1}'
            self._defineOutputs(**{outputName: newItem})
            self._defineSourceRelation(inputSet, newItem)

    # --------------------------- INFO functions ------------------------------
    def _summary(self):
        summary = []
        if self.isFinished():
            inputSet = self.inputSet.get()
            summary.append(f'Ungrouped {inputSet.getSize()} items from the input set.')
        return summary

    def _methods(self):
        return []

    def _validate(self):
        validations = []
        inputSet = self.inputSet.get()
        if inputSet is not None and inputSet.getSize() == 0:
            validations.append('The input set is empty.')
        return validations