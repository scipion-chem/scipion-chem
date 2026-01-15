# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params



class ProtChemGroupByAtt(EMProtocol):
    """Group small molecules by attribute."""
    _label = 'group by attribute'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass="EMSet", label='Input set: ',
                      help='Set of small molecules to group by pocket id.')
        form.addParam('refColumn', params.StringParam, label='Grouping column: ', default='',
                       help='Attribute for the grouping operation.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    # --------------------------- STEPS subfunctions ------------------------------
    def defineOutputStep(self):
        inputSet = self.inputSet.get()
        attr = self.refColumn.get()

        groupingValues = self.getGroupingVals()

        for j, i in enumerate(groupingValues):
            print(f' j: {j} ; i: {i}')
            outputPath = self._getPath(f'group{j}')
            os.makedirs(outputPath, exist_ok=True)

            outputSet = inputSet.getClass()().create(outputPath=outputPath)

            for item in inputSet:
                if self.getAttrValue(item, attr) == i:
                    newItem = item.clone()  # generic clone for any item type
                    outputSet.append(newItem)

                # Optional: preserve docked status if attribute exists
            if hasattr(inputSet, 'isDocked') and hasattr(outputSet, 'setDocked'):
                outputSet.setDocked(inputSet.isDocked())
            if hasattr(inputSet, 'proteinFile') and hasattr(outputSet, 'proteinFile'):
                outputSet.proteinFile.set(inputSet.proteinFile.get())

            outputName = f'outputSet_{j}'
            self._defineOutputs(**{outputName: outputSet})
            self._defineSourceRelation(inputSet, outputSet)



    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputSet.get()
        if molSet.getClass() == 'SetOfSmallMols' and not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    # --------------------------- UTILS functions -----------------------------------
    def getAttrValue(self, item, opAttr):
        try:
            opId = getattr(item, opAttr).get()
        except:
            opId = getattr(item, opAttr)
        return opId

    def getGroupingVals(self):
        attr = self.refColumn.get()
        values = set()

        for mol in self.inputSet.get():
            values.add(self.getAttrValue(mol, attr))

        return sorted(values)
