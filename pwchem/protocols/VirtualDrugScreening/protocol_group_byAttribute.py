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
    """
    AI Generated:

    This protocol is used to group elements of a Scipion EMSet based on the value
    of a selected attribute. It splits an input dataset into multiple subsets,
    where each subset contains items sharing the same value for a given field.

    The protocol is generic and can be applied to different types of EMSets
    (e.g., small molecules, particles, features), as long as they contain a
    consistent attribute used for grouping.

    Core Concepts
    -------------
    Input Set:
        Any Scipion EMSet containing objects with accessible attributes.

    Grouping Attribute:
        A user-defined field (column) of the input objects used to partition the
        dataset into groups (e.g., pocketId, chainId, score, label).

    Grouped Output Sets:
        Multiple EMSets generated from the original one, each corresponding to a
        unique value of the grouping attribute.

    Attribute Access:
        Attributes may be stored as direct values or wrapped objects (e.g., Scipion
        attribute containers), requiring flexible extraction.

    Workflow
    --------
    1. Input EMSet is provided by the user.
    2. User selects an attribute to group by (refColumn).
    3. All unique values of the selected attribute are extracted.
    4. For each unique value:
       a. A new output directory is created.
       b. A new empty EMSet is initialized.
       c. Input items are iterated and filtered by attribute value.
       d. Matching items are cloned into the corresponding group set.
    5. Each grouped EMSet is saved as an independent output.
    6. Source relationships between input and outputs are recorded.

    Grouping Logic
    --------------
    - The grouping is based on exact equality of attribute values.
    - Values are collected into a sorted unique list before processing.
    - Each group is independent and preserves original item metadata.

    Attribute Handling
    ------------------
    The protocol supports two attribute access patterns:
    - Direct attribute values (e.g., int, str, float)
    - Wrapped Scipion objects requiring `.get()` extraction

    Output
    ------
    - outputSet0, outputSet1, ...:
        Multiple EMSet objects, one per unique attribute value.

    Each output set contains:
        - Cloned items from the input set
        - Original metadata preserved via copyInfo()

    Use Cases
    ---------
    - Grouping ligands by pocket or binding site
    - Splitting datasets by classification labels
    - Organizing results by score thresholds or categories
    - Preparing subsets for downstream analysis or visualization
    """
    _label = 'group by attribute'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass="EMSet", label='Input set: ',
                      help='Set of small molecules to group by pocket id.')
        form.addParam('refColumn', params.StringParam, label='Grouping column: ', default='',
                       help='Attribute for the grouping operation.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
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
                    newItem = item.clone()
                    outputSet.append(newItem)


            outputSet.copyInfo(inputSet)

            outputName = f'outputSet{j}'
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
