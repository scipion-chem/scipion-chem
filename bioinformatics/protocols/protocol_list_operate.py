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

from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer
from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtBioinformaticsListOperate(EMProtocol):
    """Filter a set by a column value or keep just a few columns"""
    _label = 'operate set'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('operation', EnumParam, choices=['Filter', 'Keep columns'],
                      label='Operation', default=0)
        form.addParam('inputSet', PointerParam, pointerClass="EMSet",
                       label='Set to filter:', allowsNull=False)
        form.addParam('filterColumn', StringParam,
                       label='Filter column:', condition='(operation==0)',
                       help='It must exist in the input database list.')
        form.addParam('filterOp', EnumParam, choices=['==', '>', '>=', '<', '<=', '!=', 'startswith',
                                                      'endswith', 'contains', 'does not startwith',
                                                      'does not end with', 'does not contain'],
                       label='Filter operation:', condition='(operation==0)')
        form.addParam('filterValue', StringParam,
                       label='Value:', condition='(operation==0)',
                       help='Value to use in the filter')
        form.addParam('keepColumns', StringParam,
                       label='Keep columns:', condition='(operation==1)',
                       help='They must exist in the input database list. Separated by semicolons')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        outputSet = self.inputSet.get().create(self._getPath())
        if self.operation.get()==0:
            # Filter columns
            referenceValue = self.filterValue.get()
            value = self.inputSet.get().getFirstItem().getAttributeValue(self.filterColumn.get())
            if isinstance(value,float):
                referenceValue = float(referenceValue)
            elif isinstance(value,int):
                referenceValue = int(referenceValue)

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if isinstance(value, Float):
                    value = float(value)
                elif isinstance(value, Integer):
                    value = int(value)

                filterOp = self.filterOp.get()
                add = False
                if filterOp == 0: # ==
                    add = value==referenceValue
                elif filterOp == 1: # >
                    add = value>referenceValue
                elif filterOp == 2:  # >=
                    add = value > referenceValue
                elif filterOp == 3:  # <
                    add = value < referenceValue
                elif filterOp == 4:  # <=
                    add = value <= referenceValue
                elif filterOp == 5:  # !=
                    add = value != referenceValue
                elif filterOp == 6:  #startswith
                    add = value.startswith(referenceValue)
                elif filterOp == 7:  # endswith
                    add = value.endswith(referenceValue)
                elif filterOp == 8:  # contains
                    add = referenceValue in value
                elif filterOp == 9:  # does not startswith
                    add = not (value.startswith(referenceValue))
                elif filterOp == 10:  # does not endswith
                    add = not (value.endswith(referenceValue))
                elif filterOp == 11:  # does not contains
                    add = not (referenceValue in value)
                if add:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==1:
            # Keep columns
            keepList=[x.strip() for x in self.keepColumns.get().split()]

            ignoreList=[]
            for name, _ in self.inputSet.get().getFirstItem().getAttributes():
                if not name in keepList:
                    ignoreList.append(name)
            for oldEntry in self.inputSet.get():
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(oldEntry,ignoreAttrs=ignoreList)
                outputSet.append(newEntry)

        self._defineOutputs(output=outputSet)
        self._defineSourceRelation(self.inputSet, outputSet)
