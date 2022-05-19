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

from math import ceil
import numpy as np

from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer
from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam, IntParam, FloatParam

class ProtChemListOperate(EMProtocol):
    """Filter a set by a column value or keep just a few columns"""
    _label = 'operate set'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('operation', EnumParam, choices=['Filter', 'Keep columns', 'Unique', 'Top N', 'Bottom N',
                                                       'Top %', 'Bottom %', 'Count', 'Intersection', 'Sort'],
                      label='Operation', default=0,
                      help='In intersection, we keep those entries of the Set to filter whose identifier (filter column) '
                           'are in the second set')
        form.addParam('inputSet', PointerParam, pointerClass="EMSet",
                       label='Set to filter:', allowsNull=False)
        form.addParam('secondSet', PointerParam, pointerClass="EMSet", condition="operation==8", # Intersection between 2 Sets
                       label='Set with IDs:', allowsNull=True)
        form.addParam('filterColumn', StringParam,
                       label='Filter column:', condition='(operation!=1)',
                       help='It must exist in the input object.')
        form.addParam('filterOp', EnumParam, choices=['==', '>', '>=', '<', '<=', '!=', 'between', 'startswith',
                                                      'endswith', 'contains', 'does not startwith',
                                                      'does not end with', 'does not contain'],
                       label='Filter operation:', condition='(operation==0)')
        form.addParam('filterValue', StringParam,
                       label='Value:', condition='(operation==0)',
                       help='Value to use in the filter')
        form.addParam('filterValue2', StringParam,
                       label='Lower Value:', condition='(operation==0 and filterOp==6)',
                       help='Value to use in the filter')
        form.addParam('keepColumns', StringParam,
                       label='Keep columns:', condition='(operation==1)',
                       help='They must exist in the input database list. Separated by semicolons (e.g. column1 ; column2 ; ...)')
        form.addParam('N', IntParam,
                       label='N:', default=500, condition='(operation==3 or operation==4)')
        form.addParam('percentile', FloatParam,
                       label='Percentile:', default=5, condition='(operation==5 or operation==6)',
                       help='Between 0 and 100')
        form.addParam('direction', EnumParam, choices=['Ascending', 'Descending'], default=0,
                       label='Sorting direction:', condition='(operation==9)')

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

            filterOp = self.filterOp.get()
            if filterOp == 6:
                referenceValue2 = self.filterValue2.get()

                if isinstance(value, float):
                    referenceValue2 = float(referenceValue2)
                elif isinstance(value, int):
                    referenceValue2 = int(referenceValue2)

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if isinstance(value, Float):
                    value = float(value)
                elif isinstance(value, Integer):
                    value = int(value)


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
                elif filterOp == 6:  # between
                    add = (value <= referenceValue and value >= referenceValue2)
                elif filterOp == 7:  #startswith
                    add = value.startswith(referenceValue)
                elif filterOp == 8:  # endswith
                    add = value.endswith(referenceValue)
                elif filterOp == 9:  # contains
                    add = referenceValue in value
                elif filterOp == 10:  # does not startswith
                    add = not (value.startswith(referenceValue))
                elif filterOp == 11:  # does not endswith
                    add = not (value.endswith(referenceValue))
                elif filterOp == 12:  # does not contains
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

        elif self.operation.get()==2:
            # Unique
            found={}
            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if not value in found:
                    found[value] = True
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()>=3 and self.operation.get()<=6:
            # Top N, Bottom N,Top %, Bottom %
            V = []
            for entry in self.inputSet.get():
                V.append(entry.getAttributeValue(self.filterColumn.get()))
            V.sort()
            op = self.operation.get()
            if op==3:
                threshold = V[-self.N.get()]
            elif op==4:
                threshold = V[self.N.get()-1]
            elif op==5:
                threshold = V[-ceil(self.percentile.get()/100*len(V))]
            elif op==6:
                threshold = V[ceil(self.percentile.get()/100*len(V))-1]

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if (op==3 or op==5) and value>=threshold:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)
                elif (op==4 or op==6) and value<=threshold:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==7:
            # Count the number of entries that are the same
            count={}
            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())
                if not value in count:
                    count[value] = 0
                count[value] += 1

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())

                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(oldEntry)
                newEntry.count = Integer(count[value])
                outputSet.append(newEntry)

        elif self.operation.get()==8:
            # Intersection between 2 Sets
            secondSet={}
            for entry in self.secondSet.get():
                value = entry.getAttributeValue(self.filterColumn.get())
                if not value in secondSet:
                    secondSet[value] = True

            for oldEntry in self.inputSet.get():
                value = oldEntry.getAttributeValue(self.filterColumn.get())

                if value in secondSet:
                    newEntry = self.inputSet.get().ITEM_TYPE()
                    newEntry.copy(oldEntry)
                    outputSet.append(newEntry)

        elif self.operation.get()==9:
            # Sort
            V = []
            newEntries = []
            for entry in self.inputSet.get():
                V.append(entry.getAttributeValue(self.filterColumn.get()))
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(entry)
                newEntry.cleanObjId()
                newEntries.append(newEntry)
            if self.direction.get()==0:
                idxSort = np.argsort(V)
            else:
                idxSort = np.argsort(-np.asarray(V))

            for idx in idxSort:
                outputSet.append(newEntries[idx])

        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)
