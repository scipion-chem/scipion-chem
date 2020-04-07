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

import copy
import os
import sys

from pwem.protocols import EMProtocol
from pyworkflow.object import Float, Integer
from pyworkflow.protocol.params import PointerParam, EnumParam, MultiPointerParam, BooleanParam, StringParam
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtBioinformaticsListOperate(EMProtocol):
    """This protocol will remove all duplicated entries using the DbID as key"""
    _label = 'operate list'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('operation', EnumParam, choices=['Unique', 'Union', 'Intersection', 'Difference', 'Change DbID',
                                                       'Keep columns', 'Filter'],
                      label='Operation', default=0,
                      help='Unique: Remove replicated Ids.')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation!=1)')
        form.addParam('multipleInputListID', MultiPointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation==1)')
        form.addParam('inputListID2', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation==2 or operation==3)')
        form.addParam('newDb', StringParam,
                       label='New Db:', condition='(operation==4)',
                       help='New database. It can be a label or one of the columns in the table')
        form.addParam('newDbId', StringParam,
                       label='New DbID:', condition='(operation==4)',
                       help='It must be one of the existing labels in the database list')
        form.addParam('keepColumns', StringParam,
                       label='Keep columns:', condition='(operation==5)',
                       help='They must exist in the input database list. Separated by semicolons')
        form.addParam('filterColumn', StringParam,
                       label='Filter column:', condition='(operation==6)',
                       help='It must exist in the input database list.')
        form.addParam('filterOp', EnumParam, choices=['==', '>', '>=', '<', '<=', '!=', 'startswith',
                                                      'endswith', 'contains', 'does not startwith',
                                                      'does not end with', 'does not contain'],
                       label='Filter operation:', condition='(operation==6)')
        form.addParam('filterValue', StringParam,
                       label='Value:', condition='(operation==6)',
                       help='Value to use in the filter')
        form.addParam('removeDuplicates', BooleanParam, default=False,
                       label='Remove duplicates:', condition='(operation!=1)')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        outputDict = {}
        if self.operation.get()==1:
            # Union
            for database in self.multipleInputListID:
                for databaseEntry in database.get():
                    add=True
                    if self.removeDuplicates.get():
                        add=not databaseEntry.getDbId() in outputDict
                    if add:
                        dbEntry = DatabaseID()
                        dbEntry.copy(databaseEntry, copyId=False)
                        outputDict[databaseEntry.getDbId()]=dbEntry
        elif self.operation.get()==0 or self.operation.get()==2 or self.operation.get()==3:
            # Unique, Intersection, Difference
            outputList2 = []
            if self.operation.get()==2 or self.operation.get()==3:
                for databaseEntry in self.inputListID2.get():
                    outputList2.append(databaseEntry.getDbId())

            for databaseEntry in self.inputListID.get():
                add=False
                if self.operation.get()==0: # Unique
                    add = not databaseEntry.getDbId() in outputDict
                elif self.operation.get()==2: # Intersection
                    add = databaseEntry.getDbId() in outputList2
                    if self.removeDuplicates.get():
                        add=add and not databaseEntry.getDbId() in outputDict
                elif self.operation.get()==3: # Difference
                    add = not databaseEntry.getDbId() in outputList2
                    if self.removeDuplicates.get():
                        add = add and not databaseEntry.getDbId() in outputDict
                if add:
                    dbEntry = DatabaseID()
                    dbEntry.copy(databaseEntry)
                    outputDict[databaseEntry.getDbId()]=dbEntry
        elif self.operation.get()==4:
            # Change ID
            newLabel=True
            for name, _ in self.inputListID.get().getFirstItem().getAttributes():
                if self.newDb.get()==name:
                    newLabel=False
                    break

            for databaseEntry in self.inputListID.get():
                dbEntry = DatabaseID()
                dbEntry.copy(databaseEntry)
                if hasattr(dbEntry,self.newDbId.get()):
                    if newLabel:
                        dbEntry.setDatabase(self.newDb.get())
                    else:
                        dbEntry.setDatabase(dbEntry.getAttributeValue(self.newDb.get()))
                    dbEntry.setDbId(dbEntry.getAttributeValue(self.newDbId.get()))
                add = True
                if self.removeDuplicates.get():
                    add = add and not dbEntry.getDbId() in outputDict
                if add:
                    outputDict[dbEntry.getDbId()] = dbEntry
        elif self.operation.get()==5:
            # Keep columns
            keepList=[x.strip() for x in self.keepColumns.get().split()]
            keepList.append("database")
            keepList.append("dbId")

            ignoreList=[]
            for name, _ in self.inputListID.get().getFirstItem().getAttributes():
                if not name in keepList:
                    ignoreList.append(name)
            for databaseEntry in self.inputListID.get():
                dbEntry = DatabaseID()
                dbEntry.copy(databaseEntry,ignoreAttrs=ignoreList)
                add = True
                if self.removeDuplicates.get():
                    add = add and not dbEntry.getDbId() in outputDict
                if add:
                    outputDict[dbEntry.getDbId()] = dbEntry
        elif self.operation.get()==6:
            # Filter columns
            referenceValue = self.filterValue.get()
            value = self.inputListID.get().getFirstItem().getAttributeValue(self.filterColumn.get())
            if isinstance(value,float):
                referenceValue = float(referenceValue)
            elif isinstance(value,int):
                referenceValue = int(referenceValue)

            for databaseEntry in self.inputListID.get():
                dbEntry = DatabaseID()
                dbEntry.copy(databaseEntry)
                add = False

                value = dbEntry.getAttributeValue(self.filterColumn.get())
                if isinstance(value, Float):
                    value = float(value)
                elif isinstance(value, Integer):
                    value = int(value)

                filterOp = self.filterOp.get()
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

                if self.removeDuplicates.get():
                    add = add and not dbEntry.getDbId() in outputDict
                if add:
                    outputDict[dbEntry.getDbId()] = dbEntry

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for dbId in outputDict:
            outputDatabaseID.append(outputDict[dbId])
        self._defineOutputs(output=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
