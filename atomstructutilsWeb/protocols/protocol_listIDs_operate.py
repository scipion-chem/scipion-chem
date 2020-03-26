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
from pyworkflow.protocol.params import PointerParam, EnumParam, MultiPointerParam, BooleanParam, StringParam
from atomstructutilsWeb.objects import DatabaseID, SetOfDatabaseID

class ProtAtomStructListOperate(EMProtocol):
    """This protocol will remove all duplicated entries using the DbID as key"""
    _label = 'operate list'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('operation', EnumParam, choices=['Unique', 'Union', 'Intersection', 'Difference', 'Change DbID'],
                      label='Operation', default=0,
                      help='Unique: Remove replicated Ids.')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation!=1)')
        form.addParam('multipleInputListID', MultiPointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation==1)')
        form.addParam('inputListID2', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=True, condition='(operation==2 or operation==3)')
        form.addParam('newDbId', StringParam,
                       label='New DbID:', condition='(operation==4)',
                       help='It must be one of the existing labels in the database list')
        form.addParam('removeDuplicates', BooleanParam, default=False,
                       label='Remove duplicates:', condition='(operation!=1)')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        ligandDict = {}
        if self.operation.get()==1:
            # Union
            for database in self.multipleInputListID:
                for databaseEntry in database.get():
                    add=True
                    if self.removeDuplicates.get():
                        add=not databaseEntry.getDbId() in ligandDict
                    if add:
                        dbEntry = DatabaseID()
                        dbEntry.copy(databaseEntry, copyId=False)
                        ligandDict[databaseEntry.getDbId()]=dbEntry
        elif self.operation.get()==0 or self.operation.get()==2 or self.operation.get()==3:
            ligandList2 = []
            if self.operation.get()==2 or self.operation.get()==3:
                for databaseEntry in self.inputListID2.get():
                    ligandList2.append(databaseEntry.getDbId())

            for databaseEntry in self.inputListID.get():
                add=False
                if self.operation.get()==0: # Unique
                    add = not databaseEntry.getDbId() in ligandDict
                elif self.operation.get()==2: # Intersection
                    add = databaseEntry.getDbId() in ligandList2
                    if self.removeDuplicates.get():
                        add=add and not databaseEntry.getDbId() in ligandDict
                elif self.operation.get()==3: # Difference
                    add = not databaseEntry.getDbId() in ligandList2
                    if self.removeDuplicates.get():
                        add = add and not databaseEntry.getDbId() in ligandDict
                if add:
                    dbEntry = DatabaseID()
                    dbEntry.copy(databaseEntry)
                    ligandDict[databaseEntry.getDbId()]=dbEntry
        elif self.operation.get()==4: # Change ID
            for databaseEntry in self.inputListID.get():
                dbEntry = DatabaseID()
                dbEntry.copy(databaseEntry)
                if hasattr(dbEntry,self.newDbId.get()):
                    dbEntry.setDbId(dbEntry.getAttributeValue(self.newDbId.get()))
                add = True
                if self.removeDuplicates.get():
                    add = add and not dbEntry.getDbId() in ligandDict
                if add:
                    ligandDict[dbEntry.getDbId()] = dbEntry

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for dbId in ligandDict:
            outputDatabaseID.append(ligandDict[dbId])
        self._defineOutputs(output=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
