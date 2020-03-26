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
from pyworkflow.protocol.params import PointerParam
from atomstructutilsWeb.objects import DatabaseID, SetOfDatabaseID

class ProtAtomStructListUnique(EMProtocol):
    """This protocol will remove all duplicated entries using the DbID as key"""
    _label = 'unique list'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of DB Ids:', allowsNull=False,
                       help="List of identifiers to simplify")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('uniqueStep',self.inputListID.getObjId())

    def uniqueStep(self, objId):
        ligandDict = {}
        for databaseEntry in self.inputListID.get():
            if not databaseEntry.getDbId() in ligandDict:
                dbEntry = DatabaseID()
                dbEntry.copy(databaseEntry)
                ligandDict[databaseEntry.getDbId()]=dbEntry

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for dbId in ligandDict:
            outputDatabaseID.append(ligandDict[dbId])
        self._defineOutputs(output=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
