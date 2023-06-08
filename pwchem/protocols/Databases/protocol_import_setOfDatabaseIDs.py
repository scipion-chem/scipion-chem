# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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
from pyworkflow.protocol.params import PathParam, StringParam
from pwchem.objects import DatabaseID, SetOfDatabaseID

class ProtChemImportSetOfDatabaseIDs(EMProtocol):
    """Import a set of databaseIDs from a text file with each ID in a line
    """
    _label = 'Import set of DatabaseIDs'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('filePath', PathParam,
                      label='IDs file: ',
                      help='File containing a set of IDs to import')

        form.addParam('databaseName', StringParam,
                      label='Database name: ',
                      help='Name of the database the IDs come from')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        outputDatabaseIDs = SetOfDatabaseID().create(outputPath=self._getPath())
        with open(self.filePath.get()) as f:
            for line in f:
                outputDatabaseIDs.append(DatabaseID(dbId=line.strip(), database=self.databaseName.get()))

        self._defineOutputs(outputDatabaseIDs=outputDatabaseIDs)

