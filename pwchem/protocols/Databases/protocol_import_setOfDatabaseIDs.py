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
    """
    This protocol imports a set of database identifiers from a plain text file
    and converts them into a structured SetOfDatabaseID object for use in
    downstream bioinformatics and chemoinformatics workflows.

    Each line of the input file is expected to contain a single identifier
    belonging to a specified biological or chemical database.

    Inputs
    ------
    filePath:
        Path to a text file containing database identifiers.
        Each line represents one ID.

    databaseName:
        Name of the source database to assign to all imported IDs
        (e.g., UniProt, ChEMBL, PDB, BindingDB, etc.).

    Workflow
    --------
    1. Input reading
       - Opens the provided text file
       - Reads identifiers line by line

    2. Object creation
       - Initializes an empty SetOfDatabaseID container
       - Converts each line into a DatabaseID object:
         - dbId = trimmed line content
         - database = user-defined databaseName

    3. Population step
       - Appends each DatabaseID object into the collection

    4. Output generation
       - Produces a structured SetOfDatabaseID object
       - Stores it in the protocol output directory

    Output
    ------
    outputDatabaseIDs:
        SetOfDatabaseID object containing:
        - All imported identifiers
        - Associated database label for each entry

    Summary
    -------
    This protocol provides a simple interface to import external identifier
    lists into the Scipion/PWEM ecosystem, enabling interoperability with
    other protocols that require structured database ID inputs.

    Notes
    -----
    - The input file must contain one ID per line.
    - No validation of ID format is performed beyond trimming whitespace.
    - All IDs are assigned the same database name.
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

