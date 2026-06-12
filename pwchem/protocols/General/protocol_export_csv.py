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
from pyworkflow.protocol.params import PointerParam

class ProtChemExportCSV(EMProtocol):
    """
    AI Generated:

    Protocol to export a Scipion EMSet into a CSV file.

    This protocol iterates over all elements in a given EMSet and writes
    their attributes into a structured CSV file located in the protocol
    run directory.

    Inputs
    ------
    inputSet:
        Input dataset (EMSet) to be exported. Each element in the set
        must provide attributes accessible via getAttributes().

    CSV format
    ----------
    - The first row contains the header with attribute names.
    - Each subsequent row corresponds to one element in the set.
    - Values are separated by semicolons ("; ").
    - Attribute values are extracted using value.get().

    Workflow
    --------
    1. Open output CSV file in the run directory.

    2. Iterate over entries in the input set:
       - Retrieve attribute key-value pairs using getAttributes().
       - For the first entry:
         * Build header line using attribute names.
       - For all entries:
         * Extract attribute values.
         * Convert values to string.
         * Concatenate using ";" separator.

    3. Write:
       - Header line (only once, for the first element)
       - One line per entry with attribute values

    4. Close file.

    Output
    ------
    output.csv:
        CSV file stored in the protocol run directory containing
        all elements and their attributes.

    Error handling
    --------------
    - Assumes all entries share the same attribute structure.
    - Does not validate missing or inconsistent attributes across entries.
    - If attributes differ between entries, columns may become inconsistent.

    Validation
    ----------
    - Requires a valid EMSet input.
    - Each entry must implement getAttributes() returning key-value pairs.

    Summary
    -------
    This protocol provides a simple way to serialize Scipion datasets into
    a human-readable CSV format, facilitating:
    - Data inspection
    - External analysis (Excel, pandas, R, etc.)
    - Data exchange between workflows

    Notes
    -----
    - Separator is fixed to "; " (semicolon + space).
    - Output file name is always "output.csv".
    - No quoting or escaping is applied to values.
    - Recommended for simple datasets with flat attribute structures.
    """
    _label = 'export csv'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="EMSet",
                       label='Set:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('exportStep')

    def exportStep(self):
        fh = open(self._getPath("output.csv"),"w")
        lineNo = 0
        lineHdr = ""
        for entry in self.inputSet.get():
            line = ""
            for key, value in entry.getAttributes():
                if lineNo == 0:
                    if lineHdr!="":
                        lineHdr+="; "
                    lineHdr += key
                if line!="":
                    line += "; "
                line+=str(value.get())
            if lineNo==0:
                fh.write(lineHdr+"\n")
            fh.write(line+"\n")
            lineNo+=1
        fh.close()
