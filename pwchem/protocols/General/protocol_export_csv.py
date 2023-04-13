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
    """Export a set as a csv. It is located in the Run directory"""
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
