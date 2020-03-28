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

import os
import sys
import urllib.request

from pwem.protocols import EMProtocol
from pwem.convert.sequence import sequenceLength
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam
from atomstructutilsWeb.objects import DatabaseID, SetOfDatabaseID, ProteinSequenceFile

class ProtAtomStructUniprotDownload(EMProtocol):
    """Download the Fasta files of a set of uniprotId's"""
    _label = 'uniprot download'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of Uniprot Ids:', allowsNull=False,
                       help="List of atomic structures for the query")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        fnList = []
        for item in self.inputListID.get():
            newItem = DatabaseID()
            newItem.copy(item)
            newItem._uniprotFile = pwobj.String("Not available")
            newItem._unitprotSeqLength = pwobj.Integer(-1)

            uniprotId = item._uniprotId.get()
            print("Processing %s"%uniprotId)

            urlId = "https://www.uniprot.org/uniprot/%s.fasta" % uniprotId

            fnFasta=self._getExtraPath("%s.fasta"%uniprotId)
            if not os.path.exists(fnFasta):
                print("Fetching uniprot: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnFasta)
                        if not fnFasta in fnList:
                            fnList.append(fnFasta)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnFasta):
                newItem._uniprotFile = pwobj.String(fnFasta)
                newItem._unitprotSeqLength = pwobj.Integer(sequenceLength(fnFasta))

            outputDatabaseID.append(newItem)

        fnAll = self._getPath("sequences.fasta")
        with open(fnAll, 'w') as outfile:
            for fname in fnList:
                with open(fname) as infile:
                    for line in infile:
                        outfile.write(line)
                    outfile.write('\n\n')
        seqFile=ProteinSequenceFile()
        seqFile.setFileName(fnAll)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
        self._defineOutputs(outputSequence=seqFile)
        self._defineSourceRelation(self.inputListID, seqFile)

    def _validate(self):
        errors=[]
        if not hasattr(self.inputListID.get().getFirstItem(),"_uniprotId"):
            errors.append("The set does not have an _uniprotId")
        return errors