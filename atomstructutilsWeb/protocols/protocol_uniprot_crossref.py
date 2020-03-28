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

import lxml.etree as ET
import os
import sys
import urllib.request

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import PointerParam, EnumParam
from atomstructutilsWeb.objects import DatabaseID, SetOfDatabaseID

class ProtAtomStructUniprotCrossRef(EMProtocol):
    """Extract cross references from uniprot"""
    _label = 'uniprot crossref'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of Uniprot Ids:', allowsNull=False,
                       help="List of atomic structures for the query")
        form.addParam('extract', EnumParam, choices=['PDB'], default=0,
                       label='What to extract:')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractStep')

    def extractStep(self):
        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for item in self.inputListID.get():

            uniprotId = item._uniprotId.get()
            print("Processing %s"%uniprotId)

            urlId = "https://www.uniprot.org/uniprot/%s.xml" % uniprotId

            fnXML=self._getExtraPath("%s.xml"%uniprotId)
            if not os.path.exists(fnXML):
                print("Fetching uniprot: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnXML)
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnXML):
                tree = ET.parse(fnXML)

                outputId = []
                for child in tree.getroot().iter():
                    if child.tag.endswith("dbReference"):
                        if self.extract.get()==0:
                            if child.attrib['type']=='PDB':
                                outputId.append(child.attrib['id'])
                if len(outputId)>0:
                    for outId in outputId:
                        newItem = DatabaseID()
                        newItem.copy(item, copyId=False)
                        if self.extract.get()==0:
                            newItem._pdbId = pwobj.String(outId)
                            newItem._PDBLink = pwobj.String("https://www.rcsb.org/structure/%s" % outId)
                        outputDatabaseID.append(newItem)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)

    def _validate(self):
        errors=[]
        if not hasattr(self.inputListID.get().getFirstItem(),"_uniprotId"):
            errors.append("The set does not have an _uniprotId")
        return errors