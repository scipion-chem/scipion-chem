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
from pyworkflow.protocol.params import (PointerParam)
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtBioinformaticsPDBUniprot(EMProtocol):
    """Query PDB for Uniprot sequences related to these proteins"""
    _label = 'pdb -> uniprot'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of PDB Ids:', allowsNull=False,
                       help="List of atomic structures for the query")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath())
        for item in self.inputListID.get():
            newItem = DatabaseID()
            newItem.copy(item)
            newItem._uniprotId = pwobj.String("Not available")
            newItem._uniprotLink = pwobj.String("Not available")

            pdbId = item._pdbId.get()
            print("Processing %s"%pdbId)

            urlId = "https://www.rcsb.org/pdb/rest/das/pdb_uniprot_mapping/alignment?query=%s" % pdbId
            if hasattr(item,"_chain"):
                urlId+="."+item._chain.get().upper()

            fnXml=self._getExtraPath("%s.xml"%pdbId)
            if not os.path.exists(fnXml):
                print("Fetching uniprot: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnXml)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnXml):
                try:
                    tree = ET.parse(fnXml)
                    # print(ET.tostring(tree, pretty_print=True))

                    uniprotId = None
                    for child in tree.getroot().iter():
                        if child.tag.endswith("alignObject"):
                            if child.attrib['dbSource']=="UniProt":
                                uniprotId=child.attrib['dbAccessionId']
                                break
                    if uniprotId:
                        newItem._uniprotId = pwobj.String(uniprotId)
                        newItem._uniprotLink = pwobj.String("https://www.uniprot.org/uniprot/%s"%uniprotId)
                except:
                    print("    Cannot parse the Uniprot XML: %s"%fnXml)

                outputDatabaseID.append(newItem)

        self._defineOutputs(outputUniprot=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
