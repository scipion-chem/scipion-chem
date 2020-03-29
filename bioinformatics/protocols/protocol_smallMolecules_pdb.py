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
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import (PointerParam)
from bioinformatics.objects import DatabaseID, SetOfDatabaseID

class ProtAtomStructSmallMoleculesPDB(EMProtocol):
    """Extract a set of PDB Ids from the interaction list of a set of small molecules"""
    _label = 'small mols -> pdb'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of small molecules:', allowsNull=False,
                       help="List of small molecules to extract the list of PDB IDs")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractStep')

    def extractStep(self):
        listIds={}
        smallMolID = "none"
        for item in self.inputListID.get():
            if hasattr(item,"_iteractsWithPDBId"):
                tokens = item._iteractsWithPDBId.get().split(";")
                for token in tokens:
                    pdbId = token.strip()
                    if not pdbId in listIds:
                        listIds[pdbId]=[]
                    chemId=item.getDbId()
                    if hasattr(item,"_PDBChemId"):
                        chemId=item._PDBChemId.get()
                        smallMolID = "pdbchem"
                    listIds[pdbId].append(chemId)

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath(),suffix='PDBs')
        for pdbId in listIds:
            pdb = DatabaseID()
            pdb.setDatabase("pdb")
            pdb.setDbId(pdbId)
            pdb._pdbId = pwobj.String(pdbId)
            pdb._PDBLink = pwobj.String("https://www.rcsb.org/structure/%s" % pdbId)
            aux = " ; ".join(listIds[pdbId])
            if smallMolID == "none":
                pdb._interactsWithChemId = pwobj.String(aux)
            elif smallMolID == "pdbchem":
                pdb._interactsWithPDBChemId = pwobj.String(aux)
            outputDatabaseID.append(pdb)
        self._defineOutputs(outputPDBs=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
