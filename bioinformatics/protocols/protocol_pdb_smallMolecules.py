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
from bioinformatics import Plugin

class ProtBioinformaticsPDBSmallMolecules(EMProtocol):
    """Query PDB for small molecules of a list of PDB ids"""
    _label = 'pdb -> small mols'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', PointerParam, pointerClass="SetOfDatabaseID",
                       label='List of PDB Ids:', allowsNull=False,
                       help="List of atomic structures for the query")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        listIds=[]
        for item in self.inputListID.get():
            pdbId = item._pdbId.get()
            if pdbId not in listIds:
                listIds.append(pdbId)

        ligandDict = {}
        currentLigand = None
        ligandName = None
        for pdbId in listIds:
            print("Processing %s"%pdbId)
            fnXml=self._getExtraPath("%s.xml"%pdbId)
            if not os.path.exists(fnXml):
                urlId="https://www.rcsb.org/pdb/rest/ligandInfo?structureId=%s"%pdbId
                print("Fetching ligands: %s"%urlId)
                for i in range(3):
                    try:
                        urllib.request.urlretrieve(urlId,fnXml)
                        break
                    except: # The library raises an exception when the web is not found
                        pass
            if os.path.exists(fnXml):
                tree = ET.parse(fnXml)
                # print(ET.tostring(tree, pretty_print=True))

                for child in tree.getroot().iter():
                    if child.tag=="ligand":
                        if ligandName:
                            ligandDict[ligandName] = currentLigand

                        newLigandName=child.attrib["chemicalID"]

                        if not newLigandName in ligandDict:
                            currentLigand = DatabaseID()
                            currentLigand.setDatabase("pdb")
                            currentLigand.setDbId(newLigandName)
                            currentLigand._PDBChemId=pwobj.String(newLigandName)
                            currentLigand._PDBLink=pwobj.String("https://www.rcsb.org/ligand/%s"%newLigandName)
                            currentLigand._iteractsWithPDBId=pwobj.String(pdbId)
                            if "type" in child.attrib:
                                currentLigand._PDBLigandType=pwobj.String(child.attrib["type"])
                            else:
                                currentLigand._PDBLigandType = pwobj.String("Not available")
                            if "molecularWeight" in child.attrib:
                                currentLigand._PDBLigandMolWeight = pwobj.Float(child.attrib["molecularWeight"])
                            else:
                                currentLigand._PDBLigandMolWeight = pwobj.Float(-1)
                            ligandName = newLigandName
                        else:
                            ligandDict[newLigandName]._iteractsWithPDBId.set(ligandDict[newLigandName]._iteractsWithPDBId.get()+" ; "+pdbId)
                            ligandName = None # Skip this ligand as it is already in the dictionary
                    if ligandName:
                        ctl = child.tag.lower()
                        if ctl == "chemicalname":
                            currentLigand._PDBLigandChemicalName=pwobj.String(child.text)
                        elif ctl == "formula":
                            currentLigand._PDBLigandFormula = pwobj.String(child.text)
                        elif ctl == "inchi":
                            currentLigand._PDBLigandInChI = pwobj.String(child.text)
                        elif ctl == "inchikey":
                            currentLigand._PDBLigandInChiKey = pwobj.String(child.text)
                        elif ctl == "smiles":
                            currentLigand._PDBLigandSmiles = pwobj.String(child.text)

        outputDatabaseID = SetOfDatabaseID().create(path=self._getPath(),suffix='SmallMols')
        for chemId in ligandDict:
            if not hasattr(ligandDict[chemId],"_PDBLigandChemicalName"):
                ligandDict[chemId]._PDBLigandChemicalName = pwobj.String("Not available")
            if not hasattr(ligandDict[chemId],"_PDBLigandFormula"):
                ligandDict[chemId]._PDBLigandFormula = pwobj.String("Not available")
            if not hasattr(ligandDict[chemId],"_PDBLigandInChI"):
                ligandDict[chemId]._PDBLigandInChI = pwobj.String("Not available")
            if not hasattr(ligandDict[chemId],"_PDBLigandInChiKey"):
                ligandDict[chemId]._PDBLigandInChiKey = pwobj.String("Not available")
            if not hasattr(ligandDict[chemId],"_PDBLigandSmiles"):
                ligandDict[chemId]._PDBLigandSmiles = pwobj.String("Not available")
                ligandDict[chemId]._PDBLigandSmilesImage = pwobj.String("Not available")
            else:
                fnTmp = self._getTmpPath("molecule.smiles")
                fnOut = self._getExtraPath("%s.png"%chemId)
                fh = open(fnTmp,"w")
                fh.write(ligandDict[chemId]._PDBLigandSmiles.get())
                fh.close()
                args = Plugin.getPluginHome('utils/rdkitUtils.py')+" draw %s %s"%(fnTmp,fnOut)
                Plugin.runRDKit(self,"python3",args)
                ligandDict[chemId]._PDBLigandImage = pwobj.String(fnOut)
            outputDatabaseID.append(ligandDict[chemId])
        self._defineOutputs(outputSmallMols=outputDatabaseID)
        self._defineSourceRelation(self.inputListID, outputDatabaseID)
