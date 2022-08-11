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

import contextlib
import lxml.etree as ET
import os
import urllib.request

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam
from pwchem.utils import runOpenBabel

class ProtChemPubChemSearch(EMProtocol):
    """Add the best batching entry from Pubchem https://pubchem.ncbi.nlm.nih.gov/"""
    _label = 'PubChem'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set to filter:', allowsNull=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('operateStep')

    def convertInputStep(self):
        self.convMols = []
        for mol in self.inputSet.get():
            fnSmall = mol.getFileName()
            fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))
            if ext != '.smi':
                args = " -i%s %s --partialcharge none -O %s.smi" % (ext[1:], os.path.abspath(fnSmall), fnRoot)
                runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getExtraPath()))
                mol.setFileName(self._getExtraPath(fnRoot + '.smi'))
            self.convMols.append(mol.clone())


    def operateStep(self):
        outputSet = self.inputSet.get().create(self._getPath())
        for mol in self.convMols:
            fnSmall = mol.smallMoleculeFile.get()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            fnName = self._getExtraPath(fnBase+".txt")
            cid = None
            if not os.path.exists(fnName):
                if fnSmall.endswith('.smi'):
                    fhSmile = open(fnSmall)
                    smile = fhSmile.readlines()[0].split()[0].strip() # Only first line
                    fhSmile.close()
                    url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT"%smile
                    print(url)

                pubChemName = ""
                try:
                    with contextlib.closing(urllib.request.urlopen(url)) as fp:
                        mybytes = fp.read()
                        cid = mybytes.decode("utf8").split()[0]
                        fp.close()

                    if cid!="0":
                        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/%s/XML/?response_type=save&response_basename=compound_CID_%s"%(cid,cid)
                        print(url)
                        fnXml = self._getTmpPath("compound_{}.xml".format(cid))
                        urllib.request.urlretrieve(url, fnXml)
                        if os.path.exists(fnXml):
                            tree = ET.parse(fnXml)

                        pubChemName = ""
                        for child in tree.getroot().iter():
                            if "RecordTitle" in child.tag:
                                pubChemName = child.text
                                print(pubChemName)
                                fh = open(fnName,'w')
                                fh.write(pubChemName+" ;; %s"%cid)
                                fh.close()
                                break
                except Exception as e:
                    print(e)
                    print("  Could not be retrieved")
            else:
                fh = open(fnName)
                tokens = fh.readline().split(';;')
                pubChemName = tokens[0].strip()
                cid = tokens[1].strip()
                fh.close()


            mol.pubChemName = pwobj.String(pubChemName)
            if cid is not None and cid!="0":
                url = "https://pubchem.ncbi.nlm.nih.gov/compound/%s"%cid
            else:
                url = ""
            mol.pubChemURL = pwobj.String(url)
            outputSet.append(mol)

        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)

    def _validate(self):
        errors = []
        firstItem = self.inputSet.get().getFirstItem()
        if not hasattr(firstItem, "smallMoleculeFile"):
            errors.append("The input set does not contain small molecules")
        return errors
