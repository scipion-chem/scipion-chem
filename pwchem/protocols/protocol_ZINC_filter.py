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
import os
import sys
import urllib.request

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, BooleanParam, EnumParam

class ProtBioinformaticsZINCFilter(EMProtocol):
    """Filter a set of small molecules by being in all selected catalogs of ZINC.
       See https://zinc15.docking.org/substances/subsets/"""
    _label = 'ZINC filter'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set to filter:', allowsNull=False)
        form.addParam('mode', EnumParam, choices=["Remove if included", "Remove if not included"],
                      label='Mode', default=0)
        form.addParam('notForSale', BooleanParam, label='Nor for sale', default=True)
        form.addParam('agent', BooleanParam, label='Agent', default=False)
        form.addParam('forSale', BooleanParam, label='For sale', default=False)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        outputSet = self.inputSet.get().create(self._getPath())
        for oldEntry in self.inputSet.get():
            fnSmall = oldEntry.smallMoleculeFile.get()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            if "-" in fnBase:
                fnBase = fnBase.split("-")[0]
            if not "ZINC" in fnBase:
                continue
            fnAdd = self._getExtraPath(fnBase+".txt")
            if not os.path.exists(fnAdd):
                url="http://zinc15.docking.org/substances/%s"%fnBase
                print(url)
                sys.stdout.flush()
                add = True
                try:
                    with contextlib.closing(urllib.request.urlopen(url)) as fp:
                        mybytes = fp.read()
                        mystr = mybytes.decode("utf8")
                        fp.close()
                    notForSale=False
                    agent=False
                    forSale=False
                    inTitle = False
                    title = ""
                    for line in mystr.split('\n'):
                        if "/substances/subsets/not-for-sale/" in line:
                            notForSale = True
                        elif "/substances/subsets/agent/" in line:
                            agent = True
                        elif "/substances/subsets/for-sale/" in line:
                            forSale = True
                        if "</title>" in line:
                            inTitle = False
                        if inTitle and not fnBase in line:
                            title+=line.strip()
                        if "<title>" in line:
                            inTitle = True
                    print("  Title: %s"%title)
                    print("  Not for sale: %s"%str(notForSale))
                    print("  Agent: %s"%str(agent))
                    print("  For sale: %s"%str(forSale))
                    if self.mode.get()==0:
                        if notForSale and self.notForSale.get():
                            add=False
                        if agent and self.agent.get():
                            add=False
                        if forSale and self.forSale.get():
                            add=False
                    else:
                        if not notForSale and self.notForSale.get():
                            add=False
                        if not agent and self.agent.get():
                            add=False
                        if not forSale and self.forSale.get():
                            add=False
                except:
                    add = True
                    title = "Could not retrieve from ZINC"
                    print("  Could not be retrieved")
                fh = open(fnAdd,'w')
                fh.write(str(add)+" ;; "+title)
                fh.close()
            else:
                fh = open(fnAdd)
                tokens = fh.readline().split(';;')
                add = tokens[0].strip()=="True"
                if len(tokens)>1:
                    title = tokens[1].strip()
                else:
                    title = ""
                fh.close()

            if add:
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(oldEntry)
                newEntry.ZINCname = pwobj.String(title)
                outputSet.append(newEntry)

        if len(outputSet)>0:
            self._defineOutputs(output=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)

    def _validate(self):
        errors = []
        firstItem = self.inputSet.get().getFirstItem()
        if not hasattr(firstItem,"smallMoleculeFile"):
            errors.append("The input set does not contain small molecules")
        elif not "ZINC" in firstItem.smallMoleculeFile.get():
            errors.append("Cannot find the ZINC code in the small molecule filename")
        return errors