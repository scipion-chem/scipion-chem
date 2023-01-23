# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import json
from urllib.request import urlopen

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol import params

from pwchem.objects import DatabaseID, SetOfDatabaseID
from pwchem.utils import performBatchThreading


dbChoices = ['PDB', 'UniProtKB', 'EMBL', 'ChEMBL', 'DrugBank', 'BindingDB', 'GO', 'Gene3D', 'InterPro', 'Pfam',
             'DOI', 'PubMed', 'Other']

class ProtChemUniprotCrossRef(EMProtocol):
    """Extract cross references from uniprot"""
    _label = 'uniprot crossref'

    def __init__(self, *args, **kwargs):
        EMProtocol.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputListID', params.PointerParam, pointerClass="SetOfDatabaseID",
                      label='List of Uniprot Ids:', allowsNull=False,
                      help="List of atomic structures for the query")
        form.addParam('extract', params.EnumParam, default=0,
                      choices=dbChoices, label='Database info to extract: ',
                      help='Choose among some of these database examples or decide your own')
        form.addParam('extractOther', params.StringParam, default='', label='Other database info to extract: ',
                      condition='extract=={}'.format(len(dbChoices) -1),
                      help='Specify the name of the database to extract info from. Uniprot must have information from '
                            'this database to extract any information')

        form.addParam('storeProps', params.BooleanParam, default=False, label='Store properties: ',
                      help='Whether to store additional properties of the crossreferences')
        form.addParam('crossMain', params.BooleanParam, default=False, label='Cross reference as main ID: ',
                      help='Whether to store the crossreference IDs as the main ID of the output or as crossRefID')
        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('extractStep')

    def extractStep(self):
        nt = self.numberOfThreads.get()
        allProps, dbName = set([]), self.getDBName()
        outputIds = performBatchThreading(self.fetchCrossRef, self.inputListID.get(), nt, cloneItem=True,
                                          allProps=allProps)
        outputIds = {item: itemDic[item] for itemDic in outputIds for item in itemDic}

        print()
        outputDatabaseIDs = SetOfDatabaseID().create(self._getPath())
        if len(outputIds)>0:
            for item in outputIds:
                uniprotId = item.getDbId()

                if len(outputIds[item]) == 0:
                    print('Warning: {} uniprot ID has no cross references with {} database so it will not be included'
                          ' in the output'.format(uniprotId, dbName))

                for crossId in outputIds[item]:
                    newItem = DatabaseID()
                    newItem.copy(item, copyId=False)
                    if self.crossMain.get():
                        newItem._crossId = pwobj.String(uniprotId)
                        newItem._crossDatabase = pwobj.String('UniProt')
                        newItem.setDbId(crossId)
                        newItem.setDatabase(dbName)
                        prefix = ''
                    else:
                        newItem._crossId = pwobj.String(crossId)
                        newItem._crossDatabase = pwobj.String(dbName)
                        prefix = 'cross_'

                    if self.storeProps.get():
                        for prop in allProps:
                            if prop in outputIds[item][crossId]:
                                value = pwobj.String(outputIds[item][crossId][prop])
                            else:
                                value = pwobj.String('None')
                            setattr(newItem, '_{}{}'.format(prefix, prop.replace(' ', '_')), value)

                    outputDatabaseIDs.append(newItem)

        self._defineOutputs(outputDatabaseIDs=outputDatabaseIDs)
        self._defineSourceRelation(self.inputListID, outputDatabaseIDs)

    def _validate(self):
        errors=[]
        for dbId in self.inputListID.get():
            if not dbId.getDatabase().lower() == "uniprot":
                errors.append("{} is not labeled as an Uniprot ID".format(dbId.getDbId()))
                break
        return errors

    def getDBName(self):
        if self.extract.get() != len(dbChoices) - 1:
            return self.getEnumText('extract')
        else:
            return self.extractOther.get()

    def fetchCrossRef(self, dbIds, outputLists, it, allProps):
        for item in dbIds:
            nItem = item.clone()
            outputIds = {nItem: {}}
            uniprotId = item.getDbId()

            jsonDic = self.fetchUniprotJSON(uniprotId)
            for cross in jsonDic['uniProtKBCrossReferences']:
                if cross['database'] == self.getDBName():
                    crossId = cross['id']
                    outputIds[nItem][crossId] = {}
                    if self.storeProps.get():
                        for propDic in cross['properties']:
                            prop = propDic["key"]
                            allProps.add(prop)
                            outputIds[nItem][crossId][prop] = propDic["value"]

            outputLists[it].append(outputIds)
        return outputLists[it]

    def fetchUniprotJSON(self, uniprotId):
        print("Processing %s" % uniprotId)
        url = "https://www.uniprot.org/uniprot/%s.json" % uniprotId
        with urlopen(url) as response:
            jDic = json.loads(response.read())
        return jDic
