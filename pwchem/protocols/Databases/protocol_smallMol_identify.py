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

import os, time, json, glob
from urllib.request import urlopen

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem import Plugin
from pwchem.utils import performBatchThreading
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1
PUBCHEMID, PUBCHEMNAME, ZINC, CHEMBL = 1, 2, 3, 4
DbChoices = ['PubChemID', 'PubChemName', 'ZINC_ID', 'ChEMBL_ID']

class ProtChemSmallMolIdentify(EMProtocol):
    """Uses the PubChem search tools to identify the ligands from their SMILES.
    If no exact match is found, Tanimoto similarity check can also be performed
    Info: https://pubchem.ncbi.nlm.nih.gov/docs/pug-rest"""
    _label = 'Small molecule identification'

    def __init__(self, *args, **kwargs):
        EMProtocol.__init__(self, *args, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSet', params.PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Set to filter:', allowsNull=False)

        form.addParam('useManager', params.EnumParam, default=0, label='Manage structure using: ',
                       choices=['RDKit', 'OpenBabel'],
                       help='Whether to manage the structure (conversion to SMILES) using RDKit or OpenBabel')

        form.addParam('nameDatabase', params.EnumParam, default=0, label='Store mol name from: ',
                      choices=['None'] + DbChoices,
                      help='Replace the current molname with one of these options if available. All database values'
                           ' will be stored in the molecule attributes')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('identifyStep')

    def identifyStep(self):
        molDic = {}
        nt = self.numberOfThreads.get()
        for mol in self.inputSet.get():
            smi = self.getSMI(mol, nt)
            if smi in molDic:
                molDic[smi] += [mol.clone()]
            else:
                molDic[smi] = [mol.clone()]

        simBase, outDir = 'similarIds', os.path.abspath(self._getTmpPath())
        smisDics = performBatchThreading(self.identifySMI, molDic, nt, cloneItem=False,
                                         similarBase=os.path.join(outDir, simBase))

        # DEBUG
        # smiDic = self.identifySMI(molDic, [[]], 0, simBase)

        allSimFile = self._getPath('{}.txt'.format(simBase))
        with open(allSimFile, 'w') as f:
            f.write('SMI\tPubChemIds\n')
            for simFile in glob.glob(os.path.join(outDir, '{}_*'.format(simBase))):
                simFile = os.path.join(outDir, simFile)
                with open(simFile) as fcur:
                    f.write(fcur.read())

        outputSet = self.inputSet.get().create(self._getPath())
        for smiDic in smisDics:
            for smi in smiDic:
                for mol in molDic[smi]:
                    # Setting Ids
                    for dbName in DbChoices:
                        if dbName in smiDic[smi]:
                            setattr(mol, dbName, pwobj.String(smiDic[smi][dbName]))
                        else:
                            setattr(mol, dbName, pwobj.String(None))

                    # Setting molname
                    chdbName = self.getEnumText('nameDatabase')
                    if chdbName != 'None' and chdbName in smiDic[smi] and smiDic[smi][chdbName]:
                        mol.setMolName(smiDic[smi][chdbName])
                    outputSet.append(mol)

        if len(outputSet) > 0:
            self._defineOutputs(outputSmallMolecules=outputSet)
            self._defineSourceRelation(self.inputSet, outputSet)


    ########################### UTILS FUNCTIONS ##########################

    def identifySMI(self, smis, smisLists, it, similarBase):
        for molSMI in smis:
            cid = self.getCIDFromSmiles(molSMI)
            if cid == "0":
                print('No exact match found for SMILES: {}\nLooking for similar compounds'.format(molSMI))
                listKey = self.getSimilarityListKey(molSMI)
                simCids = self.getSimilarIds(listKey)
                cid = simCids[0]
            else:
                simCids = [cid]

            similarFile = '{}_{}.txt'.format(similarBase, it)
            mode = 'a' if os.path.exists(similarFile) else 'w'
            with open(similarFile, mode) as f:
                f.write('{}\t{}\n'.format(molSMI, '; '.join(simCids)))

            idsDic = self.getNamesFromCID(cid)
            smiDic = {molSMI: idsDic}

            smisLists[it].append(smiDic)

    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if not smi.lower() == 'smiles':
                    break
        return smi

    def getSMI(self, mol, nt):
        fnSmall = os.path.abspath(mol.getFileName())
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))
        if ext != '.smi':
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))
            args = ' -i "{}" -of smi -o {} --outputDir {} -nt {}'.format(fnSmall, fnOut, outDir, nt)

            if self.useManager.get() == RDKIT or fnSmall.endswith(".mae") or fnSmall.endswith(".maegz"):
                Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir)

            # Formatting with OpenBabel
            elif self.useManager.get() == OPENBABEL:
                Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

        return self.parseSMI(fnOut)

    def getCIDFromSmiles(self, smi):
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT" % smi
        with urlopen(url) as response:
            cid = response.read().decode('utf-8').split()[0]
        return cid

    def getSimilarityListKey(self, smi):
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{}' \
              '/JSON?Threshold=95&MaxRecords=5'.format(smi)
        listKey = ''
        while not listKey:
            with urlopen(url) as response:
                jDic = json.loads(response.read().decode('utf-8'))
                listKey = jDic['Waiting']['ListKey']

        return listKey

    def getSimilarIds(self, listKey):
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{}/cids/txt'.format(listKey)
        cids = []
        while len(cids) == 0:
            time.sleep(2)
            with urlopen(url) as response:
                r = response.read().decode('utf-8')
                if not 'Your request is running' in r:
                    for line in r.split('\n'):
                        if line.strip():
                            cids.append(line.strip())
        return cids

    def getNamesFromCID(self, cid):
        idsDic = {'PubChemID': cid}
        url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{}/synonyms/TXT".format(cid)
        try:
            with urlopen(url) as response:
                r = response.read().decode('utf-8')
                for i, line in enumerate(r.split('\n')):
                    if i == 0:
                        idsDic['PubChemName'] = line.strip()

                    if line.upper().startswith('ZINC'):
                        idsDic['ZINC_ID'] = line.strip()
                    elif line.upper().startswith('CHEMBL'):
                        idsDic['ChEMBL_ID'] = line.strip()
        except:
            print('No synonims found for CID: {}'.format(cid))
            idsDic['PubChemName'], idsDic['ZINC_ID'], idsDic['ChEMBL_ID'] = [None] * 3
        return idsDic

    def _validate(self):
        errors = []
        firstItem = self.inputSet.get().getFirstItem()
        if not hasattr(firstItem, "smallMoleculeFile"):
            errors.append("The input set does not contain small molecules")
        return errors
