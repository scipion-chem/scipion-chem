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
from pwchem.utils import performBatchThreading, makeSubsets, concatThreadFiles
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
        form.addParam('inputSmallMolecules', params.PointerParam, pointerClass="SetOfSmallMolecules",
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
        aSteps = []
        nt = self.numberOfThreads.get()
        inLen = len(self.inputSmallMolecules.get())

        # Ensuring there are no more subsets than input molecules
        nSubsets = min(nt - 1, inLen)
        iStep = self._insertFunctionStep(self.createInputStep, nSubsets)
        for it in range(nSubsets):
            cStep = self._insertFunctionStep(self.convertStep, it, prerequisites=[iStep])
            aSteps += [self._insertFunctionStep(self.identifyStep, it, prerequisites=cStep)]
        self._insertFunctionStep(self.createOutputStep, prerequisites=aSteps)


    def createInputStep(self, nSubsets):
        ligFiles = []
        for mol in self.inputSmallMolecules.get():
            ligFiles.append(os.path.abspath(mol.getPoseFile()) if mol.getPoseFile()
                            else os.path.abspath(mol.getFileName()))

        inputSubsets = makeSubsets(ligFiles, nSubsets, cloneItem=False)
        for it, fileSet in enumerate(inputSubsets):
            with open(self.getInputFile(it), 'w') as f:
                f.write(' '.join(fileSet))

        os.mkdir(self.getOutputSMIDir())
        os.mkdir(self.getSimilarSMIDir())

    def convertStep(self, it):
        inFile = self.getInputFile(it)
        with open(inFile) as f:
            molFiles = f.read().strip().split()

        with open(self.getSMILigandFile(it), 'w') as f:
            for molFile in molFiles:
                smi = self.getSMI(molFile)
                f.write(f'{smi} {molFile}\n')

    def identifyStep(self, it):
        smiList = []
        with open(self.getSMILigandFile(it)) as f:
            for line in f:
                smiList.append(line.split()[0])
        self.identifySMI(smiList, it)

    def createOutputStep(self):
        allSimFile = self._getPath('similarIds.txt')
        if not os.path.exists(allSimFile):
            concatThreadFiles(allSimFile, self.getSimilarSMIDir())

        allIdentFile = self._getPath('identification.txt')
        if not os.path.exists(allIdentFile):
            concatThreadFiles(allIdentFile, self.getOutputSMIDir())

        allInputSMIFile = self._getPath('inputSMI.smi')
        if not os.path.exists(allInputSMIFile):
            concatThreadFiles(allInputSMIFile, self._getExtraPath())

        identDic = self.parseIdentification(allIdentFile)
        molDic = self.parseSMIFile(allInputSMIFile)

        outputSet = self.inputSmallMolecules.get().create(self._getPath())
        for smi in identDic:
            mol = molDic[smi]
            # Setting Ids
            for dbName in DbChoices:
                if dbName in identDic[smi]:
                    setattr(mol, dbName, pwobj.String(identDic[smi][dbName]))
                else:
                    setattr(mol, dbName, pwobj.String(None))

            # Setting molname
            chdbName = self.getEnumText('nameDatabase')
            if chdbName != 'None' and chdbName in identDic[smi] and identDic[smi][chdbName]:
                mol.setMolName(identDic[smi][chdbName])
            outputSet.append(mol)

        if len(outputSet) > 0:
            self._defineOutputs(outputSmallMolecules=outputSet)
            self._defineSourceRelation(self.inputSmallMolecules, outputSet)


    ########################### UTILS FUNCTIONS ##########################
    def getInputDir(self):
      return os.path.abspath(self._getTmpPath())

    def getInputFile(self, it):
      return os.path.join(self.getInputDir(), f'inputLigandFiles_{it}.txt')

    def getInputLigandsDir(self, it):
        return os.path.abspath(self._getExtraPath(f'inputSmallMolecules_{it}'))

    def getSMILigandFile(self, it):
        return os.path.abspath(self._getExtraPath(f'inputSMI_{it}.smi'))

    def getOutputSMIDir(self):
        oDir = os.path.abspath(self._getExtraPath(f'identification'))
        return oDir

    def getOutputSMIFile(self, it):
        return os.path.join(self.getOutputSMIDir(), f'identification_{it}.txt')

    def getSimilarSMIDir(self):
        oDir = os.path.abspath(self._getExtraPath(f'similarIds'))
        return oDir

    def getSimilarSMIFile(self, it):
        return os.path.join(self.getSimilarSMIDir(), f'similarIds_{it}.txt')

    def identifySMI(self, smis, it):
        '''Function to run in performBatchThreading.
        Executes the smi identification and returns a dictionary as:
        {molSMI: {'PubChemID': cid, 'PubChemName': pubChemName, 'ZINC_ID': zincID, 'ChEMBL_ID': chemblID}}
        '''
        for molSMI in smis:
            cid = self.getCIDFromSmiles(molSMI)
            if cid == "0":
                print('No exact match found for SMILES: {}\nLooking for similar compounds'.format(molSMI))
                listKey = self.getSimilarityListKey(molSMI)
                simCids = self.getSimilarIds(listKey) if listKey else []
                cid = simCids[0] if len(simCids) > 0 else None
            else:
                simCids = [cid]

            if simCids:
                simFile = self.getSimilarSMIFile(it)
                mode = 'a' if os.path.exists(simFile) else 'w'
                with open(simFile, mode) as f:
                    simIdsStr = "\t".join(simCids)
                    f.write(f'{molSMI}\t{simIdsStr}\n')

            if cid:
                idsDic = self.getNamesFromCID(cid)
                identFile = self.getOutputSMIFile(it)
                mode = 'a' if os.path.exists(identFile) else 'w'
                with open(identFile, mode) as f:
                    f.write(f'{molSMI}')
                    for dbName, molID in idsDic.items():
                        f.write(f'\t{dbName}\t{molID}')
                    f.write('\n')

    def parseIdentification(self, identFile):
        identDic = {}
        with open(identFile) as f:
            for line in f:
                sline = line.split()
                smi = sline[0]
                identDic[smi] = {}

                for i in range(1, len(sline)-1, 2):
                    identDic[smi][sline[i]] = sline[i+1]
        return identDic

    def parseSMIFile(self, smiFile):
        inputDic = {}
        for mol in self.inputSmallMolecules.get():
            inputDic[os.path.abspath(mol.getFileName())] = mol.clone()
            if mol.getPoseFile():
                inputDic[os.path.abspath(mol.getPoseFile())] = mol.clone()

        molDic = {}
        with open(smiFile) as f:
            for line in f:
                sline = line.split()
                molDic[sline[0]] = inputDic[sline[1]]
        return molDic

    def parseSMI(self, smiFile):
        smi = None
        with open(smiFile) as f:
            for line in f:
                smi = line.split()[0].strip()
                if not smi.lower() == 'smiles':
                    break
        return smi

    def getSMI(self, fnSmall):
        fnRoot, ext = os.path.splitext(os.path.basename(fnSmall))
        if ext != '.smi':
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.abspath(self._getExtraPath(fnRoot + '.smi'))
            args = ' -i "{}" -of smi -o {} --outputDir {}'.format(fnSmall, fnOut, outDir)

            if self.useManager.get() == RDKIT or fnSmall.endswith(".mae") or fnSmall.endswith(".maegz"):
                Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir)

            # Formatting with OpenBabel
            elif self.useManager.get() == OPENBABEL:
                Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

        return self.parseSMI(fnOut)

    def getCIDFromSmiles(self, smi):
        i = 0
        while i < 5:
            try:
                url = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/%s/cids/TXT" % smi
                with urlopen(url) as response:
                    cid = response.read().decode('utf-8').split()[0]
                i = 5
            except:
                i += 1
                time.sleep(1)
                cid = "0"

        return cid

    def getSimilarityListKey(self, smi):
        url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/similarity/smiles/{}' \
              '/JSON?Threshold=95&MaxRecords=5'.format(smi)
        listKey, i = '', 0
        while i < 5:
            try:
                with urlopen(url) as response:
                    jDic = json.loads(response.read().decode('utf-8'))
                    listKey = jDic['Waiting']['ListKey']
                i = 5
            except:
                i += 1

        return listKey

    def getSimilarIds(self, listKey):
        url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/listkey/{listKey}/cids/txt'
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
        firstItem = self.inputSmallMolecules.get().getFirstItem()
        if not hasattr(firstItem, "smallMoleculeFile"):
            errors.append("The input set does not contain small molecules")
        return errors
