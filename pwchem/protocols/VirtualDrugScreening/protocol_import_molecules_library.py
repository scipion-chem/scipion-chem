# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Alberto M. Parra Pérez (amparraperez@gmail.com)
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

import os, requests, re, shutil, time, gzip
from itertools import product

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SmallMoleculesLibrary
from pwchem.utils import getBaseFileName, runInParallel, reduceNRandomLines, swapColumns, gunzipFile, flipDic, \
  downloadUrlFile

MIN_STR, MAX_STR = 'Min: ', 'Max: '
NO_DEF_LIB = 'not defLibraries'

baseZINCUrl = "https://files.docking.org/zinc22/2d-all"
pubchemUrl = 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz'

libraries = ['ZINC20', 'ZINC22', 'PubChem']
ZINC20, ZINC22, PUBCHEM = 0, 1, 2

def getMatchingFiles(hSize, logPRange):
  url = f"{baseZINCUrl}/H{hSize}/"
  try:
    response = requests.get(url)
    response.raise_for_status()
  except Exception as e:
    print(f"Failed to fetch {url}: {e}")
    return []

  pattern = re.compile(rf"H{hSize}([PM]\d+)\.smi\.gz")
  matches = list(set(pattern.findall(response.text)))

  urls = []
  for match in matches:
    let = match[0]
    sign = 1 if let == 'P' else -1
    num = int(match[1:]) * sign / 100

    if logPRange[0] <= num <= logPRange[1]:
      filename = f"H{hSize}{match}.smi.gz"
      fileUrl = f"{url}{filename}"
      urls.append(fileUrl)

  return urls

class ProtChemImportMoleculesLibrary(EMProtocol):
  """Import a small molecules library from local or a web server in a single SMI format file.
  """
  _label = 'import small molecules library'
  stepsExecutionMode = params.STEPS_PARALLEL

  zinc22Subsets = {'big-leads (13B)': ((23, 25), (-5.00, 3.40)),
                   'big-n-greasy (17B)': ((26, 30), (3.50, 9.00)),
                   'drug-like (67.7B)': ((4, 28), (-5.00, 4.90)),
                   'fragments (1.0B)': ((16, 20), (-5.00, 3.40)),
                   'fragments2 (70.0M)': ((8, 16), (-5.00, 2.40)),
                   'goldilocks (5.5B)': ((19, 24), (1.10, 2.90)),
                   'lead-like (16.1B)': ((17, 25), (-5.00, 3.40)),
                   'lugs (22.8B)': ((24, 26), (-5.00, 3.90)),
                   'medium-leads (2.6B)': ((20, 22), (-5.00, 3.40)),
                   'shards (177.2K)': ((5, 10), (-5.00, 2.40)),
                   'small-leads (522M)': ((17, 19), (-5.00, 3.40))}

  zinc20Subsets = {'big-n-greasy (40M)': ((500, ">500"), (4.5, '>5')),
                   'drug-like (1.87B)': ((250, 500), (-1, 5)),
                   'fragments (148M)': ((200, 250), (-1, 3.5)),
                   'fragments2 (684M)': ((250, 325), (-1, 3.5)),
                   'goldilocks (590M)': ((300, 350), (2, 3)),
                   'lead-like (907M)': ((300, 350), (-1, 3.5)),
                   'lugs (969M)': ((350, 450), (-1, 4.5)),
                   'shards (17M)': ((200, 200), (-1, '>5'))}

  zinc20SizeTranches = {200: 'A', 250: 'B', 300: 'C', 325: 'D', 350: 'E', 375: 'F', 400: 'G', 425: 'H', 450: 'I',
                        500: 'J', '>500': 'K'}
  zinc20LogPTranches = {-1: 'A', 0: 'B', 1: 'C', 2: 'D', 2.5: 'E', 3: 'F', 3.5: 'G', 4: 'H', 4.5: 'I',
                        5: 'J', '>5': 'K'}

  # Reactivity groups for Zinc20, 3rd letter in code. First element in tuple is the exclusive, the second not
  reactGroups = {'Anodyne': ('A', 'A'), 'Bother': ('B', 'AB'), 'Clean': ('C', 'ABC'), 'Standard': ('E', 'ABCE'),
                 'Reactive': ('G', 'ABCEG'), 'Hot': ('I', 'ABCEGI')}
  # Purchasability groups for Zinc20, 4th letter in code. First element in tuple is the exclusive, the second not
  purchGroups = {'In-Stock': ('AB', 'AB'), 'Agent': ('C', 'ABC'), 'Wait OK': ('D', 'ABCD'), 'Boutique': ('E', 'ABCDE'),
                 'Annotated': ('F', 'ABCDEF')}

  def _defineZINC20TranchesParams(self, form, condition='True'):
    form.addParam('zinc20Subset', params.EnumParam, default=0,
                  choices=['Any', 'Manual'] + list(self.zinc20Subsets.keys()),
                  label='ZINC20 subset: ', condition=condition,
                  help='ZINC20 define subsets based on Heavy atom count and logP ranges.'
                       'Use the wizard to select a predefined subset and the ranges will automatically update.'
                       'You can also define your own ranges.\n'
                       'Be aware that for subsets other than "Any" and "lead-like", ZINC22 is not prepared for the '
                       'random download, so the whole dataset will be downloaded and then the selected number of '
                       'molecules will be randomly picked.')

    # Size / logP ranges ZINC20
    form.addParam('setRanges20', params.LabelParam, label='Set ZINC20 subset ranges: ',
                  condition=f'{condition} and zinc20Subset>1',
                  help='Set the ranges for the selected subset')

    line = form.addLine('Size range (Daltons): ', condition=f'{condition} and zinc20Subset!=0',
                        help='Range of size for the downloaded molecules in Daltons, limits included.\n'
                             'If Manual, the range goes from 200 to >500')
    line.addParam('minSize20', params.StringParam, label=MIN_STR, default=250)
    line.addParam('maxSize20', params.StringParam, label=MAX_STR, default=500)

    line = form.addLine('LogP range ZINC20: ', condition=f'{condition} and zinc20Subset!=0',
                        help='Range of logP for the downloaded molecules, limits included.\n'
                             'If Manual, the range goes from -1 to >5')
    line.addParam('minLogP20', params.StringParam, label=MIN_STR, default=-1)
    line.addParam('maxLogP20', params.StringParam, label=MAX_STR, default=5)

    line = form.addLine('Reactivity: ', condition=condition,
                        help='Reactivity subsets defined in ZINC20. In the order showed, each subset includes the previous '
                       'ones unless the exclusive param is set to true, then the molecules exclusively in that group '
                       'are selected.\nFor more information visit ZINC20: https://zinc20.docking.org/ ')
    line.addParam('reactivity', params.EnumParam, default=0, choices=list(self.reactGroups.keys()),
                  label='Reactivity subset: ')
    line.addParam('reactExclusive', params.BooleanParam, default=False, label='Exclusive: ')

    line = form.addLine('Purchasability: ', condition=condition,
                        help='Purchasability subsets defined in ZINC20. In the order showed, each subset includes the '
                       'previous ones unless the exclusive param is set to true, then the molecules exclusively in '
                       'that group are selected.\nFor more information visit ZINC20: https://zinc20.docking.org/ ')
    line.addParam('purchasability', params.EnumParam, default=0, choices=list(self.purchGroups.keys()),
                  label='Purchasability subset: ')
    line.addParam('purchExclusive', params.BooleanParam, default=False, label='Exclusive: ')

    return form


  def _defineParams(self, form):
    form.addSection(label='Input')

    form.addParam('defLibraries', params.BooleanParam, default=False, label='Download a predefined library: ',
                  help='Download a predefined library from a website')
    form.addParam('nMols', params.IntParam, default=10000, label='Number of desired molecules: ',
                  condition='defLibraries',
                  help='Number of desired random molecules from the Pubchem/ZINC to download.'
                       'Be aware that for PubChem and ZINC22 subsets other than "Any" and "lead-like", the protocol is '
                       'not prepared for the random download, so the whole dataset will be downloaded and then '
                       'selected number of molecules will be randomly picked. This can take a while.')
    form.addParam('randomSeed', params.IntParam, default=44, label='Random seed: ', condition='nMols>0',
                  expertLevel=params.LEVEL_ADVANCED, help='Random seed for the random molecules')

    group = form.addGroup('Default libraries', condition='defLibraries')
    group.addParam('choicesLibraries', params.EnumParam, default=ZINC20,
                   label='Download ligand library: ', condition='defLibraries', choices=libraries,
                   help='Choose the predefined library you want to download.\n'
                        'ZINC20: SMI libraries from ZINC20'
                        'ZINC22: SMI libraries from ZINC22'
                        'PubChem: SMI libraries from PubChem')
    groupT = form.addGroup('ZINC20 tranches', condition=f'defLibraries and choicesLibraries == {ZINC20}')
    groupT = self._defineZINC20TranchesParams(groupT, condition=f'defLibraries and choicesLibraries == {ZINC20}')
    group.addParam('zinc22Subset', params.EnumParam, default=0,
                   condition=f'defLibraries and choicesLibraries == {ZINC22}',
                   label='ZINC22 subset: ', choices=['Any', 'Manual'] + list(self.zinc22Subsets.keys()),
                   help='ZINC22 define subsets based on Heavy atom count and logP ranges.'
                        'Use the wizard to select a predefined subset and the ranges will automatically update.'
                        'You can also define your own ranges.\n'
                        'Be aware that for subsets other than "Any" and "lead-like", ZINC22 is not prepared for the '
                        'random download, so the whole dataset will be downloaded and then the selected number of '
                        'molecules will be randomly picked.')

    # Size / logP ranges ZINC22
    group.addParam('setRanges22', params.LabelParam, label='Set ZINC22 subset ranges: ',
                   condition=f'defLibraries and choicesLibraries == {ZINC22} and zinc22Subset>1',
                   help='Set the ranges for the selected subset')

    line = group.addLine('Heavy atoms range: ',
                         condition=f'defLibraries and choicesLibraries == {ZINC22} and zinc22Subset!=0',
                         help='Range of heavy atoms for the downloaded molecules, limits included.')
    line.addParam('minSize22', params.IntParam, label=MIN_STR, default=4)
    line.addParam('maxSize22', params.IntParam, label=MAX_STR, default=20)

    line = group.addLine('LogP range ZINC22: ',
                         condition=f'defLibraries and choicesLibraries == {ZINC22} and zinc22Subset!=0',
                         help='Range of logP for the downloaded molecules, limits included.')
    line.addParam('minLogP22', params.FloatParam, label=MIN_STR, default=0)
    line.addParam('maxLogP22', params.FloatParam, label=MAX_STR, default=2)


    group = form.addGroup('Local files', condition=NO_DEF_LIB)
    group.addParam('filePath', params.PathParam, condition=NO_DEF_LIB,
                   label='Library file: ', help='File path to the SMI library in local.')
    group.addParam('headers', params.StringParam, condition=NO_DEF_LIB,
                   label='Library headers: ', default='',
                   help='Headers of the selected library file. Set them separated by commas. e.g: SMI,Name')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    if self.defLibraries:
      self._insertFunctionStep('downloadStep')
    self._insertFunctionStep('createOutputStep')

  def downloadStep(self):
    libChoice = self.choicesLibraries.get()
    oFile = self.getOutLibraryFile()
    nMols = self.nMols.get()

    if libChoice == ZINC20:
      self.getZINC20Range(oFile)
      if nMols > 0:
        reduceNRandomLines(oFile, nMols, oDir=self._getTmpPath(), seed=self.randomSeed.get())

    elif libChoice == ZINC22:
      if self.checkAvailableRandom():
        subset = '' if self.getEnumText('zinc22Subset') == 'Any' else 'lead-like'
        self.getZINC22Random(nMols, subset, oFile)

      else:
        self.getZINC22Ranges(oFile)
        if nMols > 0:
          reduceNRandomLines(oFile, nMols, oDir=self._getTmpPath(), seed=self.randomSeed.get())

    elif libChoice == PUBCHEM:
      gzFile = downloadUrlFile(pubchemUrl, self._getPath())
      smiFile = gunzipFile(gzFile, oFile=self.getOutLibraryFile(), remove=True)
      if nMols > 0:
        reduceNRandomLines(smiFile, nMols, oDir=self._getTmpPath(), seed=self.randomSeed.get())
      swapColumns(smiFile, oDir=self._getTmpPath())

  def createOutputStep(self):
    oFile = self.getOutLibraryFile()
    if not self.defLibraries.get():
      shutil.copy(self.filePath.get(), oFile)

    headers = self.getHeaders()
    outputLib = SmallMoleculesLibrary(libraryFilename=oFile, headers=headers)
    outputLib.calculateLength()
    self._defineOutputs(outputLibrary=outputLib)


  ################# MAIN FUNCTIONS #####################

  def getZINC22Random(self, count=1000, subset='', outFile='output.smi'):
    url = "https://cartblanche22.docking.org/substance/random.txt"
    data = {'count': count}
    if subset:
      data.update({'subset': subset})

    response = requests.post(url, data=data)
    if response.status_code != 200:
      print(f"Error: Failed to submit job. Status code: {response.status_code}")

    taskId = response.json().get("task")
    if not taskId:
      print("Error: No task ID found in the response.")

    statusUrl = f"https://cartblanche22.docking.org/substance/random/{taskId}.txt"
    while True:
      statusResponse = requests.get(statusUrl)
      if statusResponse.status_code != 200:
        print(f"Error: Failed to check status. Status code: {statusResponse.status_code}")

      statusJson = statusResponse.json()
      status = statusJson.get("status")

      if status == "SUCCESS":
        results = statusJson.get("result")
        self.writeSMIFile(results, outFile)
        break

      elif status == "FAILURE":
        print("Task failed.")

      else:
        time.sleep(3)

  def getZINC22Ranges(self, outFile):
    sizeRange = self.checkZINC22Ranges([self.minSize22.get(), self.maxSize22.get()])
    logPRange = [self.minLogP22.get(), self.maxLogP22.get()]

    nt = self.numberOfThreads.get()
    sizeStrs = [f'{size:02}' for size in range(sizeRange[0], sizeRange[1]+1)]
    allUrls = runInParallel(getMatchingFiles, logPRange, paramList=sizeStrs, jobs=nt)
    allUrls = [url for urlList in allUrls for url in urlList]

    oDir = os.path.abspath(self._getTmpPath())
    pFiles = runInParallel(downloadUrlFile, oDir, paramList=allUrls, jobs=nt)
    self.mergeZINC22Files(pFiles, outFile)

  def getZINC20Codes(self):
    sizeCodes = logPCodes = [letter for letter in 'ABCDEFGHIJK']
    if self.zinc20Subset.get() != 0:
      # Reducing the size and logP tranches to the ones defined by user
      sizeTranches, logpTranches = list(self.zinc20SizeTranches.keys()), list(self.zinc20LogPTranches.keys())
      maxSize = self.maxSize20.get() if self.maxSize20.get() != sizeTranches[-1] else sizeTranches[-2] + 1
      maxLogP = self.maxLogP20.get() if self.maxLogP20.get() != logpTranches[-1] else logpTranches[-2] + 1
      minSizeIdx, maxSizeIdx = self.getCloserTrancheIdx(sizeTranches, self.minSize20.get()), \
                               self.getCloserTrancheIdx(sizeTranches, maxSize)
      minLogPIdx, maxLogPIdx = self.getCloserTrancheIdx(logpTranches, self.minLogP20.get()), \
                               self.getCloserTrancheIdx(logpTranches, maxLogP)

      sizeCodes, logPCodes = sizeCodes[minSizeIdx: maxSizeIdx + 1], logPCodes[minLogPIdx: maxLogPIdx + 1]

    reactExclIdx = 0 if self.reactExclusive.get() else 1
    reactCodes = self.reactGroups[self.getEnumText('reactivity')][reactExclIdx]

    purchExclIdx = 0 if self.purchExclusive.get() else 1
    purchCodes = self.purchGroups[self.getEnumText('purchasability')][purchExclIdx]
    return sizeCodes, logPCodes, reactCodes, purchCodes

  def getZINC20Range(self, oFile):
    sizeCodes, logPCodes, reactCodes, purchCodes = self.getZINC20Codes()
    downTranches = [''.join(p) for p in product(sizeCodes, logPCodes, reactCodes, purchCodes)]
    allUrls = [f'https://files.docking.org/2D/{tranche[:2]}/{tranche}.smi' for tranche in downTranches]

    oDir = os.path.abspath(self._getTmpPath())
    nt = self.numberOfThreads.get()
    pFiles = runInParallel(downloadUrlFile, oDir, paramList=allUrls, jobs=nt)
    pFiles = [pFile for pFile in pFiles if os.path.exists(pFile)]
    self.mergeZINC20Files(pFiles, oFile)

  #################### UTILS FUNCTIONS #################

  def checkZINC22Ranges(self, sizeRange):
    sr = []
    for ha in sizeRange:
      if ha < 4:
        sr.append(4)
      elif ha > 49:
        sr.append(49)
      else:
        sr.append(ha)
    return sr

  def getCloserTrancheIdx(self, tranches, value):
    tr, i = tranches[0], 1
    while i < len(tranches) and float(value) > tr:
      tr = tranches[i]
      i += 1
    return i - 1

  def writeSMIFile(self, zincResults, outFile):
    if zincResults:
      with open(outFile, 'w') as f:
        for line in zincResults.split('\n')[1:]:
          if line:
            sline = line.split()
            hAtoms, logP = self.getInfoFromTranche22(sline[0])
            f.write(f'{sline[2]}\t{sline[1]}\t{hAtoms}\t{logP}\n')

  def checkAvailableRandom(self):
    sizeRange, logpRange = (self.minSize22.get(), self.maxSize22.get()), (self.minLogP22.get(), self.maxLogP22.get())
    leadSizeRange, leadLogPRange = self.zinc22Subsets['lead-like (16.1B)']
    return self.getEnumText('zinc22Subset') == 'Any' or \
           (sizeRange[0] == leadSizeRange[0] and sizeRange[1] == leadSizeRange[1] and
            logpRange[0] == leadLogPRange[0] and logpRange[1] == leadLogPRange[1] and
            self.nMols.get() > 0)

  def getDecodeZinc20Dics(self):
    sizeDic, logPDic = flipDic(self.zinc20SizeTranches), flipDic(self.zinc20LogPTranches)
    sizeDic['K'], logPDic['K'] = 600, 6
    return [sizeDic,  # Size (Daltons)
            logPDic,  # LogP
            {v[0]: k for k, v in self.reactGroups.items()},  # Reactivity
            {v1: k for k, v in self.purchGroups.items() for v1 in v[0]}  # Purchasability
            ]

  def decodeTranche(self, tranche, infoDics):
    '''Returns the info from the tranche code as: [size, logP, Reactivity, Purchasability]'''
    return [infoDics[i][tranche[i]] for i in range(len(infoDics))]

  def mergeZINC20Files(self, inFiles, oFile):
    '''Merge a set of smi ZINC20 files into a smi file, adding the size, logP, reactivity and purchasability
    '''
    # Tranche code to values dic
    infoDics = self.getDecodeZinc20Dics()

    with open(oFile, 'w') as f:
      for inFile in inFiles:
        tranche = os.path.split(inFile)[-1].split('.')[0]
        size, logP, react, purch = self.decodeTranche(tranche, infoDics)
        with open(inFile) as fIn:
          fIn.readline()
          for line in fIn:
            f.write(f'{line.strip()}\t{size}\t{logP}\t{react}\t{purch}\n')
  
  def mergeZINC22Files(self, zFiles, oFile):
    '''Merge a set of smi.gz files into a smi file, adding the heavy atoms count and logP'''
    with open(oFile, 'w') as f:
      for zFile in zFiles:
        tranche = os.path.split(zFile)[-1].split('.')[0]
        hAtoms, logP = self.getInfoFromTranche22(tranche)
        with gzip.open(zFile, 'rb') as fIn:
          for line in fIn:
            f.write(f'{line.decode().strip()}\t{hAtoms}\t{logP}\n')


############# VALIDATION FUNCTIONS ##################

  def _validate(self):
    errors = []
    return errors

  def getOutLibraryFile(self):
    if self.defLibraries.get():
      if self.getEnumText('choicesLibraries') == 'ZINC20':
        base = 'zinc20.smi'
      elif self.getEnumText('choicesLibraries') == 'ZINC22':
        base = 'zinc22.smi'
      else:
        base = 'pubchem.smi'
    else:
        base = getBaseFileName(self.filePath.get())
    return os.path.abspath(self._getPath(base))

  def getInfoFromTranche22(self, tranche):
    '''Return the heavy atoms and the logP of a ZINC22 tranch like H12P200'''
    hAtoms = int(tranche[1:3])
    sign = 1 if tranche[3] == 'P' else -1
    logP = int(tranche[4:]) * sign/100
    return hAtoms, logP

  def getHeaders(self):
    headers = ['SMI', 'Name']
    if self.defLibraries.get():
      if self.choicesLibraries.get() == ZINC22:
        headers += ['Heavy_atoms', 'LogP']
      elif self.choicesLibraries.get() == ZINC20:
        headers += ['Size(Daltons)', 'LogP', 'Reactivity', 'Purchasability']
    elif not self.defLibraries.get() and self.headers.get().strip():
        headers = [h.strip() for h in self.headers.get().split(',')]
    return headers



