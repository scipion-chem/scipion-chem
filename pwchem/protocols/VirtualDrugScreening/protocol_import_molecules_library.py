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

import os, requests, re, subprocess, shutil, time, gzip

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SmallMoleculesLibrary
from pwchem.utils import getBaseFileName, runInParallel, reduceNRandomLines, swapColumns, gunzipFile

baseZINCUrl = "https://files.docking.org/zinc22/2d-all"
pubchemUrl = 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz'

libraries = ['ZINC22', 'PubChem']
ZINC, PUBCHEM = 0, 1

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
      file_url = f"{url}{filename}"
      urls.append(file_url)

  return urls

def downloadUrlFile(url, oDir):
  '''Downloads a file from its url and return the local file downloaded'''
  coms = f"wget -P {oDir} {url}"
  subprocess.run(coms, shell=True)
  return os.path.join(oDir, getBaseFileName(url))


class ProtChemImportMoleculesLibrary(EMProtocol):
  """Import a small molecules library from local or a web server in a single SMI format file.
  """
  _label = 'import small molecules library'
  stepsExecutionMode = params.STEPS_PARALLEL

  zincSubsets = {'big-leads (13b)': ((23, 25), (-5.00, 3.40)),
                 'big-n-greasy (17b)': ((26, 30), (3.50, 9.00)),
                 'drug-like (67.7b)': ((4, 28), (-5.00, 4.90)),
                 'fragments (1.0b)': ((16, 20), (-5.00, 3.40)),
                 'fragments2 (70.0m)': ((8, 16), (-5.00, 2.40)),
                 'goldilocks (5.5b)': ((19, 24), (1.10, 2.90)),
                 'lead-like (16.1b)': ((17, 25), (-5.00, 3.40)),
                 'lugs (22.8b)': ((24, 26), (-5.00, 3.90)),
                 'medium-leads (2.6b)': ((20, 22), (-5.00, 3.40)),
                 'shards (177.2k)': ((5, 10), (-5.00, 2.40)),
                 'small-leads (522m)': ((17, 19), (-5.00, 3.40))}

  def _defineParams(self, form):
    form.addSection(label='Input')

    form.addParam('defLibraries', params.BooleanParam, default=False, label='Download a predefined library: ',
                  help='Download a predefined library from a website')
    group = form.addGroup('Default libraries', condition='defLibraries')
    group.addParam('choicesLibraries', params.EnumParam, default=ZINC,
                   label='Download ligand library: ', condition='defLibraries', choices=libraries,
                   help='Choose the predefined library you want to download.\n'
                        'ZINC: SMI libraries from ZINC22'
                        'PubChem: SMI libraries from PubChem')

    group.addParam('zincSubset', params.EnumParam, default=0, choices=['Any', 'Manual'] + list(self.zincSubsets.keys()),
                   label='ZINC subset: ', condition='defLibraries and choicesLibraries == 0',
                   help='ZINC define subsets based on Heavy atom count and logP ranges.'
                        'Use the wizard to select a predefined subset and the ranges will automatically update.'
                        'You can also define your own ranges.\n'
                        'Be aware that for subsets other than "Any" and "lead-like", ZINC is not prepared for the '
                        'random download, so the whole dataset will be downloaded and then the selected number of '
                        'molecules will be randomly picked.')

    group.addParam('setRanges', params.LabelParam, label='Set subset ranges: ',
                   condition='defLibraries and choicesLibraries == 0 and zincSubset>1',
                   help='Set the ranges for the selected subset')

    line = group.addLine('Heavy atoms range: ',
                         condition='defLibraries and choicesLibraries == 0 and zincSubset!=0',
                         help='Range of heavy atoms for the downloaded molecules, limits included.')
    line.addParam('minSize', params.IntParam, label='Min: ', default=4)
    line.addParam('maxSize', params.IntParam, label='Max: ', default=20)

    line = group.addLine('LogP range: ',
                         condition='defLibraries and choicesLibraries == 0 and zincSubset!=0',
                         help='Range of logP for the downloaded molecules, limits included.')
    line.addParam('minlogP', params.FloatParam, label='Min: ', default=0)
    line.addParam('maxlogP', params.FloatParam, label='Max: ', default=2)

    group.addParam('nMols', params.IntParam, default=10000, label='Number of desired molecules: ',
                   condition='defLibraries',
                   help='Number of desired random molecules from the Pubchem/ZINC22 to download.'
                        'Be aware that for PubChem and ZINC subsets other than "Any" and "lead-like", the protocol is '
                        'not prepared for the random download, so the whole dataset will be downloaded and then '
                        'selected number of molecules will be randomly picked. This can take a while.')

    group = form.addGroup('Local files', condition='not defLibraries')
    group.addParam('filePath', params.PathParam, condition='not defLibraries',
                   label='Library file: ', help='File path to the SMI library in local.')
    group.addParam('headers', params.StringParam, condition='not defLibraries',
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
    if libChoice == ZINC:
      if self.checkAvailableRandom():
        subset = '' if self.getEnumText('zincSubset') == 'Any' else 'lead-like'
        self.getZINC22Random(nMols, subset, oFile)

      else:
        sizeRange = [self.minSize.get(), self.maxSize.get()]
        logPRange = [self.minlogP.get(), self.maxlogP.get()]
        self.getZINC22Ranges(sizeRange, logPRange, oFile)
        if nMols > 0:
          reduceNRandomLines(oFile, nMols, oDir=self._getTmpPath())

    elif libChoice == PUBCHEM:
      gzFile = downloadUrlFile(pubchemUrl, self._getPath())
      smiFile = gunzipFile(gzFile, oFile=self.getOutLibraryFile(), remove=True)
      if nMols > 0:
        reduceNRandomLines(smiFile, nMols, oDir=self._getTmpPath())
      swapColumns(smiFile, oDir=self._getTmpPath())

  def createOutputStep(self):
    oFile = self.getOutLibraryFile()
    if not self.defLibraries.get():
      shutil.copy(self.filePath.get(), oFile)

    headers = self.getHeaders()
    outputLib = SmallMoleculesLibrary(libraryFilename=oFile, headers=headers)
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

    task_id = response.json().get("task")
    if not task_id:
      print("Error: No task ID found in the response.")

    status_url = f"https://cartblanche22.docking.org/substance/random/{task_id}.txt"
    while True:
      status_response = requests.get(status_url)
      if status_response.status_code != 200:
        print(f"Error: Failed to check status. Status code: {status_response.status_code}")

      status_json = status_response.json()
      status = status_json.get("status")

      if status == "SUCCESS":
        results = status_json.get("result")
        self.writeSMIFile(results, outFile)
        break

      elif status == "FAILURE":
        print("Task failed.")

      else:
        time.sleep(3)

  def getZINC22Ranges(self, sizeRange, logPRange, outFile):
    nt = self.numberOfThreads.get()
    sizeStrs = [f'{size:02}' for size in range(sizeRange[0], sizeRange[1]+1)]
    allUrls = runInParallel(getMatchingFiles, logPRange, paramList=sizeStrs, jobs=nt)
    allUrls = [url for urlList in allUrls for url in urlList ]

    oDir = os.path.abspath(self._getTmpPath())
    pFiles = runInParallel(downloadUrlFile, oDir, paramList=allUrls, jobs=nt)
    self.mergeZINCFiles(pFiles, outFile)


  #################### UTILS FUNCTIONS #################

  def writeSMIFile(self, zincResults, outFile):
    if zincResults:
      with open(outFile, 'w') as f:
        for line in zincResults.split('\n')[1:]:
          if line:
            sline = line.split()
            hAtoms, logP = self.getInfoFromTranche(sline[0])
            f.write(f'{sline[2]}\t{sline[1]}\t{hAtoms}\t{logP}\n')

  def checkAvailableRandom(self):
    sizeRange, logpRange = (self.minSize.get(), self.maxSize.get()), (self.minlogP.get(), self.maxlogP.get())
    leadSizeRange, leadLogPRange = self.zincSubsets['lead-like (16.1b)']
    return self.getEnumText('zincSubset') == 'Any' or \
           (sizeRange[0] == leadSizeRange[0] and sizeRange[1] == leadSizeRange[1] and
            logpRange[0] == leadLogPRange[0] and logpRange[1] == leadLogPRange[1] and self.nMols.get() > 0)

  def mergeZINCFiles(self, zFiles, oFile):
    '''Merge a set of smi.gz files into a smi file, adding the heavy atoms count and logP'''
    with open(oFile, 'w') as f:
      for zFile in zFiles:
        tranche = os.path.split(zFile)[-1].split('.')[0]
        hAtoms, logP = self.getInfoFromTranche(tranche)
        with gzip.open(zFile, 'rb') as fIn:
          for line in fIn:
            f.write(f'{line.decode().strip()}\t{hAtoms}\t{logP}\n')


############# VALIDATION FUNCTIONS ##################

  def _validate(self):
    errors = []
    return errors

  def getOutLibraryFile(self):
    if self.defLibraries.get():
      if self.getEnumText('choicesLibraries') == 'ZINC22':
        base = 'zinc22.smi'
      else:
        base = 'pubchem.smi'
    else:
        base = getBaseFileName(self.filePath.get())
    return os.path.abspath(self._getPath(base))

  def getInfoFromTranche(self, tranche):
    '''Return the heavy atoms and the logP of a ZINC tranch like H12P200'''
    hAtoms = int(tranche[1:3])
    sign = 1 if tranche[3] == 'P' else -1
    logP = int(tranche[4:]) * sign/100
    return hAtoms, logP

  def getHeaders(self):
    headers = ['SMI', 'Name']
    if self.defLibraries.get() and self.choicesLibraries.get() == ZINC:
        headers += ['Heavy_atoms', 'LogP']
    elif not self.defLibraries.get() and self.headers.get().strip():
        headers = [h.strip() for h in self.headers.get().split(',')]
    return headers



