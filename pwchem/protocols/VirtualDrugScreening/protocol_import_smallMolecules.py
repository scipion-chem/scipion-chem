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

import os, requests, glob, re, json, gzip
from html.parser import HTMLParser
from urllib.parse import urljoin
import random as rd

from pwem.protocols import EMProtocol
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam, LEVEL_ADVANCED, EnumParam, \
  STEPS_PARALLEL, IntParam

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin
from pwchem.utils import performBatchThreading, splitFile, reduceNRandomLines, runInParallel, swapColumns, makeSubsets
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1
DEFAULT_FORMAT = 'sdf'

libraries = ['ECBL', 'ZINC', 'PubChem']

urlECBLDic = {
  'Bioactive (~5000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/Pilot_08_09_2021.sdf',
  'Diversity (~10⁵)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/ECBL_08_09_2021.sdf',
  'Fragments (~1000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Fragments.sdf'}

urlZINCJsonDic = {'FDA (~1600)': 'https://zinc.docking.org/substances/subsets/fda.json?count=all',
                  'World Not FDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.json?count=all',
                  'World (~3600)': 'https://zinc.docking.org/substances/subsets/world.json?count=all',
                  'Investigational only (~3900)': 'https://zinc.docking.org/substances/subsets/investigational-only.json?count=all',
                  # 'In trials (~10k)': 'https://zinc.docking.org/substances/subsets/in-trials.json?count=all',
                  # 'In man (~100k)': 'https://zinc.docking.org/substances/subsets/in-man.json?count=all',
                  # 'In man only (~90k)': 'https://zinc.docking.org/substances/subsets/in-man-only.json?count=all',
                  # 'In vivo (~115k)': 'https://zinc.docking.org/substances/subsets/in-vivo.json?count=all',
                  # 'In vivo only (~16k)': 'https://zinc.docking.org/substances/subsets/in-vivo-only.json?count=all',
                  # 'In cells (~115k)': 'https://zinc.docking.org/substances/subsets/in-cells.json?count=all',
                  # 'In cells only (~129)': 'https://zinc.docking.org/substances/subsets/in-cells-only.json?count=all',
                  # 'In vitro (~276k)': 'https://zinc.docking.org/substances/subsets/in-vitro.json?count=all',
                  # 'In vitro only (~162k)': 'https://zinc.docking.org/substances/subsets/in-vitro-only.json?count=all',
                  # 'Endogenous (~51k)': 'https://zinc.docking.org/substances/subsets/endogenous.json?count=all',
                  # 'Metabolites (~53k)': 'https://zinc.docking.org/substances/subsets/metabolites.json?count=all',
                  'Non-human metabolites (~2k)': 'https://zinc.docking.org/substances/subsets/nonhuman-metabolites.json?count=all',
                  # 'Natural products (~81k)': 'https://zinc.docking.org/substances/subsets/natural-products.json?count=all',
                  # 'Biogenic (~136k)': 'https://zinc.docking.org/substances/subsets/biogenic.json?count=all'
                  }

urlPubChemDic = {'Compound_3D': 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/',
                 'All': 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz'}


class LinkParser(HTMLParser):
  def __init__(self, baseUrl):
    super().__init__()
    self.links = []
    self.baseUrl = baseUrl

  def handle_starttag(self, tag, attrs):
    if tag == 'a':
      attrs = dict(attrs)
      if 'href' in attrs:
        href = attrs['href']
        if not href.startswith('?'):
          fullUrl = urljoin(self.baseUrl, href)
          self.links.append(fullUrl)

def extractLinks(html_content, baseUrl, pattern=r'\.sdf\.gz$'):
  parser = LinkParser(baseUrl)
  parser.feed(html_content)

  # Filtrar los enlaces según el patrón
  filteredLinks = [link for link in parser.links if re.search(pattern, link)]
  return filteredLinks

def downloadSDFZINCBatch(zIds, oFileBase, oFormat='sdf'):
  it, zIds = zIds
  url = 'https://zinc20.docking.org/substances/resolved/'
  data = {
    'paste': "\n".join(zIds), 'output_format': oFormat,
    'identifiers': 'y', 'structures': 'y', 'names': 'y',
  }

  response = requests.post(url, data=data)
  oFile = f'{oFileBase}_{it}.{oFormat}'
  with open(oFile, 'wb') as f:
    f.write(response.content)

  return splitZINCSDF(oFile)

def splitZINCSDF(sdfFile, oDir=None, remove=True):
  oFiles = []
  if oDir is None:
    oDir = os.path.dirname(sdfFile)

  with open(sdfFile) as f:
    sdfText = f.read().strip()

  zincPattern = r'<zinc_id>\s*(\(\d+\))\s*\n(.*)\n'
  for molText in sdfText.split('$$$$'):
    molText = molText.strip()
    if molText:
      matches = re.search(zincPattern, molText)
      if matches:
        number, zincId = matches.group(1), matches.group(2)
        molText = molText.replace(number, '(1)')

        oFile = os.path.join(oDir, f'{zincId}.sdf')
        with open(oFile, 'w') as f:
          f.write(molText)
        oFiles.append(oFile)

  if remove:
    os.remove(sdfFile)
  return oFiles


def downloadSDFPubChem(url, oDir, nMols=0, splitSize=None):
  response = requests.get(url, stream=True)  # Stream to handle large files
  if response.status_code == 200:
    filename = os.path.split(url)[-1]
    gzFile = os.path.join(oDir, filename)
    with open(gzFile, "wb") as f:
      for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
        f.write(chunk)

    sdfFile = os.path.abspath(gunzipFile(gzFile))

    if 'CID-SMILES.smi' in sdfFile:
      if nMols > 0:
        reduceNRandomLines(sdfFile, nMols, oDir=oDir)
      swapColumns(sdfFile)
      splitFile(sdfFile, b=splitSize, pref='smis')

  else:
    raise Exception(f"Failed to download {url}")

def gunzipFile(gzFile, remove=True):
  if 'CID-SMILES.gz' in gzFile:
    ngzFile = gzFile.replace('.gz', '.smi.gz')
    os.rename(gzFile, ngzFile)
    gzFile = ngzFile

  sdfFile = gzFile.replace('.gz', '')
  with gzip.open(gzFile, 'rb') as f:
    sdfText = f.read()
    with open(sdfFile, 'wb') as fo:
      fo.write(sdfText)

  if remove:
    os.remove(gzFile)
  return sdfFile

class ProtChemImportSmallMolecules(EMProtocol):
  """Import small molecules from a directory. Each molecule should be in a separate file.
     Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
     are accepted.
  """
  _label = 'import small mols'
  supportedExt = ['smi', 'mol2', 'sdf', 'pdb', 'mae', 'maegz', 'mol', 'pdbqt']

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')

    form.addParam('defLibraries', BooleanParam, default=False, label='Download a predefined library: ',
                  help='Download a predefined library from a website')
    group = form.addGroup('Default libraries', condition='defLibraries')
    group.addParam('choicesLibraries', EnumParam, default=0,
                   label='Download ligand library: ', condition='defLibraries',
                   choices=libraries,
                   help='Choose the predefined library you want to download.\n'
                        'ECBL: https://www.eu-openscreen.eu/services/compound-collection.html')

    group.addParam('choicesECBL', EnumParam, default=0,
                   label='Download ECBL library: ', condition='defLibraries and choicesLibraries == 0',
                   choices=list(urlECBLDic.keys()),
                   help='Choose the predefined library you want to download.\n'
                        'ECBL: https://www.eu-openscreen.eu/services/compound-collection.html')

    group.addParam('choicesZINC', EnumParam, default=0,
                   label='Download ZINC library: ', condition='defLibraries and choicesLibraries == 1',
                   choices=list(urlZINCJsonDic.keys()),
                   help='Choose the predefined library you want to download.\n'
                        'ZINC: https://zinc.docking.org/substances/subsets/')
    group.addParam('fromSmiles', BooleanParam, default=False, condition='defLibraries and choicesLibraries == 1',
                   label='Format molecule from Smiles: ', expertLevel=LEVEL_ADVANCED,
                   help='Download just the smiles from ZINC and then optimize their structure using '
                        'openbabel ort rdkit instead of downloading the structure.')

    group.addParam('choicesPubChem', EnumParam, default=0,
                   label='Download PubChem library: ', condition='defLibraries and choicesLibraries == 2',
                   choices=list(urlPubChemDic.keys()),
                   help='Choose the molecules library you want to download.\n'
                        'Compound_3D: around 15M molecules with 3D coordinates generated in SDF files.'
                        'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/Release_notes.txt\n'
                        'All: around 120M molecules (2025) with SMIs of all stored molecules in PubChem.'
                        'This is a listing of all CIDs with their isomeric SMILES.')
    group.addParam('nChunks', IntParam, default=10, label='Number of desired molecule chunks: ',
                   condition='defLibraries and choicesLibraries == 2 and choicesPubChem == 0',
                   help='PubChem database contains around 15M different molecules, stored in 25000 molecule IDs '
                        'chunks. This protocol will sample the desired number of these chunks.')
    group.addParam('nMolsPubChem', IntParam, default=10000, label='Number of desired molecules: ',
                   condition='defLibraries and choicesLibraries == 2 and choicesPubChem == 1',
                   help='From the whole set of PubChem molecules, save only this number of random ones.')
    group.addParam('pubChemSeed', IntParam, default=44, expertLevel=LEVEL_ADVANCED,
                   label='Seed for random sampling: ', condition='defLibraries and choicesLibraries == 2',
                   help='The desired number of molecules will be sampled from the PubChem database using this seed')

    group = form.addGroup('Local files', condition='not defLibraries')
    group.addParam('multipleFiles', BooleanParam, default=True, condition='not defLibraries',
                   label='Each file is a molecule: ',
                   help='Whether to import each molecule from a file or all molecules stored in a single file')
    group.addParam('filesPath', PathParam, condition='multipleFiles and not defLibraries',
                   label="Files directory: ",
                   help="Directory with the files you want to import.\n\n"
                        "The path can also contain wildcards to select"
                        "from several folders. \n\n"
                        "Examples:\n"
                        "  ~/project/data/day??_files/\n"
                        "Each '?' represents one unknown character\n\n"
                        "  ~/project/data/day*_files/\n"
                        "'*' represents any number of unknown characters\n\n"
                        "  ~/project/data/day##_files/\n"
                        "'##' represents two digits that will be used as "
                        "file ID\n\n"
                        "NOTE: wildcard characters ('*', '?', '#') "
                        "cannot appear in the actual path.)")
    group.addParam('filesPattern', StringParam, condition='multipleFiles and not defLibraries',
                   label='Pattern', default="*",
                   help="Pattern of the files to be imported.\n\n"
                        "The pattern can contain standard wildcards such as\n"
                        "*, ?, etc, or special ones like ### to mark some\n"
                        "digits in the filename as ID.\n\n"
                        "NOTE: wildcards and special characters "
                        "('*', '?', '#', ':', '%') cannot appear in the "
                        "actual path.\n\n"
                        "You may create small molecules from Smiles (.smi), Tripos Mol2 (.mol2), "
                        "SDF (.sdf), Maestro (.mae, .maegz), or PDB blocks (.pdb)")

    group.addParam('filePath', PathParam, condition='not multipleFiles and not defLibraries',
                   label='File', help='Allowed formats: \n '
                                      ' - SMILES ("SMILES ID" per line) \n'
                                      ' - Mol2 (Multiple mol2 file) \n'
                                      ' - SDF (Multiple sdf file)')

    group = form.addGroup('Output', condition='(defLibraries and choicesLibraries == 2 and choicesPubChem == 1) or (not defLibraries and not multipleFiles)')

    group.addParam('splitSize', IntParam, label='Format files size (MB): ', expertLevel=LEVEL_ADVANCED, default=1,
                   condition='defLibraries and choicesLibraries == 2 and choicesPubChem == 1',
                   help='Maximum size (MB) of the files sent to format to avoid memory issues')
    group.addParam('maxFormatThreads', IntParam, label='Maximum number of format threads: ', expertLevel=LEVEL_ADVANCED,
                   condition='defLibraries and choicesLibraries == 2 and choicesPubChem == 1',
                   help='Maximum number of format threads to use to avoid memory issues', default=10)

    group = form.addGroup('Molecules manager')
    group.addParam('useManager', EnumParam, default=0, label='Manage structure using: ',
                   choices=['RDKit', 'OpenBabel'],
                   help='Whether to manage the structure (parse and optimize if needed) using RDKit or OpenBabel')
    group.addParam('make3d', BooleanParam, default=False, label='Optimize 3D structure: ',
                   help='Optimize 3D structure of the molecules using manager. '
                        'It is automatically done for smi, and very recommendable for 2D structures.\n'
                        'ECBL provide 2D structures.\nZINC provide 3D structures, so this would not be needed')
    group.addParam('keepSMI', BooleanParam, default=False, label='Output SMI molecules: ',
                   expertLevel=LEVEL_ADVANCED, condition='not make3d',
                   help='Usually, SMIs imported are converted to a 3D representation in any case.'
                        'Set this to True if you want the output molecules to be just SMIs '
                        '(will affect visualization and maybe other features)')
    group.addParam('nameKey', StringParam, default='', label='Molecule name key: ',
                   condition='not multipleFiles', expertLevel=LEVEL_ADVANCED,
                   help='Key of the section in the file that defines the name of the molecule. If present,'
                        'it will be used to name molecules')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    if self.defLibraries:
      self._insertFunctionStep('downloadStep')
    if not self.isDirectDownload():
      self._insertFunctionStep('formatStep')
    self._insertFunctionStep('createOutputStep')

  def downloadStep(self):
    oDir = self._getExtraPath() if self.isDirectDownload() else self._getTmpPath()
    libName = self.getEnumText('choicesLibraries')
    if libName == 'ECBL':
      url = urlECBLDic[self.getEnumText('choicesECBL')]
      r = requests.get(url, allow_redirects=True)
      inFile = self.getDownloadFiles(oDir)[0]
      open(inFile, 'wb').write(r.content)

    elif libName == 'ZINC':
      url = urlZINCJsonDic[self.getEnumText('choicesZINC')]
      r = requests.get(url, allow_redirects=True)
      zIds = self.downloadJsonZINC(r, oDir)
      if not self.fromSmiles.get():
        self.downloadSDFThread(zIds, oDir, db=libName)
        os.remove(self.getDownloadFiles(oDir, True)[0])

    elif libName == 'PubChem':
      url = self.getPubChemUrls()
      self.downloadSDFThread(url, oDir, db=libName)

  def formatStep(self):
    outDir = os.path.abspath(self._getExtraPath())
    make3d, nameKey = self.make3d.get(), self.getNameKey()
    nt = min(self.numberOfThreads.get() - 1, self.maxFormatThreads.get())
    keepSMI = self.keepSMI.get()

    if self.isLocalMultiple():
      filesPath, filesPattern = self.filesPath.get(), self.filesPattern.get()
      fnSmalls = [os.path.join(filesPath, filename) for filename in
                  glob.glob(os.path.join(os.getcwd(), filesPath, filesPattern))]
      kwargs = {"make3d": make3d, "nameKey": False, "keepSMI": keepSMI, "outDir": outDir,
                "setOutName": True, "multiInput": False}

    else:
      fnSmalls = self.getDownloadFiles(self._getTmpPath())
      kwargs = {"make3d": make3d, "nameKey": nameKey, "keepSMI": keepSMI, "outDir": outDir,
                "setOutName": False, "multiInput": True}
    # self.formatMolecule(fnSmalls, [], 1, **kwargs)
    performBatchThreading(self.formatMolecule, fnSmalls, nt, cloneItem=False, **kwargs)

  def createOutputStep(self):
    outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
    for fnSmall in glob.glob(self._getExtraPath("*")):
      smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
      smallMolecule.setMolName(os.path.splitext(os.path.basename(fnSmall))[0])
      outputSmallMolecules.append(smallMolecule)

    outputSmallMolecules.updateMolClass()
    self._defineOutputs(outputSmallMolecules=outputSmallMolecules)


  ################# MAIN FUNCTIONS #####################

  def getPubChemUrls(self):
    rd.seed(self.pubChemSeed.get())
    pubChemChoice = self.getEnumText('choicesPubChem')
    pubChemUrl = urlPubChemDic[pubChemChoice]

    if pubChemChoice == 'Compound_3D':
      response = requests.get(pubChemUrl)
      if response.status_code == 200:
        files = extractLinks(response.text, pubChemUrl)
        nChunks = min(self.nChunks.get(), len(files))
        urls = rd.sample(files, nChunks)
      else:
        raise Exception("Failed to access PubChem url")
    else:
      urls = [pubChemUrl]

    return urls

  def getDownloadFiles(self, oDir, json=False):
    inFiles = []
    if self.defLibraries:
      libName = self.getEnumText('choicesLibraries')
      if libName == 'ZINC' and json:
        name, ext = self.getEnumText('choicesZINC'), '.smi'
        inFiles = [os.path.join(oDir, name + ext)]
      elif libName == 'ECBL':
        name, ext = self.getEnumText('choicesECBL'), '.sdf'
        inFiles = [os.path.join(oDir, name + ext)]
      elif libName == 'PubChem' or libName == 'ZINC' and not json:
        for file in os.listdir(oDir):
          inFiles.append(os.path.join(oDir, file))

    else:
      inFiles = [self.filePath.get()]
    return [os.path.abspath(inFile) for inFile in inFiles]

  def downloadSDFThread(self, ids, oDir, db='ZINC'):
    nt = self.numberOfThreads.get()
    if db == 'ZINC':
      oFileBase = os.path.join(oDir, 'ZINC_structures')
      nChunks = len(ids)//200 + 2
      idLists = makeSubsets(ids, nChunks, cloneItem=False)
      idLists = [(i, iList) for i, iList in enumerate(idLists)]
      runInParallel(downloadSDFZINCBatch, oFileBase, paramList=idLists, jobs=nt)
    else:
      runInParallel(downloadSDFPubChem, oDir, self.nMolsPubChem.get(), self.splitSize.get(), paramList=ids, jobs=nt)

  def formatMolecule(self, fnSmalls, oLists, it, make3d, nameKey, keepSMI, outDir, setOutName, multiInput):
    for fnSmall in fnSmalls:
      inFormat = os.path.splitext(fnSmall)[1]
      if make3d or multiInput or inFormat in ['.mae', '.maegz'] or (inFormat == '.smi' and not keepSMI):
        outFormat = os.path.splitext(fnSmall)[1][1:]
        if outFormat in ['smi', 'smiles'] and (not keepSMI or make3d):
          outFormat = DEFAULT_FORMAT

        args = f' -i "{fnSmall}" -of {outFormat} --outputDir {outDir} '
        if setOutName:
          outName = os.path.splitext(os.path.basename(fnSmall))[0]
          args += f'--outputName {outName} '

        if make3d:
          args += '--make3D '
        if nameKey:
          args += f'--nameKey {nameKey} '

        if self.useManager.get() == RDKIT and inFormat != '.mol2':
          # Formatting with RDKit
          Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir, popen=True)

        elif self.useManager.get() == OPENBABEL or inFormat == '.mol2':
          # Formatting with OpenBabel
          Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir, popen=True)

        if make3d:
          self.downloadErrors(outDir)

      # No need of formatting, just copying
      else:
        copyFile(fnSmall, os.path.join(outDir, os.path.basename(fnSmall)))

  #################### UTILS FUNCTIONS #################

  def downloadErrors(self, outDir):
    zIdErrors = []
    for errFile in glob.glob(os.path.join(outDir, 'errors3D_*')):
      errFile = os.path.join(outDir, errFile)
      with open(errFile) as f:
        zIdErrors += f.read().split('\n')[:-1]
      os.remove(errFile)

    if len(zIdErrors) > 0:
      self.downloadSDFThread(zIdErrors)

  def downloadJsonZINC(self, r, oDir):
    content, zIds = '', []
    for zjson in eval(r.content):
      jDic = json.loads(str(zjson).replace("'", '"'))
      content += '{} {}\n'.format(jDic['smiles'], jDic['zinc_id'])
      zIds.append(jDic['zinc_id'])

    inFile = self.getDownloadFiles(oDir, True)[0]
    open(inFile, 'w').write(content)
    return zIds

  def isDirectDownload(self):
    return self.defLibraries and self.getEnumText('choicesLibraries') == 'ZINC' and not self.fromSmiles

  def isLocalMultiple(self):
    '''Whether it is local multiple files (each a single molecule)
    '''
    return not self.defLibraries and self.multipleFiles

  def getNameKey(self):
    nameKey = None
    if not self.multipleFiles.get() and self.nameKey.get():
      nameKey = self.nameKey.get().strip()
    elif self.defLibraries:
      if self.getEnumText('choicesLibraries') == 'ZINC':
        nameKey = 'zinc_id'
    return nameKey


############# VALIDATION FUNCTIONS ##################

  def _validate(self):
    errors = []
    if not self.defLibraries:
      if not self.multipleFiles.get():
        ext = os.path.splitext(self.filePath.get())[1][1:]
        if not ext in self.supportedExt:
          errors.append('Unknown input file format {}\n'
                        'Recognized formats: {}'.format(ext, self.supportedExt))
      else:
        for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
          ext = os.path.splitext(filename)[1][1:]
          if not ext in self.supportedExt:
            errors.append('Unknown input file format {}\n'
                          'Recognized formats: {}'.format(ext, self.supportedExt))
            break
    return errors