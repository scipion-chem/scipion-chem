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

import os, requests, glob, re, json, gzip, sys
from html.parser import HTMLParser
from urllib.parse import urljoin
import random as rd
from itertools import product

from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin
from pwchem.utils import performBatchThreading, runInParallel, makeSubsets, downloadUrlFile, getBaseName
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC
from pwchem.protocols.VirtualDrugScreening.protocol_import_molecules_library import ProtChemImportMoleculesLibrary

RDKIT, OPENBABEL = 0, 1
DEFAULT_FORMAT = 'sdf'

libraries = ['ECBL', 'ZINC', 'PubChem']
ECBL, ZINC, PUBCHEM = 0, 1, 2

urlECBLDic = {
  'Bioactive (~5000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/Pilot_08_09_2021.sdf',
  'Diversity (~10⁵)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/ECBL_08_09_2021.sdf',
  'Fragments (~1000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Fragments.sdf'}

urlZINCJsonDic = {'FDA (~1600)': 'https://zinc.docking.org/substances/subsets/fda.json?count=all',
                  'World Not FDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.json?count=all',
                  'World (~3600)': 'https://zinc.docking.org/substances/subsets/world.json?count=all',
                  'Investigational only (~3900)': 'https://zinc.docking.org/substances/subsets/investigational-only.json?count=all',
                  'Non-human metabolites (~2000)': 'https://zinc.docking.org/substances/subsets/nonhuman-metabolites.json?count=all',
                  }

pubChemUrl = 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/'

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

def extractLinks(htmlContent, baseUrl, pattern=r'\.sdf\.gz$'):
  parser = LinkParser(baseUrl)
  parser.feed(htmlContent)

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
        molText = f'{zincId}\n\t{molText}\n\n$$$$'

        oFile = os.path.join(oDir, f'{zincId}.sdf')
        with open(oFile, 'w') as f:
          f.write(molText)
        oFiles.append(oFile)

  if remove:
    os.remove(sdfFile)
  return oFiles


def downloadSDFPubChem(url, oDir):
  response = requests.get(url, stream=True)  # Stream to handle large files
  if response.status_code == 200:
    filename = os.path.split(url)[-1]
    gzFile = os.path.join(oDir, filename)
    with open(gzFile, "wb") as f:
      for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
        f.write(chunk)

    os.path.abspath(gunzipFile(gzFile))

  else:
    raise Exception(f"Failed to download {url}")

def gunzipFile(gzFile, remove=True):
  sdfFile = gzFile.replace('.gz', '')
  with gzip.open(gzFile, 'rb') as f:
    sdfText = f.read()
    with open(sdfFile, 'wb') as fo:
      fo.write(sdfText)

  if remove:
    os.remove(gzFile)
  return sdfFile

class ProtChemImportSmallMolecules(ProtChemImportMoleculesLibrary):
  """Import small molecules from a directory. Each molecule should be in a separate file.
     Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
     are accepted.
  """
  _label = 'import small mols'
  supportedExt = ['smi', 'mol2', 'sdf', 'pdb', 'mae', 'maegz', 'mol', 'pdbqt']
  stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')

    form.addParam('defLibraries', params.BooleanParam, default=False, label='Download a predefined library: ',
                  help='Download a predefined library from a website')
    group = form.addGroup('Default libraries', condition='defLibraries')
    group.addParam('choicesLibraries', params.EnumParam, default=ECBL,
                   label='Download ligand library: ', condition='defLibraries', choices=libraries,
                   help='Choose the predefined library you want to download.\n'
                        'ECBL: https://www.eu-openscreen.eu/services/compound-collection.html')

    group.addParam('choicesECBL', params.EnumParam, default=0, choices=list(urlECBLDic.keys()),
                   label='Download ECBL library: ', condition=f'defLibraries and choicesLibraries == {ECBL}',
                   help='Choose the predefined library you want to download.\n'
                        'ECBL: https://www.eu-openscreen.eu/services/compound-collection.html')

    group.addParam('fromTranches', params.BooleanParam, default=False,
                   label='Download ZINC tranches: ', condition=f'defLibraries and choicesLibraries == {ZINC}',
                   help='Download base on ZINC tranches (Based on molecular weight and logP) instead of predefined '
                        'small subsets.')

    group.addParam('choicesZINC', params.EnumParam, default=0, choices=list(urlZINCJsonDic.keys()),
                   label='Download ZINC library: ',
                   condition=f'defLibraries and choicesLibraries == {ZINC} and not fromTranches',
                   help='Choose the predefined library you want to download.\n'
                        'ZINC: https://zinc.docking.org/substances/subsets/')

    group.addParam('fromSmiles', params.BooleanParam, default=False,
                   condition=f'defLibraries and choicesLibraries == {ZINC}',
                   label='Format molecule from Smiles: ', expertLevel=params.LEVEL_ADVANCED,
                   help='Download just the smiles from ZINC and then optimize their structure using '
                        'openbabel ort rdkit instead of downloading the structure.')

    groupT = form.addGroup('ZINC20 tranches', condition=f'defLibraries and choicesLibraries == {ZINC} and fromTranches')
    groupT = self._defineZINC20TranchesParams(groupT,
                                             condition=f'defLibraries and choicesLibraries == {ZINC} and fromTranches')

    groupT.addParam('repPH', params.StringParam, label='Representation PH(s): ', expertLevel=params.LEVEL_ADVANCED, default='R, M',
                   condition=f'defLibraries and choicesLibraries == {ZINC} and fromTranches and not fromSmiles',
                   help='Representation at determined PH, comma separated. Possible values: R, M, H, L\n'
                        'Reference (R), Mid (M), High (H) and Low (L)')
    groupT.addParam('repCharge', params.StringParam, label='Allowed charges: ', expertLevel=params.LEVEL_ADVANCED, default='L, M, N, O, P',
                   condition=f'defLibraries and choicesLibraries == {ZINC} and fromTranches and not fromSmiles',
                   help='Allowed charges, comma separated. Possible values: L, M, N, O, P\n'
                        '-2 (L), -1 (M), 0 (O), P (+1) and Q (+2)')

    group.addParam('nChunks', params.IntParam, default=10, label='Number of desired molecule chunks: ',
                   condition=f'defLibraries and choicesLibraries == {PUBCHEM}',
                   help='PubChem database contains around 15M different molecules, stored in 25000 molecule IDs '
                        'chunks. This protocol will sample the desired number of these chunks.')
    group.addParam('pubChemSeed', params.IntParam, default=44, expertLevel=params.LEVEL_ADVANCED,
                   label='Seed for random sampling: ', condition=f'defLibraries and choicesLibraries == {PUBCHEM}',
                   help='The desired number of molecules will be sampled from the PubChem database using this seed')

    group = form.addGroup('Local files', condition='not defLibraries')
    group.addParam('singleFiles', params.BooleanParam, default=True, condition='not defLibraries',
                   label='Each file is a molecule: ',
                   help='Whether to import each molecule from a file or all molecules stored in a single file')
    group.addParam('filesPath', params.PathParam, condition='singleFiles and not defLibraries',
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
    group.addParam('filesPattern', params.StringParam, condition='singleFiles and not defLibraries',
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

    group.addParam('filePath', params.PathParam, condition='not singleFiles and not defLibraries',
                   label='File', help='Allowed formats: \n '
                                      ' - SMILES ("SMILES ID" per line) \n'
                                      ' - Mol2 (Multiple mol2 file) \n'
                                      ' - SDF (Multiple sdf file)')

    group = form.addGroup('Molecules manager')
    group.addParam('useManager', params.EnumParam, default=0, label='Manage structure using: ',
                   choices=['RDKit', 'OpenBabel'],
                   help='Whether to manage the structure (parse and optimize if needed) using RDKit or OpenBabel')
    group.addParam('make3d', params.BooleanParam, default=False, label='Optimize 3D structure: ',
                   help='Optimize 3D structure of the molecules using manager. '
                        'It is automatically done for smi, and very recommendable for 2D structures.\n'
                        'ECBL provide 2D structures.\nZINC provide 3D structures, so this would not be needed')
    group.addParam('nameKey', params.StringParam, default='', label='Molecule name key: ',
                   condition='not singleFiles', expertLevel=params.LEVEL_ADVANCED,
                   help='Key of the section in the file that defines the name of the molecule. If present,'
                        'it will be used to name molecules')

    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- INSERT steps functions --------------------
  def _insertAllSteps(self):
    if self.defLibraries:
      self._insertFunctionStep('downloadStep')
    self._insertFunctionStep('formatStep')
    self._insertFunctionStep('createOutputStep')

  def downloadStep(self):
    oDir = self._getTmpPath()
    libName = self.getEnumText('choicesLibraries')
    if libName == 'ECBL':
      url = urlECBLDic[self.getEnumText('choicesECBL')]
      r = requests.get(url, allow_redirects=True)
      inFile = self.getDownloadFiles(oDir)[0]
      open(inFile, 'wb').write(r.content)

    elif libName == 'ZINC':
      if not self.fromTranches.get():
        url = urlZINCJsonDic[self.getEnumText('choicesZINC')]
        r = requests.get(url, allow_redirects=True)
        zIds = self.downloadJsonZINC(r, oDir)
        if not self.fromSmiles.get():
          self.downloadSDFThread(zIds, oDir, db=libName)
          os.remove(self.getDownloadFiles(oDir, True)[0])

      else:
        self.getZINC20Range(oDir)

    elif libName == 'PubChem':
      url = self.getPubChemUrls()
      self.downloadSDFThread(url, oDir, db=libName)

  def formatStep(self):
    outDir = os.path.abspath(self._getExtraPath())
    make3d, nameKey = self.make3d.get(), self.getNameKey()
    nt = self.numberOfThreads.get() - 1

    if self.isLocalSingleFiles():
      filesPath, filesPattern = self.filesPath.get(), self.filesPattern.get()
      fnSmalls = [os.path.join(filesPath, filename) for filename in
                  glob.glob(os.path.join(os.getcwd(), filesPath, filesPattern))]
      kwargs = {"make3d": make3d, "nameKey": False, "outDir": outDir,
                "setOutName": True, "multiInput": False, 'dirPerFile': False}

    else:
      fnSmalls = self.getDownloadFiles(self._getTmpPath())
      kwargs = {"make3d": make3d, "nameKey": nameKey, "outDir": outDir,
                "setOutName": False, "multiInput": True, 'dirPerFile': False}
      if self.defLibraries.get() and self.choicesLibraries.get() == ZINC:
        if self.fromTranches.get():
          kwargs['dirPerFile'] = True
        else:
          kwargs['multiInput'] = False
    performBatchThreading(self.formatMolecule, fnSmalls, nt, cloneItem=False, **kwargs)

  def createOutputStep(self):
    zincDic = {}
    if self.defLibraries.get() and self.choicesLibraries.get() == ZINC and self.fromTranches.get():
      zincDic = self.getOutZincDic()

    outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
    for fnSmall in glob.glob(self._getExtraPath("*")):
      smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
      smallMolecule.setMolName(os.path.splitext(os.path.basename(fnSmall))[0])
      if fnSmall in zincDic:
        smallMolecule._size = pwobj.String(zincDic[fnSmall][0])
        smallMolecule._logP = pwobj.String(zincDic[fnSmall][1])
        smallMolecule._reactivity = pwobj.String(zincDic[fnSmall][2])
        smallMolecule._purchasability = pwobj.String(zincDic[fnSmall][3])

      outputSmallMolecules.append(smallMolecule)

    outputSmallMolecules.updateMolClass()
    self._defineOutputs(outputSmallMolecules=outputSmallMolecules)


  ################# MAIN FUNCTIONS #####################

  def getPubChemUrls(self):
    rd.seed(self.pubChemSeed.get())
    response = requests.get(pubChemUrl)

    files = extractLinks(response.text, pubChemUrl)
    nChunks = min(self.nChunks.get(), len(files))
    urls = rd.sample(files, nChunks)
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
      runInParallel(downloadSDFPubChem, oDir, paramList=ids, jobs=nt)

  def formatMolecule(self, fnSmalls, oLists, it, make3d, nameKey, outDir, setOutName, multiInput, dirPerFile):
    for fnSmall in fnSmalls:
      wDir = os.path.join(outDir, getBaseName(fnSmall)) if dirPerFile else outDir
      if not os.path.exists(wDir):
        os.mkdir(wDir)

      inFormat = os.path.splitext(fnSmall)[1]
      if make3d or multiInput or inFormat in ['.mae', '.maegz'] or inFormat == '.smi':
        outFormat = os.path.splitext(fnSmall)[1][1:]
        if outFormat in ['smi', 'smiles'] or make3d:
          outFormat = DEFAULT_FORMAT

        args = f' -i "{fnSmall}" -of {outFormat} --outputDir {wDir} '
        if setOutName:
          outName = os.path.splitext(os.path.basename(fnSmall))[0]
          args += f'--outputName {outName} '

        if make3d:
          args += '--make3D '
        if nameKey:
          args += f'--nameKey {nameKey} '

        if self.useManager.get() == RDKIT and inFormat != '.mol2':
          # Formatting with RDKit
          Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=wDir, popen=True)

        elif self.useManager.get() == OPENBABEL or inFormat == '.mol2':
          # Formatting with OpenBabel
          Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=wDir, popen=True)

        if make3d:
          self.downloadErrors(wDir)

      # No need of formatting, just copying
      else:
        copyFile(fnSmall, os.path.join(wDir, os.path.basename(fnSmall)))

  def getZINC20Range(self, oDir):
    sizeCodes, logPCodes, reactCodes, purchCodes = self.getZINC20Codes()

    phCodes = [c.strip() for c in self.repPH.get().split(',')]
    chargeCodes = [c.strip() for c in self.repCharge.get().split(',')]

    downTranches = [''.join(p) for p in product(sizeCodes, logPCodes, reactCodes, purchCodes, phCodes, chargeCodes)]
    allUrls = [f'https://files.docking.org/3D/{tranche[:2]}/{tranche[2:]}/{tranche}.xaa.sdf.gz' for tranche in downTranches]

    nt = self.numberOfThreads.get()
    gzFiles = runInParallel(downloadUrlFile, oDir, 2, paramList=allUrls, jobs=nt)
    gzFiles = [gzFile for gzFile in gzFiles if os.path.exists(gzFile)]
    sdfFiles = [gunzipFile(gzFile, remove=True) for gzFile in gzFiles]
    return sdfFiles

  #################### UTILS FUNCTIONS #################

  def downloadErrors(self, outDir):
    zIdErrors = []
    for errFile in glob.glob(os.path.join(outDir, 'errors3D_*')):
      errFile = os.path.join(outDir, errFile)
      with open(errFile) as f:
        zIdErrors += f.read().split('\n')[:-1]
      os.remove(errFile)

    if len(zIdErrors) > 0:
      self.downloadSDFThread(zIdErrors, outDir)

  def downloadJsonZINC(self, r, oDir):
    content, zIds = '', []
    for zjson in eval(r.content):
      jDic = json.loads(str(zjson).replace("'", '"'))
      content += f'{jDic["smiles"]} {jDic["zinc_id"]}\n'
      zIds.append(jDic['zinc_id'])

    inFile = self.getDownloadFiles(oDir, True)[0]
    open(inFile, 'w').write(content)
    return zIds

  def isLocalSingleFiles(self):
    '''Whether it is local multiple files (each a single molecule)
    '''
    return not self.defLibraries and self.singleFiles

  def getNameKey(self):
    nameKey = None
    if not self.singleFiles.get() and self.nameKey.get():
      nameKey = self.nameKey.get().strip()
    return nameKey

  def getOutZincDic(self):
    '''Reformats the directory structure in the extra dir to match the expected output, parsing first the tranches
    information for each output molecule as: {oFile: [size, logP, Reactivity, Purchasability]}
    '''
    zincDic = {}
    decodeDics = self.getDecodeZinc20Dics()
    for codeDir in os.listdir(self._getExtraPath()):
      code = codeDir.split('.')[0]
      for oFile in os.listdir(self._getExtraPath(codeDir)):
        zFile = self._getExtraPath(os.path.join(codeDir, oFile))
        oFile = self._getExtraPath(oFile)
        os.rename(zFile, oFile)
        zincDic[oFile] = self.decodeTranche(code, decodeDics)

      os.removedirs(self._getExtraPath(codeDir))
    return zincDic


############# VALIDATION FUNCTIONS ##################

  def _validate(self):
    errors = []
    if not self.defLibraries:
      if not self.singleFiles.get():
        ext = os.path.splitext(self.filePath.get())[1][1:]
        if not ext in self.supportedExt:
          errors.append(f'Unknown input file format {ext}\n'
                        f'Recognized formats: {self.supportedExt}')
      else:
        for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
          ext = os.path.splitext(filename)[1][1:]
          if not ext in self.supportedExt:
            errors.append(f'Unknown input file format {ext}\n'
                          f'Recognized formats: {self.supportedExt}')
            break

    else:
      if self.choicesLibraries.get() == ZINC and self.fromTranches.get():
        phCodes = set([c.strip() for c in self.repPH.get().split(',')])
        if phCodes.difference(set('RMHL')):
          errors.append('PH codes must be specified in comma separated and the only allowed values are: R,M,H,L')

        chargeCodes = set([c.strip() for c in self.repCharge.get().split(',')])
        if chargeCodes.difference(set('LMNOP')):
          errors.append('Charge codes must be specified in comma separated and the only allowed values are: L, M, N, O, P')

    return errors