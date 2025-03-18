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

import os, requests, glob, sys, json, gzip, subprocess
from bs4 import BeautifulSoup
import random as rd

from pwem.protocols import EMProtocol
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam, LEVEL_ADVANCED, EnumParam, \
  STEPS_PARALLEL, IntParam

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin
from pwchem.utils import performBatchThreading
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC

RDKIT, OPENBABEL = 0, 1
DEFAULT_FORMAT = 'sdf'

libraries = ['ECBL', 'ZINC', 'PubChem']

urlECBLDic = {
  'Bioactive (~5000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/Pilot_08_09_2021.sdf',
  'Diversity (~10⁵)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/ECBL_08_09_2021.sdf',
  'Fragments (~1000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Fragments.sdf'}

urlZINCJsonDic = {'FDA (~1300)': 'https://zinc.docking.org/substances/subsets/fda.json?count=all',
                  'Drugs-NotFDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.json?count=all'}

urlZINCDic = {'FDA (~1300)': 'https://zinc.docking.org/substances/subsets/fda.sdf?count=all',
              'Drugs-NotFDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.sdf?count=all'}

urlPubChemDic = {'Compound_3D': 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound_3D/01_conf_per_cmpd/SDF/',
                 'All': 'https://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Extras/CID-SMILES.gz'}


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
                   label='Download ECBL library: ', condition='defLibraries and choicesLibraries == 1',
                   choices=list(urlZINCDic.keys()),
                   help='Choose the predefined library you want to download.\n'
                        'ZINC: https://zinc.docking.org/substances/subsets/')
    group.addParam('fromSmiles', BooleanParam, default=True, condition='defLibraries and choicesLibraries == 1',
                   label='Format molecule from Smiles: ',
                   help='Downloads from ZINC database are typically not feasible, so it might be easier to download'
                        ' just the smiles and the optimize their structure using openbabel. If not, the protocol '
                        'will download the sdf structure of the molecules in the json one by one, which will be '
                        'slower.')

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
                   help='PubChem database contains around 15M different molecules, stored in 25000 molecule IDs'
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
                                      ' - CSV smiles (ID, compound; this is downloaded from ZINC) \n'
                                      ' - Mol2 (Multiple mol2 file) \n'
                                      ' - SDF (Multiple sdf file)')

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
    libName = self.getEnumText('choicesLibraries')
    if libName == 'ECBL':
      url = urlECBLDic[self.getEnumText('choicesECBL')]
    elif libName == 'ZINC':
      url = urlZINCJsonDic[self.getEnumText('choicesZINC')]
    elif libName == 'PubChem':
      url = self.getPubChemUrls()

    print('Importing molecules from ', url)
    sys.stdout.flush()

    if libName == 'ZINC':
      r = requests.get(url, allow_redirects=True)
      zIds = self.downloadJsonZINC(r)
      if not self.fromSmiles.get():
        self.downloadSDFThread(zIds, db=libName)

    elif libName == 'ECBL':
      r = requests.get(url, allow_redirects=True)
      inFile = self.getDownloadFiles()[0]
      open(inFile, 'wb').write(r.content)

    elif libName == 'PubChem':
      self.downloadSDFThread(url, db=libName)

  def formatStep(self):
    outDir = os.path.abspath(self._getExtraPath())
    make3d, nameKey = self.make3d.get(), self.getNameKey()
    nt = self.numberOfThreads.get()
    keepSMI = self.keepSMI.get()

    if self.isLocalMultiple():
      filesPath, filesPattern = self.filesPath.get(), self.filesPattern.get()
      fnSmalls = [os.path.join(filesPath, filename) for filename in
                  glob.glob(os.path.join(os.getcwd(), filesPath, filesPattern))]
      kwargs = {"make3d": make3d, "nameKey": False, "keepSMI": keepSMI, "outDir": outDir,
                "setOutName": True, "multiInput": False}

    else:
      fnSmalls = self.getDownloadFiles()
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
        soup = BeautifulSoup(response.text, "html.parser")
        files = [a['href'] for a in soup.find_all('a', href=True) if not a['href'].startswith('?')]
        files = [file for file in files if file.endswith('.sdf.gz')]
        nChunks = min(self.nChunks.get(), len(files))
        files = rd.sample(files, nChunks)
        urls = [os.path.join(pubChemUrl, file) for file in files]
      else:
        raise Exception("Failed to access PubChem url")
    else:
      urls = [pubChemUrl]

    return urls

  def getDownloadFiles(self):
    inFiles = []
    if self.defLibraries:
      libName = self.getEnumText('choicesLibraries')
      if libName == 'ZINC':
        name, ext = self.getEnumText('choicesZINC'), '.smi'
        inFiles = [self._getTmpPath(name + ext)]
      elif libName == 'ECBL':
        name, ext = self.getEnumText('choicesECBL'), '.sdf'
        inFiles = [self._getTmpPath(name + ext)]
      elif libName == 'PubChem':
        for file in os.listdir(self._getTmpPath()):
          inFiles.append(self._getTmpPath(file))

    else:
      inFiles = [self.filePath.get()]
    return [os.path.abspath(inFile) for inFile in inFiles]

  def downloadSDFThread(self, ids, db='ZINC'):
    nt = self.numberOfThreads.get()
    downFunc = self.downloadSDFZINC if db == 'ZINC' else self.downloadSDFPubChem
    performBatchThreading(downFunc, ids, nt, cloneItem=False)


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

  def downloadSDFZINC(self, zIds, outLists, it):
    for zId in zIds:
      urlSdf = 'https://zinc.docking.org/substances/{}.sdf'.format(zId)
      print('Importing molecules from ', urlSdf)
      rSdf = requests.get(urlSdf, allow_redirects=True)
      oFile = self._getExtraPath('{}.sdf'.format(zId))
      open(oFile, 'wb').write(rSdf.content)
    return zIds

  def downloadSDFPubChem(self, urls, oList, it):
    for url in urls:
      response = requests.get(url, stream=True)  # Stream to handle large files
      if response.status_code == 200:
        filename = os.path.split(url)[-1]
        gzFile = self._getTmpPath(filename)
        with open(gzFile, "wb") as f:
          for chunk in response.iter_content(chunk_size=1024):  # Download in chunks
            f.write(chunk)

        sdfFile = self.gunzipFile(gzFile)

        if 'CID-SMILES.smi' in sdfFile:
          nMols = self.nMolsPubChem.get()
          if nMols > 0:
            self.reduceNRandomLines(sdfFile, nMols)
          self.swapColumns(sdfFile)

      else:
        raise Exception(f"Failed to download {url}")

  def reduceNRandomLines(self, inFile, n):
    tFile = os.path.abspath(self._getTmpPath('temp.smi'))
    subprocess.check_call(f'shuf -n {n} {os.path.abspath(inFile)} > {tFile}', shell=True)
    os.rename(tFile, inFile)

  def swapColumns(self, smiFile):
    smiFile = os.path.abspath(smiFile)
    tFile = os.path.abspath(self._getTmpPath('smis.smi'))
    subprocess.check_call(f"awk '{{ t=$1; $1=$2; $2=t; print }}' {smiFile} > {tFile}", shell=True)
    os.rename(tFile, smiFile)

  def gunzipFile(self, gzFile, remove=True):
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

  def downloadJsonZINC(self, r):
    content, zIds = '', []
    for zjson in eval(r.content):
      jDic = json.loads(str(zjson).replace("'", '"'))
      content += '{} {}\n'.format(jDic['smiles'], jDic['zinc_id'])
      zIds.append(jDic['zinc_id'])

    inFile = self.getDownloadFiles()[0]
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