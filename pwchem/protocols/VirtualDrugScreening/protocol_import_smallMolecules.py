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

import os, requests, glob, sys, json


from pwem.protocols import EMProtocol
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam, LEVEL_ADVANCED, EnumParam, STEPS_PARALLEL

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin
from pwchem.utils import performBatchThreading

RDKIT, OPENBABEL = 0, 1
DEFAULT_FORMAT = 'sdf'

libraries = ['ECBL', 'ZINC']

urlECBLDic = {'Bioactive (~5000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/Pilot_08_09_2021.sdf',
              'Diversity (~10⁵)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/ECBL_08_09_2021.sdf',
              'Fragments (~1000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Fragments.sdf'}

urlZINCJsonDic = {'FDA (~1300)': 'https://zinc.docking.org/substances/subsets/fda.json?count=all',
                  'Drugs-NotFDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.json?count=all'}

urlZINCDic = {'FDA (~1300)': 'https://zinc.docking.org/substances/subsets/fda.sdf?count=all',
              'Drugs-NotFDA (~2000)': 'https://zinc.docking.org/substances/subsets/world-not-fda.sdf?count=all'}

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
        group.addParam('filesPattern', StringParam,  condition='multipleFiles and not defLibraries',
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

        print('Importing molecules from ', url)
        sys.stdout.flush()
        r = requests.get(url, allow_redirects=True)

        if libName == 'ZINC':
            zIds = self.downloadJsonZINC(r)
            if not self.fromSmiles.get():
                self.downloadSDFZINCThread(zIds)

        else:
            inFile = self.getDownloadFile()
            open(inFile, 'wb').write(r.content)

    def formatStep(self):
        make3d, mulFiles = self.make3d.get(), self.multipleFiles.get()
        filesPath = self.filesPath.get()
        nameKey = self.getNameKey()

        if self.isLocalMultiple():
            # Format multiple local files placed in filesPath with filesPattern
            outDir = os.path.abspath(self._getExtraPath())
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnSmall = os.path.join(self.filesPath.get(), filename)

                # Files need to be converted if they are maestro or they are asked to be optimized 3D
                if fnSmall.endswith(".mae") or fnSmall.endswith(".maegz") or make3d:
                  outName = os.path.splitext(os.path.basename(fnSmall))[0]
                  outFormat = os.path.splitext(fnSmall)[1][1:]
                  if outFormat in ['smi', 'smiles']:
                    outFormat = DEFAULT_FORMAT

                  args = ' -i "{}" --outputName {} -of {} --outputDir "{}"'. format(fnSmall, outName, outFormat, outDir)
                  if make3d:
                    args += ' --make3D -nt {}'.format(self.numberOfThreads.get())

                  # Formatting with RDKit (neccessary if they are maestro)
                  if self.useManager.get() == RDKIT or fnSmall.endswith(".mae") or fnSmall.endswith(".maegz"):
                      Plugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=outDir)

                  # Formatting with OpenBabel
                  elif self.useManager.get() == OPENBABEL:
                      Plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

                  if make3d:
                      self.downloadErrors(outDir)

                # No need of formatting, just copying
                else:
                    copyFile(fnSmall, os.path.join(outDir, os.path.basename(fnSmall)))

        else:
            # Format single local file or downloaded
            fnSmall = self.getDownloadFile()
            outDir = os.path.abspath(self._getExtraPath())

            outFormat = os.path.splitext(fnSmall)[1][1:]
            if outFormat in ['smi', 'smiles']:
              outFormat = DEFAULT_FORMAT

            args = ' -i "{}" -of {} --outputDir {}'.format(fnSmall, outFormat, outDir)
            if make3d:
                args += ' --make3D -nt {}'.format(self.numberOfThreads.get())
            if nameKey:
                args += ' --nameKey {}'.format(nameKey)

            if self.useManager.get() == RDKIT:
            # Formatting with RDKit
                Plugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=outDir)

            elif self.useManager.get() == OPENBABEL:
            # Formatting with OpenBabel
                Plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

            if make3d:
              self.downloadErrors(outDir)

    def createOutputStep(self):
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')
        for fnSmall in glob.glob(self._getExtraPath("*")):
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
            smallMolecule.setMolName(os.path.splitext(os.path.basename(fnSmall))[0])

            outputSmallMolecules.append(smallMolecule)
        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)


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

    def downloadErrors(self, outDir):
        zIdErrors = []
        for errFile in glob.glob(os.path.join(outDir, 'errors3D_*')):
          errFile = os.path.join(outDir, errFile)
          with open(errFile) as f:
            zIdErrors += f.read().split('\n')[:-1]
          os.remove(errFile)

        if len(zIdErrors) > 0:
            self.downloadSDFZINCThread(zIdErrors)

    def downloadSDFZINCThread(self, zIDs):
        nt = self.numberOfThreads.get()
        performBatchThreading(self.downloadSDFZINC, zIDs, nt, cloneItem=False)

    def downloadSDFZINC(self, zIds, outLists, it):
        for zId in zIds:
            urlSdf = 'https://zinc.docking.org/substances/{}.sdf'.format(zId)
            print('Importing molecules from ', urlSdf)
            rSdf = requests.get(urlSdf, allow_redirects=True)
            oFile = self._getExtraPath('{}.sdf'.format(zId))
            open(oFile, 'wb').write(rSdf.content)
        return zIds

    def downloadJsonZINC(self, r):
        content, zIds = '', []
        for zjson in eval(r.content):
          jDic = json.loads(str(zjson).replace("'", '"'))
          content += '{} {}\n'.format(jDic['smiles'], jDic['zinc_id'])
          zIds.append(jDic['zinc_id'])

        inFile = self.getDownloadFile()
        open(inFile, 'w').write(content)
        return zIds

    def getDownloadFile(self):
      if self.defLibraries:
          libName = self.getEnumText('choicesLibraries')
          if libName == 'ZINC':
              name, ext = self.getEnumText('choicesZINC'), '.smi'
          else:
              name, ext = self.getEnumText('choicesECBL'), '.sdf'
          inFile = self._getTmpPath(name + ext)
      else:
          inFile = self.filePath.get()
      return os.path.abspath(inFile)

    def isDirectDownload(self):
        return self.defLibraries and self.getEnumText('choicesLibraries') == 'ZINC' and not self.fromSmiles

    def isLocalMultiple(self):
        return not self.defLibraries and self.multipleFiles

    def getNameKey(self):
      nameKey = None
      if not self.multipleFiles.get() and self.nameKey.get():
          nameKey = self.nameKey.get().strip()
      elif self.defLibraries:
          if self.getEnumText('choicesLibraries') == 'ZINC':
              nameKey = 'zinc_id'
          elif self.getEnumText('choicesLibraries') == 'ECBL':
              nameKey = 'Supplier_ID'
      return nameKey
