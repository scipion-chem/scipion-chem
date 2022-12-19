# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Alberto M. Parra Pérez (amparraperez@gmail.com)
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

import os, requests, glob


from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam, LEVEL_ADVANCED, EnumParam

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin
from pwchem.utils import runOpenBabel

urlDic = {'ECBL bioactive (~5000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/Pilot_08_09_2021.sdf',
          'ECBL diversity (~10⁵)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Video-Content/other_uploads/ECBL_08_09_2021.sdf',
          'ECBL fragments (~1000)': 'https://www.eu-openscreen.eu/fileadmin/user_upload/Fragments.sdf'}

class ProtChemImportSmallMolecules(EMProtocol):
    """Import small molecules from a directory. Each molecule should be in a separate file.
       Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
       are accepted.
    """
    _label = 'import small mols'
    supportedExt = ['smi', 'mol2', 'sdf', 'pdb', 'mae', 'maegz', 'mol', 'pdbqt']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('defLibraries', BooleanParam, default=True, label='Download a predefined library: ',
                      help='Download a predefined library from a website')
        form.addParam('choicesLibraries', EnumParam, default=0,
                      label='Download ligand library: ', condition='defLibraries',
                      choices=list(urlDic.keys()),
                      help='Choose the predefined library you want to download.\n'
                           'ECBL: https://www.eu-openscreen.eu/services/compound-collection.html')

        form.addParam('multipleFiles', BooleanParam, default=True, condition='not defLibraries',
                      label='Each file is a molecule: ',
                      help='Whether to import each molecule from a file or all molecules stored in a single file')
        form.addParam('filesPath', PathParam, condition='multipleFiles and not defLibraries',
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
        form.addParam('filesPattern', StringParam,  condition='multipleFiles and not defLibraries',
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

        form.addParam('filePath', PathParam, condition='not multipleFiles and not defLibraries',
                      label='File', help='Allowed formats: \n '
                                         ' - CSV smiles (ID, compound; this is downloaded from ZINC) \n'
                                         ' - Mol2 (Multiple mol2 file) \n'
                                         ' - SDF (Multiple sdf file)')

        # todo: make 3D opt actually everytime is true
        form.addParam('make3d', BooleanParam, default=False, label='Optimize 3D structure',
                      help='Optimize 3D structure of the molecules using the pybel fucntion. '
                           'It is automatically done for smi, and very recommendable for 2D structures')
        form.addParam('draw', BooleanParam, default=False, label='Draw molecule to jpg', expertLevel=LEVEL_ADVANCED,
                      help='Draw molecule using OpenBabel and store it in the output object')
        form.addParam('nameKey', StringParam, default='', label='Molecule name key: ',
                      condition='not multipleFiles', expertLevel=LEVEL_ADVANCED,
                      help='Key of the section in the file that defines the name of the molecule. If present,'
                           'it will be used to name molecules')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        if self.defLibraries:
            libName = self.getEnumText('choicesLibraries')
            url = urlDic[libName]
            r = requests.get(url, allow_redirects=True)
            inFile = self._getTmpPath(os.path.split(url)[-1])
            open(inFile, 'wb').write(r.content)
            mulFiles = False
            nameKey = 'Supplier_ID'

        else:
          mulFiles = self.multipleFiles.get()
          inFile = self.filePath.get()
          nameKey = self.nameKey.get().strip()



        make3d = self.make3d.get()
        if not mulFiles:
            fnSmall = os.path.abspath(inFile)
            if fnSmall.endswith(".pdb"):  # Multiple pdb
                self.parseMultiPDBFile(fnSmall)

            elif fnSmall.endswith(".mae") or fnSmall.endswith(".maegz"):  # Multiple mae
                args = ' -i "{}" --outputDir {}'.format(fnSmall, os.path.abspath(self._getExtraPath()))
                Plugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=self._getExtraPath())

            else:
                outFormat = os.path.splitext(fnSmall)[1][1:]
                if outFormat in ['smi', 'smiles']:
                    outFormat = 'mol2'

                args = ' -i "{}" -of {} --outputDir "{}"'.format(fnSmall, outFormat,
                                                             os.path.abspath(self._getExtraPath()))
                if make3d:
                    args += ' --make3D'
                if nameKey != '':
                    args += ' --nameKey {}'.format(nameKey)
                Plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=self._getExtraPath())

        else:
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnSmall = os.path.join(self.filesPath.get(), filename)
                if fnSmall.endswith(".mae") or fnSmall.endswith(".maegz"):
                    outName = os.path.splitext(os.path.basename(fnSmall))[0]
                    args = ' -i "{}" --outputName {} --outputDir "{}"'.format(fnSmall, outName,
                                                                       os.path.abspath(self._getExtraPath()))
                    Plugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=self._getExtraPath())
                else:
                    copyFile(fnSmall, self._getExtraPath(os.path.basename(fnSmall)))

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),suffix='SmallMols')
        for fnSmall in glob.glob(self._getExtraPath("*")):
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
            smallMolecule.setMolName(os.path.splitext(os.path.basename(fnSmall))[0])

            if self.draw: # costly
                if not fnSmall.endswith('.mae') and not fnSmall.endswith('.maegz'):
                    fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                    fnOut = self._getExtraPath("%s.png" % fnRoot)
                    args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (fnSmall, fnOut)
                    try:
                        Plugin.runRDKit(self, "python3", args)
                        smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                    except:
                        smallMolecule._PDBLigandImage = pwobj.String("Not available")

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

    def parseMultiPDBFile(self, file, outDir=None):
      if not outDir:
          outDir = self._getExtraPath()

      with open(fnSmall) as fIn:
          pdbs = fIn.read().split('\nEND\n')[:-1]

      names = []
      for i, pdb in enumerate(pdbs):
        molName = "molecule_%s" % (i+1)
        for line in pdb.strip().split('\n'):
          if line.startswith('COMPND'):
            preName = '_'.join(line.split()[1:])
            if not preName in names:
              molName = preName
              names.append(molName)
              break

        with open(os.path.join(outDir, molName + '.pdb'), 'w') as f:
          f.write(pdb)