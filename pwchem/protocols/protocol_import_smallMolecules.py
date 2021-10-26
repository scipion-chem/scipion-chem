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

import glob
import os


from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem import Plugin

class ProtChemImportSmallMolecules(EMProtocol):
    """Import small molecules from a directory. Each molecule should be in a separate file.
       Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
       are accepted.
    """
    _label = 'import small mols'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('multiple', BooleanParam, default=True,
                      label='Each file is a molecule')
        form.addParam('filesPath', PathParam, condition='multiple',
                      label="Files directory",
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
        form.addParam('filesPattern', StringParam,  condition='multiple',
                      label='Pattern',
                      default="*",
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "You may create small molecules from Smiles (.smi), Tripos Mol2 (.mol2), "
                           "SDF (.sdf), Maestro (.mae, .maegz), or PDB blocks (.pdb)")

        form.addParam('filePath', PathParam, condition='not multiple',
                      label='File', help='Allowed formats: \n '
                                         ' - CSV smiles (ID, compound; this is downloaded from ZINC) \n'
                                         ' - Mol2 (Multiple mol2 file) \n'
                                         ' - SDF (Multiple sdf file)')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):

        if not self.multiple.get():

            fnSmall = self._getExtraPath(os.path.split(self.filePath.get())[1])
            copyFile(self.filePath.get(), fnSmall)

            if fnSmall.endswith(".mol2"): # Multiple mol2
                with open(fnSmall) as f:
                    lines = f.readlines()
                    i = 0
                    for line in lines:
                        if line.startswith("#"):
                            continue

                        if line.startswith("@<TRIPOS>MOLECULE"):
                            zincID = lines[i+1].split()[0]
                            fname_small = self._getExtraPath("%s.mol2" %zincID)
                            f_small = open(fname_small, 'w+')
                            f_small.write(line)

                        f_small.write(line)
                        try:
                            if lines[i+1].startswith("@<TRIPOS>MOLECULE"):
                                f_small.close()
                        except:
                            f_small.close()

                        i += 1

                os.remove(fnSmall)


            elif fnSmall.endswith(".sdf"): # Multiple sdf
                with open(fnSmall) as f:
                    lines = f.readlines()
                    lines2write = []
                    i = 0
                    for line in lines:
                        lines2write.append(i)

                        if line.startswith("$$$$"):
                            zincID = lines[i-5].split()[0]
                            fname_small = self._getExtraPath("%s.sdf" % zincID)
                            f_small = open(fname_small, 'w+')
                            for l in lines2write:
                                f_small.write(lines[l])

                            lines2write = []
                            f_small.close()

                        i += 1

                os.remove(fnSmall)

            else:
                fh = open(self.filePath.get())
                for line in fh.readlines():
                    tokens = line.split(',')
                    if len(tokens)==2:
                        fhSmile = open(self._getExtraPath(tokens[0].strip()+".smi"),'w')
                        fhSmile.write(tokens[1].strip()+"\n")
                        fhSmile.close()
                fh.close()

        else:
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnSmall = self._getExtraPath(os.path.split(filename)[1])
                copyFile(filename, fnSmall)

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),suffix='SmallMols')


        for fnSmall in glob.glob(self._getExtraPath("*")):
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)

            if len(os.listdir(self._getExtraPath())) <= 100: # costly
                if not fnSmall.endswith('.mae') or not fnSmall.endswith('.maegz'):
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
