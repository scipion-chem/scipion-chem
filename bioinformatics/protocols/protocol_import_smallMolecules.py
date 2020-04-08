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

import glob
import lxml.etree as ET
import os
import sys
import urllib.request

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam
from bioinformatics.objects import SmallMolecule, SetOfSmallMolecules
from bioinformatics import Plugin

class ProtBioinformaticsImportSmallMolecules(EMProtocol):
    """Import small molecules from a directory. Each molecule should be in a separate file.
       Smiles (.smi), SDF (.sdf), Tripos Mol2 (.mol2), Maestro (.mae, .maegz), PDB blocks (.pdb)
       are accepted.
    """
    _label = 'import small mols'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('filesPath', PathParam,
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
        form.addParam('filesPattern', StringParam,
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

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(),suffix='SmallMols')
        for filename in glob.glob(os.path.join(self.filesPath.get(),self.filesPattern.get())):
            fnSmall = self._getExtraPath(os.path.split(filename)[1])
            copyFile(filename,fnSmall)
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)

            fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
            fnOut = self._getExtraPath("%s.png" % fnRoot)
            args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (fnSmall, fnOut)
            try:
                Plugin.runRDKit(self, "python3", args)
                smallMolecule._PDBLigandImage = pwobj.String(fnOut)
            except:
                smallMolecule._PDBLigandImage = pwobj.String("Not available")

            outputSmallMolecules.append(smallMolecule)
        self._defineOutputs(outputSmallMols=outputSmallMolecules)
