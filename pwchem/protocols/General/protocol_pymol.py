# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              James Krieger (jmkrieger@cnb.csic.es)
# *              Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * L'Institut de genetique et de biologie moleculaire et cellulaire (IGBMC)
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

import inspect
import os

from pyworkflow import VERSION_3_0

try:
    from pwem.objects import AtomStruct
except ImportError:
    from pwem.objects import PdbFile as AtomStruct

from pwem.protocols import EMProtocol

from pyworkflow.protocol import params
from pyworkflow.utils.properties import Message

from pwchem import Plugin
from pwchem.constants import OPENBABEL_DIC

pymolScriptFileName = 'script.pml'

class ProtPymolOperate(EMProtocol):
    """This protocol provides access to Pymol and allows one to save the result 
    in the Scipion framework. Files saved in the extra path will be registered as output"""
    _version = VERSION_3_0
    _label = 'Pymol operate'

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPdbFiles', params.MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True, label='Input atomic structures: ',
                      help="In case you need to load more PDBx/mmCIF files, you can load them here and save them after "
                           "operating with them.")
        form.addParam('extraCommands', params.StringParam, default='', condition='False',
                      label='Extra commands for pymol viewer', expertLevel=params.LEVEL_ADVANCED,
                      help="Add extra commands in cmd file. Use for testing")

        return form  # DO NOT remove this return

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runPymolStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------

    def runPymolStep(self):
        # building script file
        pymolScript = self.writePymolScript()

        # run in the background
        self._log.info('Launching: ' + self._getPymol() + pymolScript)
        self.runJob(self._getPymol(), pymolScript, cwd=self._getExtraPath())

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        i = 1
        oDir = self._getExtraPath()
        for filename in sorted(os.listdir(oDir)):
            if os.path.splitext(filename)[1] in ['.pdb', '.cif']:
                path = os.path.join(oDir, filename)
                pdb = AtomStruct(path)

                kwargs = {f'outputStructure_{i}': pdb}
                self._defineOutputs(**kwargs)
                i += 1


    # --------------------------- INFO functions ----------------------------
    def _validate(self):
        errors = []
        # Check that the program exists
        program = self._getPymol()
        if not os.path.exists(program):
            errors.append("Binary '%s' does not exists.\n" % program)
        return errors

    def _summary(self):
        summary = []
        if self.getOutputsSize() > 0:
            directory = self._getExtraPath()
            summary.append("Produced files: ")
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".pdb") or filename.endswith(".cif"):
                    summary.append(f'\t- {filename}')
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary
    
    def _getPymol(self):
        return Plugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')

    def getPymolScript(self):
        return os.path.abspath(self._getTmpPath(pymolScriptFileName))

    def writePymolScript(self):
        scriptFn = self.getPymolScript()
        with open(scriptFn, 'w') as f:
            for pdb in self.inputPdbFiles:
                pdbFile = os.path.abspath(pdb.get().getFileName())
                f.write(f"load {pdbFile}\n")

            if self.extraCommands.get().strip():
                f.write(self.extraCommands.get().strip())
        return scriptFn

