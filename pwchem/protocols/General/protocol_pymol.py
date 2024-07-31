# **************************************************************************
# *
# * Authors:     Marta Martinez (mmmtnez@cnb.csic.es)
# *              Roberto Marabini (roberto@cnb.csic.es)
# *              James Krieger (jmkrieger@cnb.csic.es)
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

import os

from pyworkflow import VERSION_3_0

try:
    from pwem.objects import AtomStruct
except ImportError:
    from pwem.objects import PdbFile as AtomStruct

from pwem.protocols import EMProtocol

from pyworkflow.protocol.params import (MultiPointerParam,
                                        StringParam)
from pyworkflow.utils.properties import Message

from pwchem import Plugin
from pwchem.constants import PYMOL_DIC

class ProtPymolOperate(EMProtocol):
    """This protocol provides access to Pymol and allows one to save the result 
    in the Scipion framework. Files saved in the extra path will be registered as output"""
    _version = VERSION_3_0
    _label = 'Pymol operate'

    # --------------------------- DEFINE param functions --------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputPdbFiles', MultiPointerParam,
                      pointerClass="AtomStruct", allowsNull=True,
                      label='Other atomic structures',
                      help="In case you need to load more PDBx/mmCIF files, "
                           "you can load them here and save them after "
                           "operating with them.")
        form.addParam('extraCommands', StringParam,
                      default='',
                      condition='False',
                      label='Extra commands for chimera viewer',
                      help="Add extra commands in cmd file. Use for testing")

        return form  # DO NOT remove this return

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('runPymolStep')
        self._insertFunctionStep('createOutput')

    # --------------------------- STEPS functions ---------------------------

    def runPymolStep(self):
        # building script file
        pymolScriptFileName = 'script.pml'
        f = open(self._getTmpPath(pymolScriptFileName), "w")
        #f.write('from pymol import cmd\n')

        pdbModelCounter = 1
        if hasattr(self, 'inputPdbFiles'):
            for pdb in self.inputPdbFiles:
                pdbModelCounter += 1
                f.write("load %s\n" % os.path.abspath(pdb.get(
                ).getFileName()))

        # run the text:
        _pymolScriptFileName = os.path.abspath(
            self._getTmpPath(pymolScriptFileName))
        if len(self.extraCommands.get()) > 2:
            f.write(self.extraCommands.get())
        else:
            args = " " + _pymolScriptFileName

        f.close()

        self._log.info('Launching: ' + self._getPymol() + ' ' + args)

        # run in the background
        cwd = os.path.abspath(self._getExtraPath())
        self.runJob(self._getPymol(), args, cwd=cwd)

    def createOutput(self):
        """ Copy the PDB structure and register the output object.
        """
        # Check vol and pdb files
        directory = self._getExtraPath()
        for filename in sorted(os.listdir(directory)):
            if filename.endswith(".pdb") or filename.endswith(".cif"):
                path = os.path.join(directory, filename)
                pdb = AtomStruct()
                pdb.setFileName(path)
                if filename.endswith(".cif"):
                    keyword = filename.split(".cif")[0].replace(".","_")
                else:
                    keyword = filename.split(".pdb")[0].replace(".", "_")
                kwargs = {keyword: pdb}
                self._defineOutputs(**kwargs)

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
            summary.append("Produced files:")
            for filename in sorted(os.listdir(directory)):
                if filename.endswith(".pdb"):
                    summary.append(filename)
            summary.append("we have some result")
        else:
            summary.append(Message.TEXT_NO_OUTPUT_FILES)
        return summary
    
    def _getPymol(self):
        return Plugin.getProgramHome(PYMOL_DIC, 'pymol/bin/pymol')
