# **************************************************************************
# *
# * Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
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
import os, json
import pickle as pkl

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.utils import cleanPDB, getBaseName, removeStartWithLines
from pwchem import Plugin as pwchemPlugin
from pwchem.objects import MDSystem
from pwchem.constants import MDTRAJ_DIC

class ProtocolProlif(EMProtocol):
    """Run ProLIF Target-Ligand trajectory analysis"""
    _label = 'ProLIF analysis'
    _program = "ProLIF"

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMDSystem', params.PointerParam, pointerClass="MDSystem",
                      label='MD simulation: ', allowsNull=False,
                      help='It must be a Molecular Dynamics System with trajectory.')
        form.addParam('type', params.EnumParam, choices=['Protein-Ligand', 'Water bridges'],
                      label='Type of analysis: ', default=0,
                      help='Analyize protein-ligand interactions over each frame or the water bridges.')
        form.addParam('frameNum', params.IntParam, label='Frame frequency', default=10,
                      help='Every how many frames interactions are analyzed.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.runProlif)
        self._insertFunctionStep(self.createOutputStep)

    def runProlif(self):
        inputSystem = self.inputMDSystem.get()
        trajFile = inputSystem.getTrajectoryFile()
        topoFile = inputSystem.getTopologyFile()
        outputPath = self._getExtraPath()
        outputName = inputSystem.getSystemName()
        args = f'-i {topoFile} -t {trajFile} -o {outputPath} -n {outputName} -f {self.frameNum.get()}'
        if self.getEnumText('type') == 'Water bridges':
            args += ' -wb'

        args += ' 2>&1'

        pwchemPlugin.runScript(self, 'prolif_analysis.py', args, env=MDTRAJ_DIC, wait=False)

    def createOutputStep(self):
        inputSystem = self.inputMDSystem.get()
        baseName = inputSystem.getSystemName()
        fpPath = self._getExtraPath(f'{baseName}_fingerprint.pkl')
        outSystem = inputSystem.clone()
        outSystem.setProlifFile(fpPath)
        self._defineOutputs(outputSystem=outSystem)

