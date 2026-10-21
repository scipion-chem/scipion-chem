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
import os, json, shutil
import pickle as pkl

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.utils import cleanPDB, getBaseName, removeStartWithLines
from pwchem import Plugin as pwchemPlugin
from pwchem.objects import MDSystem
from pwchem.constants import MDTRAJ_DIC

class ProtocolProlif(EMProtocol):
    """
    AI Generated:

    Protocol for Protein-Ligand interaction analysis using ProLIF.

    This protocol runs ProLIF on a Molecular Dynamics (MD) simulation to compute
    interaction fingerprints between a protein and a ligand (or water bridges)
    across trajectory frames.

    Overview
    --------
    The protocol performs trajectory-based interaction analysis using:

    - Protein–ligand contact fingerprints
    - Optional water-bridge interaction analysis
    - Frame-sampled evaluation of interactions

    Input
    -----
    inputMDSystem:
        Molecular Dynamics system containing:
        - topology file
        - trajectory file
        - optional ligand identifier

    type:
        Type of interaction analysis:
        - Protein-Ligand
        - Water bridges

    frameNum:
        Frame sampling frequency.
        Example: 10 → analyze every 10 frames.

    Workflow
    --------
    1. Extract MD system components:
       - trajectory file
       - topology file (or system file for GROMACS)

    2. Configure ligand:
       - Uses ligand ID if available
       - Defaults to 'LIG'

    3. Build ProLIF execution command:
       - Input topology
       - Trajectory
       - Output directory
       - System name
       - Ligand name
       - Frame sampling frequency

    4. Optional analysis mode:
       - Adds "-wb" flag for water-bridge analysis

    5. Execute ProLIF script:
       - Run via pwchem Plugin
       - Uses MDTraj environment

    Output
    ------
    outputSystem:
        Cloned MD system enriched with:
        - ProLIF fingerprint file (.pkl)

    Output file:
        <system_name>_fingerprint.pkl

    Internal Methods
    ----------------
    runProlif():
        Executes ProLIF analysis over the MD trajectory.

    createOutputStep():
        Creates output system and attaches fingerprint file.

    _summary():
        Reports analysis type, frame frequency, and output fingerprint.

    Key Features
    ------------
    - Protein-ligand interaction fingerprints
    - Water bridge detection option
    - Frame-sampled trajectory analysis
    - Integration with MD system objects
    - Automated output attachment

    Notes
    -----
    - Requires valid trajectory + topology MD system
    - Uses MDTraj-based execution environment
    - Output is stored as serialized pickle fingerprint
    """
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
        gmxTopArgs = ''
        if topoFile.lower().endswith(('.top', '.tpr')):
            if topoFile.lower().endswith('.top'):
                # Bonds are read from the .top; the .gro (no bonds) is the fallback if it cannot be read
                gmxTopArgs = self.getGromacsTopArgs(inputSystem)
            topoFile = inputSystem.getSystemFile() # needed for gromacs topology
        outputPath = self._getExtraPath()
        outputName = inputSystem.getSystemName()
        if inputSystem.getLigandID():
            ligName =  inputSystem.getLigandID()
        else:
            ligName = 'LIG'
        args = f'-i {topoFile} -t {trajFile} -o {outputPath} -n {outputName} -l {ligName} -f {self.frameNum.get()}'
        args += gmxTopArgs
        if self.getEnumText('type') == 'Water bridges':
            args += ' -wb'

        args += ' 2>&1'

        pwchemPlugin.runScript(self, 'prolif_analysis.py', args, env=MDTRAJ_DIC, wait=False)

    def getGromacsTopArgs(self, inputSystem):
        '''Script arguments to read the bonds from a GROMACS topology: a copy of the .top together with the files it
        #includes by relative path (e.g. the ligand itp, missing next to the .top of older runs) and the directory of
        the force fields shipped with GROMACS'''
        topFile = os.path.abspath(inputSystem.getTopologyFile())
        workDir = os.path.abspath(self._getExtraPath('topology'))
        os.makedirs(workDir, exist_ok=True)
        shutil.copy(topFile, workDir)

        searchDirs = [os.path.dirname(topFile)]
        if inputSystem.getLigTopologyFile():
            searchDirs.append(os.path.dirname(os.path.abspath(inputSystem.getLigTopologyFile())))
        for incFile in self.getRelativeIncludes(topFile):
            for sDir in searchDirs:
                if os.path.isfile(os.path.join(sDir, incFile)):
                    os.makedirs(os.path.dirname(os.path.join(workDir, incFile)), exist_ok=True)
                    shutil.copy(os.path.join(sDir, incFile), os.path.join(workDir, incFile))
                    break

        args = f' --gmxTop {os.path.join(workDir, os.path.basename(topFile))}'
        gmxTopDir = self.getGromacsTopDir()
        if gmxTopDir:
            args += f' --includeDir {gmxTopDir}'
        return args

    @staticmethod
    def getGromacsTopDir():
        '''Force fields directory of the GROMACS installed by the gromacs plugin, if any'''
        try:
            from gromacs import Plugin as gromacsPlugin
            return gromacsPlugin.getTopDir()
        except Exception:
            return ''

    def createOutputStep(self):
        inputSystem = self.inputMDSystem.get()
        baseName = inputSystem.getSystemName()
        fpPath = self._getExtraPath(f'{baseName}_fingerprint.pkl')
        outSystem = inputSystem.clone()
        outSystem.setProlifFile(fpPath)
        self._defineOutputs(outputSystem=outSystem)

    def _summary(self):
        summary = []
        if hasattr(self, 'outputSystem'):
            outSystem = self.outputSystem
            summary.append('ProLIF analysis performed')
            summary.append(f'Analysis type  : {self.getEnumText("type")}')
            if outSystem.getProlifFile():
                summary.append(f'Fingerprint    : {outSystem.getProlifFile()}')
            summary.append(f'Fingerprint frequency: every {self.frameNum.get()} frames')
        return summary

    def getRelativeIncludes(self, topFile):
        '''Files #included by relative path in a GROMACS topology'''
        includes = []
        with open(topFile) as f:
            for line in f:
                if line.strip().startswith('#include'):
                    incFile = line.split('"')[1] if '"' in line else line.split()[1]
                    if not os.path.isabs(incFile):
                        includes.append(incFile)
        return includes
