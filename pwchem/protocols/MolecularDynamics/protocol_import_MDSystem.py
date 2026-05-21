# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

import glob, os, urllib

from parmed.tools.actions import outPDB
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params

from pwem.protocols import EMProtocol
from pwchem.objects import MDSystem
class ProtocolImportMDSystem(EMProtocol):
    """
    Protocol to import a Molecular Dynamics (MD) system into the Scipion/pwchem framework.

    This protocol creates a unified MDSystem object (or backend-specific variant such as OpenMM,
    Amber, or Gromacs systems) from coordinate, topology, and optional trajectory files.

    It acts as an interface layer between external MD simulation outputs and the Scipion data model.

    Overview
    --------
    The protocol performs the following steps:

    1. Input file preparation
       - Copies input coordinate and topology files into the protocol workspace
       - Optionally copies trajectory and system-specific files

    2. System instantiation
       - Creates the appropriate MD system object depending on selected backend:
         * MDSystem (generic)
         * OpenMMSystem
         * AmberSystem
         * GromacsSystem

    3. System configuration
       - Registers file paths inside the MD object
       - Optionally assigns ligand information

    4. Output definition
       - Produces a single structured MDSystem-like object

    Input
    -----
    outputType:
        Type of MD system to generate:
        - MDSystem (generic representation)
        - OpenMMSystem (OpenMM-specific system with XML support)
        - AmberSystem (Amber RST7-based system)
        - GromacsSystem (GROMACS-based system)

    inputCoords:
        Coordinate file of the system.
        Supported formats: PDB, CIF, GRO.

    inputTopology:
        Topology file defining the system structure.
        Supported formats: PDB, TOP, TPR, PRMTOP, PARM7, CMS, XML.

    inputSerie:
        (OpenMM only) XML series file required for OpenMM systems.

    inputCrd:
        (Amber only) Initial coordinate file in RST7 format.

    addTraj:
        Boolean flag indicating whether a trajectory file is provided.

    inputTrajectory:
        Optional trajectory file.
        Supported formats: XTC, NC, NETCDF, DCD.

    hasLig:
        Boolean flag indicating presence of a ligand in the system.

    ligID:
        Identifier of the ligand residue in the structure (e.g. "LIG").

    Workflow
    --------
    1. Copy input files to protocol working directory.
    2. Determine system type from outputType.
    3. Instantiate appropriate MD system class.
    4. Load backend-specific files:
       - OpenMM → XML + CIF
       - Amber → RST7
       - Gromacs → standard files
    5. Assign common system files:
       - Coordinates file
       - Topology file
    6. Optionally assign trajectory file.
    7. Optionally assign ligand ID.
    8. Export final MD system object.

    Output
    ------
    outputSystem:
        MDSystem (or subclass) containing:
        - system coordinates file
        - topology file
        - optional trajectory file
        - optional ligand annotation

    Supported backends
    ------------------
    MDSystem:
        Generic representation without backend dependency.

    OpenMMSystem:
        Includes OpenMM XML serialization and CIF linkage.

    AmberSystem:
        Supports Amber coordinate (RST7) integration.

    GromacsSystem:
        Minimal wrapper for GROMACS simulation files.

    Summary
    -------
    This protocol provides a standardized way to import molecular dynamics systems
    from multiple simulation engines into Scipion, enabling downstream analysis,
    visualization, and protocol chaining in a unified format.

    Key features:
    - Multi-backend MD support
    - Automatic file management inside protocol workspace
    - Optional ligand annotation
    - Optional trajectory integration
    """
    _label = 'Import Molecular Dynamics System'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('outputType', params.EnumParam, label='Output object type', default=0,
                      choices=['MDSystem', 'OpenMMSystem', 'AmberSystem', 'GromacsSystem'],
                      help='What type of system to produce. Base molecular dynamics system (MDSystem) '
                           'or Plugin specific system.')
        form.addParam('inputCoords', params.PathParam, label='Coordinate file',
                      help='Coordinates file - Formats can be PDB, CIF and GRO')
        form.addParam('inputSerie', params.PathParam, label='OpenMM Serie file',
                      help='OpenMM serie file - Format XML needed.',
                      condition='outputType == 1')
        form.addParam('inputCrd', params.PathParam, label='Initial Amber coordinate file',
                      help='Initial Amber coordinate file - Format RST7 needed.',
                      condition='outputType == 2')
        form.addParam('inputTopology', params.PathParam, label='Topology file',
                       help='Topology file - Formats can be PDB, TOP, TPR, PRMTOP, PARM7, CMS and XML')
        form.addParam('addTraj', params.BooleanParam, default=False, label='Trajectory file')
        form.addParam('inputTrajectory', params.PathParam, label='Trajectory file', condition='addTraj',
                       help='Trajectory file - Formats can be XTC, NC, NETCDF and DCD')

        form.addParam('hasLig', params.BooleanParam, label='Has ligand?', default=False,
                       help='Does the system has a ligand?')
        form.addParam('ligID', params.StringParam, label='Ligand name in PDB', condition='hasLig',
                       help='Does the system has a ligand?', default='LIG')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.importStep)

    def importStep(self):
        objType = self.getEnumText('outputType')
        fnTopo = self._getPath(os.path.split(self.inputTopology.get())[1])
        copyFile(self.inputTopology.get(), fnTopo)
        fnCoords = self._getPath(os.path.split(self.inputCoords.get())[1])
        copyFile(self.inputCoords.get(), fnCoords)

        if objType == 'MDSystem':
            outSystem = MDSystem()
        elif objType == 'OpenMMSystem':
            from openmm.objects import OpenMMSystem
            outSystem = OpenMMSystem()
            fnSerie = self._getPath(os.path.split(self.inputSerie.get())[1])
            copyFile(self.inputSerie.get(), fnSerie)
            outSystem.setSerieFile(fnSerie)
            outSystem.setCifFile(fnCoords)
        elif objType == 'AmberSystem':
            from amber.objects import AmberSystem
            outSystem = AmberSystem()
            fnCrd = self._getPath(os.path.split(self.inputCrd.get())[1])
            copyFile(self.inputCrd.get(), fnCrd)
            outSystem.setCrdFile(fnCrd)
        elif objType == 'GromacsSystem':
            from gromacs.objects import GromacsSystem
            outSystem = GromacsSystem()
        outSystem.setSystemFile(fnCoords)
        outSystem.setTopologyFile(fnTopo)
        if self.inputTrajectory.get():
            fnTraj = self._getPath(os.path.split(self.inputTrajectory.get())[1])
            copyFile(self.inputTrajectory.get(), fnTraj)
            outSystem.setTrajectoryFile(fnTraj)
        if self.hasLig.get():
            outSystem.setLigandID(self.ligID.get())

        self._defineOutputs(outputSystem=outSystem)

    # --------------------------- INFO --------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputSystem'):
            outSystem = self.outputSystem
            summary.append(f'Output type    : {self.getEnumText("outputType")}')
            summary.append(f'Coordinate file: {outSystem.getSystemFile()}')
            summary.append(f'Topology file  : {outSystem.getTopologyFile()}')
            if outSystem.hasTrajectory():
                summary.append(f'Trajectory file: {outSystem.getTrajectoryFile()}')
            else:
                summary.append('Trajectory file: not provided')
            if self.hasLig.get():
                summary.append(f'Ligand ID      : {self.ligID.get()}')
        return summary
