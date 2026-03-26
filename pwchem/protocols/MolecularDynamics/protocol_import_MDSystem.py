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
    """Import a Molecular Dynamics system from a system file and a topology file. Trajectory file and extra files are optional.
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
        form.addParam('ligID', params.StringParam, label='Has ligand?', condition='hasLig',
                       help='Does the system has a ligand?', default='LIG')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        obj_type = self.getEnumText('outputType')
        fnTopo = self._getPath(os.path.split(self.inputTopology.get())[1])
        copyFile(self.inputTopology.get(), fnTopo)
        fnCoords = self._getPath(os.path.split(self.inputCoords.get())[1])
        copyFile(self.inputCoords.get(), fnCoords)

        if obj_type == 'MDSystem':
            outSystem = MDSystem()
        elif obj_type == 'OpenMMSystem':
            from openmm.objects import OpenMMSystem
            outSystem = OpenMMSystem()
            fnSerie = self._getPath(os.path.split(self.inputSerie.get())[1])
            copyFile(self.inputSerie.get(), fnSerie)
            outSystem.setSerieFile(fnSerie)
            outSystem.setCifFile(fnCoords)
        elif obj_type == 'AmberSystem':
            from amber.objects import AmberSystem
            outSystem = AmberSystem()
            fnCrd = self._getPath(os.path.split(self.inputCrd.get())[1])
            copyFile(self.inputCrd.get(), fnCrd)
            outSystem.setCrdFile(fnCrd)
        elif obj_type == 'GromacsSystem':
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
