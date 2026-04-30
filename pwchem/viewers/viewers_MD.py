# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params

# pwchem imports
from .. import Plugin
from ..constants import MDTRAJ_DIC, TCL_MD_STR, PML_MD_STR, PML_MD_STR_AMBER, TCL_MD_LIG_STR
from ..viewers import PyMolViewer, PyMolView, VmdViewPopen
from ..objects import MDSystem

_ANAL_RMSD      = 0
_ANAL_RMSF      = 1
_ANAL_RG        = 2
_ANAL_SS        = 3
_ANAL_SASA      = 4
_ANAL_PCA       = 5
_ANAL_DISTANCE  = 6

_ANAL_CHOICES = ['RMSD', 'RMSF', 'Rg', 'Secondary Structure', 'SASA', 'PCA', 'Atom distance']

# Analyses that use the atom-selection dropdown
_USES_SEL_ATOMS   = [_ANAL_RMSD, _ANAL_RMSF]

class MDSystemViewer(pwviewer.Viewer):
  _label = 'Viewer Molecular Dynamics system'
  _environments = [pwviewer.DESKTOP_TKINTER]
  _targets = []

  def _visualize(self, obj, onlySystem=False, trjFile=None, **kwargs):
    systemFile = os.path.abspath(obj.getSystemFile())
    if not trjFile:
        trjFile = obj.hasTrajectory()

    if not trjFile or onlySystem:
        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(systemFile, cwd=os.path.dirname(systemFile))

    else:
        trjFile = os.path.abspath(obj.getTrajectoryFile())
        topoFile = os.path.abspath(obj.getTopologyFile())
        if topoFile.lower().endswith(('.top', '.tpr')):
            topoFile = obj.getSystemFile() # needed for gromacs topology
        outPml = os.path.join(os.path.dirname(trjFile), 'pymolSimulation.pml')
        _, trjExt = os.path.splitext(trjFile)
        if trjExt in ['.nc', '.netcdf']:
            template = PML_MD_STR_AMBER
        else:
            template = PML_MD_STR
        with open(outPml, 'w') as f:
          f.write(template.format(os.path.abspath(topoFile),
                                    os.path.abspath(trjFile)))

        return [PyMolView(os.path.abspath(outPml), cwd=os.path.dirname(trjFile))]


class MDSystemPViewer(pwviewer.ProtocolViewer):
    """Visualize the output of a Molecular Dynamics simulation."""
    _label = 'Viewer Molecular Dynamics System'
    _targets = [MDSystem]
    _mdtrajScript   = 'mdtraj_analysis.py'
    _prolifViewScript = 'prolif_viewer.py'

    def __init__(self, **args):
        super().__init__(**args)

    # ------------------------------------------------------------------
    # Form definition
    # ------------------------------------------------------------------

    def _defineSimParams(self, form):
        group = form.addGroup('Open MD simulation')
        group.addParam('displayMdPymol', params.LabelParam,
                       label='Display trajectory with PyMol: ',
                       help='Display trajectory with PyMol.\n'
                            'Protein as NewCartoon, waters as sticks.')
        group.addParam('displayMdVMD', params.LabelParam,
                       label='Display trajectory with VMD: ',
                       help='Display trajectory with VMD.\n'
                            'Protein as NewCartoon, waters as dots.')

    def _defineMDTrajParams(self, form):
        form.addSection(label='Trajectory analysis')
        group = form.addGroup('MDTraj analysis')

        group.addParam('mdAnalChoices', params.EnumParam,
                       label='Analysis type: ',
                       choices=_ANAL_CHOICES, default=_ANAL_RMSD,
                       help='Select the MDTraj analysis to run.\n'
                            'Relevant parameters will appear below.')

        # ── Shared: atom selection (RMSD / RMSF) ────────
        group.addParam('selAtoms', params.EnumParam,
                       label='Selection of atoms: ', default=0,
                       choices=['Protein', 'Backbone', 'CA', 'Sidechain', 'Ligand'],
                       condition=f'mdAnalChoices in {_USES_SEL_ATOMS}',
                       help='Atom selection used in the analysis.')

        # ── Shared: heavy atoms (RMSD / RMSF) ──────────────────
        group.addParam('heavyAtoms', params.BooleanParam,
                       label='Use only heavy atoms: ', default=True,
                       condition=f'mdAnalChoices in {_USES_SEL_ATOMS}',
                       help='Restrict analysis to non-hydrogen atoms.')

        # ── Secondary Structure sub-options ─────────────────────────
        group.addParam('ssDisplayType', params.EnumParam, default=0,
                       label='Display mode: ',
                       choices=['Per residue', 'Per frame', 'Heatmap'],
                       condition=f'mdAnalChoices == {_ANAL_SS}',
                       help='How to present the DSSP secondary-structure results.')

        # ── SASA sub-options ─────────────────────────────────────────
        group.addParam('sasaScope', params.EnumParam, default=0,
                       label='SASA scope: ',
                       choices=['All', 'Ligand'],
                       condition=f'mdAnalChoices == {_ANAL_SASA}',
                       help='Plot SASA for all atoms or the ligand only.\n'
                            'Computation time can be long.')

        # ── PCA sub-options ──────────────────────────────────────────
        group.addParam('pcaType', params.EnumParam, default=0,
                       label='PCA input: ',
                       choices=['Coordinates', 'Pairwise distances'],
                       condition=f'mdAnalChoices == {_ANAL_PCA}',
                       help='Run PCA on Cartesian coordinates or backbone pairwise distances.')

        # ── Distance sub-options ─────────────────────────────────────
        distLine = group.addLine('Atom indices: ',
                                 condition=f'mdAnalChoices == {_ANAL_DISTANCE}',
                                 help='Indices of the two atoms (from the PDB) whose distance '
                                      'will be tracked across frames.')
        distLine.addParam('atom1', params.IntParam, allowsNull=True, label='Atom 1: ')
        distLine.addParam('atom2', params.IntParam, allowsNull=True, label='Atom 2: ')

        # ── Run button ────────────────────────────────────────
        group.addParam('displayMDTrajAnalysis', params.LabelParam,
                       label='Display MDTraj analysis: ',
                       help='Run and display the selected analysis.')

        if self.getMDSystem().getProlifFile():
            self._defineProlifParams(form)

    def _defineProlifParams(self, form):
        form.addSection('Receptor-ligand interactions')
        group = form.addGroup('ProLIF analysis')
        group.addParam('displayFingerprint', params.LabelParam, label='Show interaction fingerprint: ',
                       help='Display interaction fingerprint generated by ProLIF')
        group.addParam('displayInterNetwork', params.LabelParam, label='Display interaction network: ',
                      help='Display ligand interaction network generated by ProLIF (web browser is opened).\n'
                           'Only the interactions that are present in more than 30% of the frames are shown.')
        group.addParam('displayProlifMatrix', params.LabelParam, label='Display similarity matrix: ',
                      help='Display Tanimoto similarity matrix calculated using ligand-target interaction fingerprint.')

    def _defineParams(self, form):
        form.addSection(label='Visualization of MD System')
        group = form.addGroup('Open MD system')
        group.addParam('displayPymol', params.LabelParam,
                       label='Open system in PyMol: ',
                       help='Display the system in the PyMol GUI.')

        if self.getMDSystem().hasTrajectory():
            self._defineSimParams(form)
            self._defineMDTrajParams(form)

    def getMDSystem(self, objType=MDSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        return self.protocol.outputSystem

    def _getVisualizeDict(self):
        return {
            'displayPymol':           self._showPymol,
            'displayMdPymol':         self._showMdPymol,
            'displayMdVMD':           self._showMdVMD,
            'displayMDTrajAnalysis':  self._showMDTrajAnalysis,
            'displayFingerprint':     self._showProlifFp,
            'displayInterNetwork':    self._showProlifNetwork,
            'displayProlifMatrix':    self._showProlifMatrix,
        }

    def _showMDTrajAnalysis(self, paramName=None):
        dispatch = {
            _ANAL_RMSD:     self._showMDTrajRMSDRMSF,
            _ANAL_RMSF:     self._showMDTrajRMSDRMSF,
            _ANAL_RG:       self._showMDTrajRGAnalysis,
            _ANAL_SS:       self._showMDTrajSSAnalysis,
            _ANAL_SASA:     self._showMDTrajSASAAnalysis,
            _ANAL_PCA:      self._showMDTrajPCAAnalysis,
            _ANAL_DISTANCE: self._showDistance,
        }
        return dispatch[self.mdAnalChoices.get()]()

    # ------------------------------------------------------------------
    # Visualize methods
    # ------------------------------------------------------------------

    def _showPymol(self, paramName=None):
        system = self.getMDSystem()
        return MDSystemViewer(project=self.getProject())._visualize(system, onlySystem=True)

    def _showMdPymol(self, paramName=None):
        system = self.getMDSystem()
        return MDSystemViewer(project=self.getProject())._visualize(system)

    def writeTCL(self, outTcl, sysFile, sysExt, sysTrj, trjExt):
        system  = self.getMDSystem()
        vmdStr  = TCL_MD_STR % (sysFile, sysExt, sysTrj, trjExt)
        vmdStr += TCL_MD_LIG_STR.format(system.getLigandID())
        with open(outTcl, 'w') as f:
            f.write(vmdStr)

    def _showMdVMD(self, paramName=None):
        system = self.getMDSystem()
        outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
        sysExt = os.path.splitext(system.getFileName())[1][1:]
        trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
        self.writeTCL(outTcl, system.getFileName(), sysExt,
                      system.getTrajectoryFile(), trjExt)
        return [VmdViewPopen('-e {}'.format(outTcl))]

    def _showMDTrajRMSDRMSF(self, paramName=None):
        """Handles both RMSD and RMSF (distinguished by mdAnalChoices text)."""
        system   = self.getMDSystem()
        analFlag = self.getEnumText('mdAnalChoices').lower()   # 'rmsd' or 'rmsf'
        selAtoms = self.getEnumText('selAtoms')
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -{analFlag} -sa {selAtoms} ')
        if self.heavyAtoms.get():
            args += '-ha '
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajRGAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -rg ')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSSAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        ssFlags  = ['--per-residue', '--per-frame', '--heatmap']
        flag     = ssFlags[self.ssDisplayType.get()]
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} {flag} ')
        Plugin.runScript(self, 'mdtraj_SS.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSASAAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        selAtoms = self.getEnumText('sasaScope')
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -sa {selAtoms} -sasa ')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajPCAAnalysis(self, paramName=None):
        system   = self.getMDSystem()
        pcaFlags = ['--pca-coord', '--pca-dist']
        flag     = pcaFlags[self.pcaType.get()]
        args = (f'-i {system.getFileName()} -t {system.getTrajectoryFile()} {flag} ')
        Plugin.runScript(self, 'mdtraj_PCA.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showDistance(self, paramName=None):
        system       = self.getMDSystem()
        atom1, atom2 = self.atom1.get(), self.atom2.get()
        args = (f' -distance -i {system.getFileName()} -t {system.getTrajectoryFile()} '
                f'-o {system.getSystemName()} -a1 {atom1} -a2 {atom2}')
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifFp(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode barcode',
                         env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifNetwork(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode network',
                         env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifMatrix(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        Plugin.runScript(self, self._prolifViewScript, f'"{fpPkl}" --mode matrix',
                         env=MDTRAJ_DIC, popen=True, wait=False)


