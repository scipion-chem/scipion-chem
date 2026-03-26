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
    """ Visualize the output of Molecular Dynamics simulation """
    _label = 'Viewer Molecular Dynamics System'
    _targets = [MDSystem]
    _mdtrajScript = 'mdtraj_analysis.py'
    _prolifViewScript = 'prolif_viewer.py'

    def __init__(self, **args):
      super().__init__(**args)

    def _defineSimParams(self, form):
        group = form.addGroup('Open MD simulation')
        group.addParam('displayMdPymol', params.LabelParam,
                       label='Display trajectory with PyMol: ',
                       help='Display trajectory with Pymol. \n'
                            'Protein represented as NewCartoon and waters as sticks'
                       )
        group.addParam('displayMdVMD', params.LabelParam,
                       label='Display trajectory with VMD: ',
                       help='Display trajectory with VMD. \n'
                            'Protein represented as NewCartoon and waters as dots')

    def _defineMDTrajParams(self, form):
        form.addSection(label='Trajectory analysis')
        group = form.addGroup('MDTraj analysis')
        group.addParam('mdAnalChoices', params.EnumParam, label='Display trajectory analysis: ',
                       choices=['RMSD', 'RMSF'], default=0,
                       help='Uses MDTraj to display this analysis of the trajectory')
        group.addParam('selAtoms', params.EnumParam, label='Selection of atoms: ', default=0,
                       choices=['Protein', 'Backbone', 'CA', 'Sidechain', 'Ligand'],
                       help='Selection of atoms to use in the analysis')
        group.addParam('heavyAtoms', params.BooleanParam, label='Use only heavy atoms: ', default=True,
                       help='Uses only heavy atoms in the analysis of the trajectory')
        group.addParam('displayMDTrajAnalysis', params.LabelParam,
                       label='Display MDTraj analysis: ', help='Display MDTraj defined analysis')

        group.addParam('displayMDTrajRGAnalysis', params.LabelParam,
                       label='Display Rg analysis: ',
                       help='Display evolution of the Radius of Gyration over the trajectory ')
        group.addParam('displayMDTrajSSAnalysis', params.EnumParam, default=0,
                       choices=['Per residue', 'Per frame', 'Heatmap'],
                       label='Display Secondary Structure analysis: ',
                       help='Display secondary structure analysis with estimations obtained with DSSP')
        group.addParam('displayMDTrajSASAAnalysis', params.EnumParam, default=0,
                       label='Display SASA analysis: ', choices=['All', 'Ligand'],
                       help='Display evolution of the SASA over the trajectory. Plot SASA of all atoms or only the Ligand atoms'
                            'Computation time is long.')
        group.addParam('displayMDTrajPCAAnalysis', params.EnumParam, default=0,
                       label='Display PCA analysis: ', choices=['Coordinates', 'Pairwise distances'],
                       help='Display a PCA analysis of the coordinates of the trajectory or the pairwise distances of the backbone.')

        distLine = group.addLine('Distance: ', help='Specify 2 atoms of the system PDB file to compute the distance for every frame'
                                ' in the trajectory')
        distLine.addParam('atom1', params.IntParam, allowsNull=True,
                       label='Atom num 1:')
        distLine.addParam('atom2', params.IntParam, allowsNull=True,
                       label='Atom num 2: ')
        group.addParam('displayDistance', params.LabelParam, help='Display distance over time',
                          label='Display distance analysis: ')

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
                     help='Display System in Pymol GUI.')

      if self.getMDSystem().hasTrajectory():
          self._defineSimParams(form)
          self._defineMDTrajParams(form)

    def getMDSystem(self, objType=MDSystem):
        if isinstance(self.protocol, objType):
            return self.protocol
        else:
            return self.protocol.outputSystem

    def _getVisualizeDict(self):
      return {
        'displayPymol': self._showPymol,
        'displayMdPymol': self._showMdPymol,
        'displayMdVMD': self._showMdVMD,

        'displayMDTrajAnalysis': self._showMDTrajAnalysis,
        'displayMDTrajRGAnalysis': self._showMDTrajRGAnalysis,
        'displayMDTrajSSAnalysis': self._showMDTrajSSAnalysis,
        'displayMDTrajSASAAnalysis': self._showMDTrajSASAAnalysis,
        'displayMDTrajPCAAnalysis': self._showMDTrajPCAAnalysis,
        'displayFingerprint': self._showProlifFp,
        'displayInterNetwork': self._showProlifNetwork,
        'displayProlifMatrix': self._showProlifMatrix,
        'displayDistance': self._showDistance
      }

    def _showPymol(self, paramName=None):
      system = self.getMDSystem()
      return MDSystemViewer(project=self.getProject())._visualize(system, onlySystem=True)

    def _showMdPymol(self, paramName=None):
      system = self.getMDSystem()
      return MDSystemViewer(project=self.getProject())._visualize(system)

    def writeTCL(self, outTcl, sysFile, sysExt, sysTrj, trjExt):
      system = self.getMDSystem()
      vmdStr = TCL_MD_STR % (sysFile, sysExt, sysTrj, trjExt)
      vmdStr += TCL_MD_LIG_STR.format(system.getLigandID())
      with open(outTcl, 'w') as f:
        f.write(vmdStr)


    def _showMdVMD(self, paramName=None):
      system = self.getMDSystem()

      outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
      sysExt = os.path.splitext(system.getFileName())[1][1:]
      trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
      self.writeTCL(outTcl, system.getFileName(), sysExt, system.getTrajectoryFile(), trjExt)

      args = '-e {}'.format(outTcl)
      return [VmdViewPopen(args)]

    def _showMDTrajAnalysis(self, paramName=None):
      system = self.getMDSystem()
      selAtoms = self.getEnumText("selAtoms")

      args = f'-i {system.getFileName()} -t {system.getTrajectoryFile()} -o {system.getSystemName()} ' \
             f'-{self.getEnumText("mdAnalChoices").lower()} -sa {selAtoms} '
      if self.heavyAtoms.get():
        args += '-ha '
      Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSSAnalysis(self, paramName=None):
        system = self.getMDSystem()
        showType = self.displayMDTrajSSAnalysis.get()
        ssFlags = ['--per-residue', '--per-frame', '--heatmap']

        args = f'-i {system.getFileName()} -t {system.getTrajectoryFile()} {ssFlags[showType]}'

        Plugin.runScript(self, 'mdtraj_SS.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showDistance(self, paramName=None):
        system = self.getMDSystem()
        atom1, atom2 = self.atom1.get(), self.atom2.get()

        args = f' -distance  -i {system.getFileName()} -t {system.getTrajectoryFile()} -o {system.getSystemName()} ' \
               f' -a1 {atom1} -a2 {atom2}'
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajRGAnalysis(self, paramName=None):
        system = self.getMDSystem()
        selAtoms = self.getEnumText("selAtoms")

        args = f'-i {system.getFileName()} -t {system.getTrajectoryFile()} -o {system.getSystemName()} ' \
               f' -sa {selAtoms} -rg '
        if self.heavyAtoms.get():
            args += '-ha '
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajSASAAnalysis(self, paramName=None):
        system = self.getMDSystem()
        selAtoms = self.getEnumText("displayMDTrajSASAAnalysis")

        args = f'-i {system.getFileName()} -t {system.getTrajectoryFile()} -o {system.getSystemName()} ' \
               f' -sa {selAtoms} -sasa '
        Plugin.runScript(self, self._mdtrajScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showMDTrajPCAAnalysis(self, paramName=None):
        system = self.getMDSystem()
        pcaType = self.displayMDTrajPCAAnalysis.get()
        pcaFlags = ['--pca-coord', '--pca-dist']
        args = f'-i {system.getFileName()} -t {system.getTrajectoryFile()} -o {system.getSystemName()} ' \
               f' {pcaFlags[pcaType]} '
        Plugin.runScript(self, 'mdtraj_PCA.py', args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifFp(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        args = f'"{fpPkl}" --mode barcode'
        Plugin.runScript(self, self._prolifViewScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifNetwork(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        args = f'"{fpPkl}" --mode network'
        Plugin.runScript(self, self._prolifViewScript, args, env=MDTRAJ_DIC, popen=True, wait=False)

    def _showProlifMatrix(self, paramName=None):
        fpPkl = self.getMDSystem().getProlifFile()
        args = f'"{fpPkl}" --mode matrix'
        Plugin.runScript(self, self._prolifViewScript, args, env=MDTRAJ_DIC, popen=True, wait=False)


