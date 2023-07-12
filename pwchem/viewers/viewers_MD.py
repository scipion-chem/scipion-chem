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

from pwchem.viewers import PyMolViewer, PyMolView, VmdViewPopen
from pwchem.objects import MDSystem
from pwchem.constants import *


from pwchem.viewers import PyMolViewer, PyMolView

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
        outPml = os.path.join(os.path.dirname(trjFile), 'pymolSimulation.pml')
        with open(outPml, 'w') as f:
          f.write(PML_MD_STR.format(os.path.abspath(systemFile),
                                    os.path.abspath(trjFile)))

        return [PyMolView(os.path.abspath(outPml), cwd=os.path.dirname(trjFile))]

class MDSystemPViewer(pwviewer.ProtocolViewer):
    """ Visualize the output of Molecular Dynamics simulation """
    _label = 'Viewer Molecular Dynamics System'
    _targets = [MDSystem]

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

    def _defineParams(self, form):
      form.addSection(label='Visualization of MD System')
      group = form.addGroup('Open MD system')
      group.addParam('displayPymol', params.LabelParam,
                     label='Open system in PyMol: ',
                     help='Display System in Pymol GUI.')

      if self.getMDSystem().hasTrajectory():
          self._defineSimParams(form)

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
      }

    def _showPymol(self, paramName=None):
      system = self.getMDSystem()
      return MDSystemViewer(project=self.getProject())._visualize(system, onlySystem=True)

    def _showMdPymol(self, paramName=None):
      system = self.getMDSystem()
      return MDSystemViewer(project=self.getProject())._visualize(system)

    def _showMdVMD(self, paramName=None):
      system = self.getMDSystem()

      outTcl = os.path.join(os.path.dirname(system.getTrajectoryFile()), 'vmdSimulation.tcl')
      systExt = os.path.splitext(system.getSystemFile())[1][1:]
      trjExt = os.path.splitext(system.getTrajectoryFile())[1][1:]
      with open(outTcl, 'w') as f:
        f.write(TCL_MD_STR % (system.getSystemFile(), systExt, system.getTrajectoryFile(), trjExt))
      args = '-e {}'.format(outTcl)

      return [VmdViewPopen(args)]