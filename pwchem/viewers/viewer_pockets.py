# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
from subprocess import Popen

from pwchem.objects import ProteinPocket, SetOfPockets
from pwchem import Plugin as pwchem_plugin
import pyworkflow.utils as pwutils
import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
from pwem.viewers import Vmd, VmdView
from pwchem.viewers import BioinformaticsDataViewer

class PyMol:
  """ Help class to run PyMol and manage its environment. """

  @classmethod
  def getEnviron(cls):
    """ Return the proper environ to launch PyMol.
    PyMol_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = pwutils.Environ(os.environ)
    environ.set('PATH', os.path.join(os.environ['PYMOL_HOME'], 'bin'),
                position=pwutils.Environ.BEGIN)
    return environ


class PyMolView(pwviewer.CommandView):
  """ View for calling an external command. """

  def __init__(self, pymolArgs, cwd, **kwargs):
    pwviewer.CommandView.__init__(self, [pwchem_plugin.getPyMolPath(), *pymolArgs.split()],
                                  cwd=cwd,
                                  env=PyMol.getEnviron(), **kwargs)

  def show(self):
    Popen(self._cmd, cwd=self._cwd, env=PyMol.getEnviron())


class PyMolViewer(pwviewer.Viewer):
  """ Wrapper to visualize pml objects with PyMol viewer. """
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **args):
    pwviewer.Viewer.__init__(self, **args)

  def visualize(self, pymolFile, cwd, **args):
    PyMolView(pymolFile, cwd).show()


class PocketPointsViewer(pwviewer.Viewer):
  _label = 'Viewer pocket points'
  _environments = [pwviewer.DESKTOP_TKINTER]
  #_targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    pmlFile = obj.createPML(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pmlFile, cwd=os.path.dirname(pmlFile))

class ContactSurfaceViewer(pwviewer.Viewer):
  _label = 'Viewer contact surface'
  _environments = [pwviewer.DESKTOP_TKINTER]
  #_targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    pmlFile = obj.createSurfacePml(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pmlFile, cwd=os.path.dirname(pmlFile))


class VmdViewFpocket(VmdView):
  def __init__(self, vmdArgs, **kwargs):
    pwviewer.CommandView.__init__(self, ['vmd', *vmdArgs.split()],
                                  env=Vmd.getEnviron(), **kwargs)

  def show(self):
    Popen(self._cmd, cwd=self._cwd, env=Vmd.getEnviron())

VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1

class ViewerGeneralPockets(pwviewer.ProtocolViewer):
  _label = 'Viewer pockets'
  _targets = [SetOfPockets]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of predicted pockets')
    form.addParam('displayAtomStruct', params.EnumParam,
                  choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=VOLUME_PYMOL,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output AtomStruct with',
                  help='*PyMol*: display AtomStruct and pockets as points / surface.'
                  )
    form.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display pocket bounding boxes',
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    form.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes',
                  help='The radius * n of each pocket will be used as grid radius')

    form.addSection(label='Table view')
    form.addParam('displayTable', params.LabelParam,
                  label='Display pockets set and attributes in table format: ',
                  help='Display the pockets set in the set in table format with their respective attributes')

  def _getVisualizeDict(self):
    return {
      'displayAtomStruct': self._showAtomStruct,
      'displayTable': self._viewSet,
    }

  def _viewSet(self, e=None):
    if type(self.protocol) == SetOfPockets:
        molSet = self.protocol
    elif hasattr(self.protocol, 'outputPockets'):
        molSet = getattr(self.protocol, 'outputPockets')
    else:
        print('Cannot find outputPockets')

    setV = BioinformaticsDataViewer(project=self.getProject())
    views = setV._visualize(molSet)
    views[0].show()

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def getOutputAtomStructFile(self):
    return os.path.abspath(self.protocol.outputAtomStruct.getFileName())

  def _showAtomStruct(self, paramName=None):
    if self.displayAtomStruct == VOLUME_PYMOL:
      return self._showAtomStructPyMolPoints()

    elif self.displayAtomStruct == VOLUME_PYMOL_SURF:
      return self._showAtomStructPyMolSurf()

  def _showAtomStructPyMolPoints(self):
    bBox = self.displayBBoxes.get()
    if bBox:
        bBox = self.pocketRadiusN.get()

    pymolV = PocketPointsViewer(project=self.getProject())
    if type(self.protocol) == SetOfPockets:
        pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputPockets'):
        pymolV._visualize(getattr(self.protocol, 'outputPockets'), bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = ContactSurfaceViewer(project=self.getProject())
    if type(self.protocol) == SetOfPockets:
        pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputPockets'):
        pymolV._visualize(getattr(self.protocol, 'outputPockets'), bBox=bBox)
    
  def _showAtomStructPyMol(self, pmlFile, outDir):
    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(pmlFile, cwd=outDir)
