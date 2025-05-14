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
from subprocess import Popen

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol.params import EnumParam
import pyworkflow.utils as pwutils

import pwem.viewers.views as pwemViews
from pwem.viewers import Chimera, ChimeraView
import pwem.viewers.showj as showj
from pwem.protocols import EMProtocol
from pwem.objects import SetOfSequences, AtomStruct, SetOfAtomStructs

import pwchem.objects
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import *

class PyMol:
  """ Help class to run PyMol and manage its environment. """

  @classmethod
  def getEnviron(cls):
    """ Return the proper environ to launch PyMol.
    PyMol_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = pwutils.Environ(os.environ)
    environ.set('PATH', pwchemPlugin.getProgramHome(OPENBABEL_DIC, path='bin'),
                position=pwutils.Environ.BEGIN)
    return environ


class PyMolView(pwviewer.CommandView):
  """ View for calling an external command. """

  def __init__(self, pymolArgs, cwd, **kwargs):
    print('command: ', [self.getPymolBin(), *pymolArgs.split()])
    pwviewer.CommandView.__init__(self, [self.getPymolBin(), *pymolArgs.split()],
                                  cwd=cwd, **kwargs)
    
  def getPymolBin(self):
    return pwchemPlugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')

  def show(self):
    Popen(self._cmd, cwd=self._cwd, env=PyMol.getEnviron())


class PyMolViewer(pwviewer.Viewer):
  """ Wrapper to visualize pml objects with PyMol viewer. """
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **args):
    pwviewer.Viewer.__init__(self, **args)

  def _visualize(self, pymolFile, cwd=None, **args):
    view = PyMolView(pymolFile, cwd)
    return [view]

class VmdViewPopen(pwviewer.CommandView):
  def __init__(self, vmdArgs, **kwargs):
    pwviewer.CommandView.__init__(self, 'vmd ' + vmdArgs, **kwargs)

  def show(self):
      fullProgram = '%s && %s' % (pwchemPlugin.getEnvActivationCommand(VMD_DIC), self._cmd)
      Popen(fullProgram, cwd=self._cwd, shell=True)

class AtomStructViewer(pwviewer.ProtocolViewer):
    _label = 'Viewer AtomStruct'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [AtomStruct]
    _viewerOptions = ['PyMol', 'ChimeraX']

    def _defineParams(self, form):
      form.addSection(label='Visualization of AtomStruct')
      group = form.addGroup('AtomStruct General Viewer')
      group.addParam('displaySoftware', EnumParam,
                     choices=self._viewerOptions, default=0,
                     label='Display AtomStruct with: ',
                     help='Display the AtomStruct object with which software.\nAvailable: PyMol, ChimeraX')

    def _getVisualizeDict(self):
      return {
        'displaySoftware': self._viewAtomStruct,
      }

    def _viewAtomStruct(self, e=None):
      if self.displaySoftware.get() == 0:
        pymolViewer = AtomStructPymolViewer(project=self.getProject())
        return pymolViewer._visualize(self.getAtomStruct())
      elif self.displaySoftware.get() == 1:
        return self._viewChimera(self.getAtomStruct())

    def getAtomStruct(self):
      obj = self.protocol
      # If the input is a protocol (Analyze results was used), extract the AtomStruct obj
      if issubclass(type(obj), EMProtocol):
        for output in self.protocol.iterOutputAttributes(outputClass=AtomStruct):
          obj = output[1]
      return obj

    def _viewChimera(self, obj):
      fnCmd = os.path.abspath(self._getPath("chimera_output.cxc"))
      with open(fnCmd, 'w') as f:
        f.write('cd %s\n' % os.getcwd())
        f.write("cofr 0,0,0\n")  # set center of coordinates
        f.write("style stick\n")

        _inputVol = obj.getVolume()
        if _inputVol is not None:
          volID = 1
          dim, sampling = _inputVol.getDim()[0], _inputVol.getSamplingRate()

          f.write("open %s\n" % _inputVol.getFileName())
          x, y, z = _inputVol.getOrigin(force=True).getShifts()
          f.write("volume #%d style surface voxelSize %f\nvolume #%d origin "
                  "%0.2f,%0.2f,%0.2f\n"
                  % (volID, sampling, volID, x, y, z))
        else:
          dim, sampling = 150., 1.

        bildFileName = self._getPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim, bildFileName=bildFileName, sampling=sampling)
        f.write("open %s\n" % bildFileName)

        f.write("open %s\n" % obj.getFileName())

        view = ChimeraView(fnCmd)
        return [view]


class AtomStructPymolViewer(PyMolViewer):
    _label = 'Pymol viewer AtomStruct'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = []

    def _visualize(self, obj, **args):
      pymolV = PyMolViewer(project=self.getProject())
      return pymolV._visualize(obj.getFileName())


class SetOfDatabaseIDView(pwemViews.ObjectView):
    """ Customized ObjectView for SetOfDatabaseID. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.MODE: 'metadata',
                             showj.RENDER: '_PDBLigandImage'}
        defaultViewParams.update(viewParams)
        pwemViews.ObjectView.__init__(self, project, inputid, path, other,
                                  defaultViewParams, **kwargs)

class BioinformaticsDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        SetOfSequences,
        pwchem.objects.SetOfDatabaseID,
        pwchem.objects.SetOfSmallMolecules,
        pwchem.objects.SetOfBindingSites,
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return pwemViews.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        # For now handle both types of SetOfTiltSeries together
        if issubclass(cls, pwchem.objects.SetOfDatabaseID):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, SetOfSequences):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfSmallMolecules):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfBindingSites):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfStructROIs):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfSequenceROIs):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))

        return views
