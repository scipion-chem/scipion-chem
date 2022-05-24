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
import pwem.viewers.views as pwemViews
import pwem.viewers.showj as showj
from pwem.objects import SetOfSequences, AtomStruct
import pyworkflow.utils as pwutils

import pwchem.objects
from pwchem import Plugin as pwchem_plugin
from ..constants import *

class PyMol:
  """ Help class to run PyMol and manage its environment. """

  @classmethod
  def getEnviron(cls):
    """ Return the proper environ to launch PyMol.
    PyMol_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = pwutils.Environ(os.environ)
    environ.set('PATH', pwchem_plugin.getProgramHome(PYMOL_DIC, path='bin'),
                position=pwutils.Environ.BEGIN)
    return environ


class PyMolView(pwviewer.CommandView):
  """ View for calling an external command. """

  def __init__(self, pymolArgs, cwd, **kwargs):
    print('command: ', [pwchem_plugin.getProgramHome(PYMOL_DIC), *pymolArgs.split()])
    pwviewer.CommandView.__init__(self, [pwchem_plugin.getProgramHome(PYMOL_DIC, 'bin/pymol'), *pymolArgs.split()],
                                  cwd=cwd,
                                  env=PyMol.getEnviron(), **kwargs)

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

class AtomStructPymolViewer(PyMolViewer):
    _label = 'Viewer AtomStruct'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [AtomStruct]

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
        pwchem.objects.SetOfSequenceROIs
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
        elif issubclass(cls, pwchem.objects.SetOfPockets):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfSequenceROIs):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))

        return views
