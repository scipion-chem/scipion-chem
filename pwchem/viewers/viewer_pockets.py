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
import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
from pwem.viewers import Vmd, VmdView
from pwchem.viewers.viewers_data import BioinformaticsDataViewer, PyMolViewer
from ..constants import *


class PocketPointsViewer(pwviewer.Viewer):
  _label = 'Viewer pocket points'
  _environments = [pwviewer.DESKTOP_TKINTER]

  # _targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    outHETMFile = obj.buildPDBhetatmFile()
    pmlFile = obj.createPML(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=os.path.dirname(pmlFile))


class ContactSurfaceViewer(pwviewer.Viewer):
  _label = 'Viewer contact surface'
  _environments = [pwviewer.DESKTOP_TKINTER]

  # _targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    pmlFile = obj.createSurfacePml(bBox=bBox)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=os.path.dirname(pmlFile))


class VmdViewPopen(VmdView):
  def __init__(self, vmdArgs, **kwargs):
    pwviewer.CommandView.__init__(self, 'vmd ' + vmdArgs,
                                  env=Vmd.getEnviron(), **kwargs)

  def show(self):
    Popen(self._cmd, cwd=self._cwd, env=Vmd.getEnviron(), shell=True)


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
    return setV._visualize(molSet)

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
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputPockets'):
      return pymolV._visualize(getattr(self.protocol, 'outputPockets'), bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = ContactSurfaceViewer(project=self.getProject())
    if type(self.protocol) == SetOfPockets:
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputPockets'):
      return pymolV._visualize(getattr(self.protocol, 'outputPockets'), bBox=bBox)

  def _showAtomStructPyMol(self, pmlFile, outDir):
    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=outDir)
