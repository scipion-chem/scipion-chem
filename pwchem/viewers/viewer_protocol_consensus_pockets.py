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

from ..protocols import ProtocolConsensusPockets
import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer
from ..viewers import PyMolViewer, PocketPointsViewer, ContactSurfaceViewer, VmdViewPopen

from subprocess import Popen

MIXED, FPOCKET, P2RANK, AUTOLIGAND, SITEMAP = 'Mixed', 'FPocket', 'P2Rank', 'AutoLigand', 'Sitemap'
VOLUME_VMD, VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1, 2

class ViewerConsensusPockets(pwviewer.ProtocolViewer):
  _label = 'Viewer pockets'
  _targets = [ProtocolConsensusPockets]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of consensus pockets')
    form.addParam('outputSet', params.EnumParam,
                  choices=self.getChoices(), default=0,
                  label='Set of pockets output: ',
                  help='Set of pockets output to visualize'
                  )
    form.addParam('setClass', params.StringParam,
                  label='Set of Pockets Class: ', default='-',
                  help='Pocket class in chosen set')

    form.addParam('displayFPocket', params.EnumParam,
                  choices=['VMD', 'PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=0, condition='setClass=="{}"'.format(FPOCKET),
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output Set of Pockets with: ',
                  help='*PyMol*: display Set of Pockets and pockets as points / surface.\n '
                       '*VMD*: display Set of Pockets with VMD.'
                  )

    form.addParam('displayPyMol', params.EnumParam,
                  choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                  default=0, condition='setClass!="{}"'.format(FPOCKET),
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output Set of Pockets with: ',
                  help='*PyMol*: display Set of Pockets and pockets as points / surface'
                  )
    form.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display pocket bounding boxes',
                  condition='not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    form.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes and not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                  help='The radius * n of each pocket will be used as grid radius')

  def getChoices(self):
      outputLabels = []
      for oAttr in self.protocol.iterOutputAttributes():
        outputLabels.append(oAttr[0])
      outputLabels.sort()
      return outputLabels

  def _getVisualizeDict(self):
    return {
      'displayFPocket': self._showFPocketPockets,
      'displayPyMol': self._showStandardPockets,
    }

  def _validate(self):
    return []

  # =========================================================================
  # ShowAtomStructs
  # =========================================================================

  def getOutputAtomStructFile(self):
    return os.path.abspath(self.protocol.outputAtomStruct.getFileName())

  def _showFPocketPockets(self, paramName=None):
    if self.displayFPocket == VOLUME_PYMOL:
      return self._showAtomStructPyMolPoints()

    elif self.displayFPocket == VOLUME_PYMOL_SURF:
      return self._showAtomStructPyMolSurf()

    elif self.displayFPocket == VOLUME_VMD:
      return self._showAtomStructVMD()

  def _showStandardPockets(self, paramName=None):
    if self.displayPyMol == VOLUME_PYMOL-1:
      return self._showAtomStructPyMolPoints()

    elif self.displayPyMol == VOLUME_PYMOL_SURF-1:
      return self._showAtomStructPyMolSurf()

  #Display functions
  def _showAtomStructPyMolPoints(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
    pymolV = PocketPointsViewer(project=self.getProject())
    return pymolV._visualize(outPockets, bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
    pymolV = ContactSurfaceViewer(project=self.getProject())
    return pymolV._visualize(outPockets, bBox=bBox)

  def _showAtomStructVMD(self):
    outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
    outFile = outPockets.getProteinHetatmFile().split('/')[-1]
    tclFile = outFile.replace('_out.pdb', '.tcl')
    outDir = os.path.abspath(outPockets.getSetDir())
    args = '{} -e {}'.format(outFile, tclFile)

    return [VmdViewPopen(args, cwd=outDir)]
