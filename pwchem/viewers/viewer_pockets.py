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

import pyworkflow.protocol.params as params
import pyworkflow.viewer as pwviewer

from pwem.viewers.mdviewer.viewer import MDViewer

from pwchem.objects import SetOfStructROIs
from pwchem.viewers.viewers_data import BioinformaticsDataViewer, PyMolViewer, VmdViewPopen
from pwchem.constants import *
from pwchem.protocols import ProtocolConsensusStructROIs


class StructROIPointsViewer(pwviewer.Viewer):
  _label = 'Viewer structural ROIs'
  _environments = [pwviewer.DESKTOP_TKINTER]

  # _targets = [SetOfPockets]

  def _visualize(self, obj, bBox=False, **kwargs):
    hetatmFile = obj.buildPDBhetatmFile()
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


VOLUME_PYMOL, VOLUME_PYMOL_SURF = 0, 1


class ViewerGeneralStructROIs(pwviewer.ProtocolViewer):
  _label = 'Viewer structural ROIs'
  _targets = [SetOfStructROIs]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of structural ROIs')
    group = form.addGroup('Pymol General Viewer')
    group.addParam('displayAtomStruct', params.EnumParam,
                  choices=['PyMol (ROI Points)', 'PyMol (Contact Surface)'],
                  default=VOLUME_PYMOL,
                  display=params.EnumParam.DISPLAY_HLIST,
                  label='Display output AtomStruct with',
                  help='*PyMol*: display AtomStruct and structural ROIs as points / surface.'
                  )
    group.addParam('displayBBoxes', params.BooleanParam,
                  default=False, label='Display ROIs bounding boxes',
                  help='Display the bounding boxes in pymol to check the size for the localized docking')
    group.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                  default=1.1, condition='displayBBoxes',
                  help='The radius * n of each ROI will be used as grid radius')

    form.addSection(label='Table view')
    form.addParam('displayTable', params.LabelParam,
                  label='Display ROIs set and attributes in table format: ',
                  help='Display the ROIs set in the set in table format with their respective attributes')

  def _getVisualizeDict(self):
    return {
      'displayAtomStruct': self._showAtomStruct,
      'displayTable': self._viewSet,
    }

  def _viewSet(self, e=None):
    if type(self.protocol) == SetOfStructROIs:
      molSet = self.protocol
    elif hasattr(self.protocol, 'outputStructROIs'):
      molSet = getattr(self.protocol, 'outputStructROIs')
    else:
      print('Cannot find outputStructROIs')

    try:
      setV = MDViewer(project=self.getProject())
    except:
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

    pymolV = StructROIPointsViewer(project=self.getProject())
    if type(self.protocol) == SetOfStructROIs:
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputStructROIs'):
      return pymolV._visualize(getattr(self.protocol, 'outputStructROIs'), bBox=bBox)

  def _showAtomStructPyMolSurf(self):
    bBox = self.displayBBoxes.get()
    if bBox:
      bBox = self.pocketRadiusN.get()

    pymolV = ContactSurfaceViewer(project=self.getProject())
    if type(self.protocol) == SetOfStructROIs:
      return pymolV._visualize(self.protocol, bBox=bBox)
    elif hasattr(self.protocol, 'outputStructROIs'):
      return pymolV._visualize(getattr(self.protocol, 'outputStructROIs'), bBox=bBox)

  def _showAtomStructPyMol(self, pmlFile, outDir):
    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(pmlFile, cwd=outDir)


MIXED, FPOCKET, P2RANK, AUTOLIGAND, SITEMAP = 'Mixed', 'FPocket', 'P2Rank', 'AutoLigand', 'Sitemap'
VOLUME_VMD = 2

class ViewerConsensusStructROIs(pwviewer.ProtocolViewer):
    _label = 'Viewer consensus structural ROIs'
    _targets = [ProtocolConsensusStructROIs]

    def __init__(self, **kwargs):
      pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
      form.addSection(label='Visualization of consensus Structural ROIs')
      pGroup = form.addGroup('Pymol General Viewer')
      pGroup.addParam('outputSet', params.EnumParam,
                    choices=self.getChoices(), default=0,
                    label='Set of pockets output: ',
                    help='Set of pockets output to visualize'
                    )
      pGroup.addParam('setClass', params.StringParam,
                    label='Set of Pockets Class: ', default='-',
                    help='Pocket class in chosen set')

      pGroup.addParam('displayFPocket', params.EnumParam,
                    choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)', 'VMD'],
                    default=0, condition='setClass=="{}"'.format(FPOCKET),
                    display=params.EnumParam.DISPLAY_HLIST,
                    label='Display output Set of Pockets with: ',
                    help='*PyMol*: display Set of Pockets and pockets as points / surface.\n '
                         '*VMD*: display Set of Pockets with VMD.'
                    )

      pGroup.addParam('displayPyMol', params.EnumParam,
                    choices=['PyMol (Pocket Points)', 'PyMol (Contact Surface)'],
                    default=0, condition='setClass!="{}"'.format(FPOCKET),
                    display=params.EnumParam.DISPLAY_HLIST,
                    label='Display output Set of Pockets with: ',
                    help='*PyMol*: display Set of Pockets and pockets as points / surface'
                    )
      pGroup.addParam('displayBBoxes', params.BooleanParam,
                    default=False, label='Display pocket bounding boxes',
                    condition='not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                    help='Display the bounding boxes in pymol to check the size for the localized docking')
      pGroup.addParam('pocketRadiusN', params.FloatParam, label='Grid radius vs pocket radius: ',
                    default=1.1, condition='displayBBoxes and not (setClass=="{}" and displayFPocket==0)'.format(FPOCKET),
                    help='The radius * n of each pocket will be used as grid radius')
      
      tGroup = form.addGroup('Table view')
      tGroup.addParam('displayTable', params.LabelParam,
                      label='Display table view of set of Pockets: ',
                      help='Table view')

    def getChoices(self):
        outputLabels = []
        for oAttr in self.protocol.iterOutputAttributes():
          outputLabels.append(oAttr[0])
        outputLabels.sort()
        return outputLabels

    def _getVisualizeDict(self):
      return {
        'displayTable': self._showTable,
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
      if self.displayPyMol == VOLUME_PYMOL:
        return self._showAtomStructPyMolPoints()

      elif self.displayPyMol == VOLUME_PYMOL_SURF:
        return self._showAtomStructPyMolSurf()

    #Display functions
    def _showTable(self, paramName=None):
      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      try:
        setV = MDViewer(project=self.getProject())
      except:
        setV = BioinformaticsDataViewer(project=self.getProject())
      return setV._visualize(outPockets)

    def _showAtomStructPyMolPoints(self):
      bBox = self.displayBBoxes.get()
      if bBox:
        bBox = self.pocketRadiusN.get()

      outPockets = getattr(self.protocol, self.getEnumText('outputSet'))
      pymolV = StructROIPointsViewer(project=self.getProject())
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
