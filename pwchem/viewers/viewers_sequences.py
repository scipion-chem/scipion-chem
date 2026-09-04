# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
import json
from subprocess import Popen
from tkinter.messagebox import askokcancel
import numpy as np
import matplotlib.pyplot as plt

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params, Protocol

from pwem.objects import SetOfSequences, Sequence, EMSet
from pwem.protocols import ProtSubSet
from pwem.viewers import EmPlotter
from pwem.viewers.mdviewer.viewer import MDViewer

from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SequenceVariants, SetOfSequenceROIs, SetOfSequencesChem, SequenceChem, MultiEpitope
from pwchem.constants import *
from pwchem.viewers.viewers_data import BioinformaticsDataViewer, plotLocalizationHistogram, plotLocalizationHistogramFromDataFrame
from pwchem.viewers.viewer_interaction import BaseInteractionViewer



seqTargets = [SequenceChem, Sequence, SetOfSequences,
              SequenceVariants, SetOfSequenceROIs, MultiEpitope]



class SequenceAliView(pwviewer.CommandView):
    """ View for calling an external command. """

    def __init__(self, seqFiles, cwd, **kwargs):
      pwviewer.CommandView.__init__(self, f'{pwchem_plugin.getProgramHome(ALIVIEW_DIC)}/aliview/aliview '
                                          f'{" ".join(seqFiles)}',
                                    cwd=cwd, **kwargs)

    def show(self):
        Popen(self._cmd, cwd=self._cwd, shell=True)

class SequenceAliViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of sequence objects
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [SetOfSequencesChem, SequenceChem, Sequence, SetOfSequences,
                SequenceVariants]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def getProtocol(self):
        if hasattr(self, 'protocol') and isinstance(self.protocol, Protocol):
            return self.protocol

    def showDefView(self, obj, outDir):
        outPath = os.path.join(outDir, 'viewSequences.fasta')
        if os.path.exists(outPath):
            os.remove(outPath)
        obj.exportToFile(outPath)
        return outPath

    def getOutDir(self):
        return os.path.abspath(self.getProtocol()._getExtraPath()
                               if self.getProtocol() else self.getProject().getTmpPath())

    def _visualize(self, obj, **kwargs):
        outDir = self.getOutDir()
        views, seqFiles = [], []
        if isinstance(obj, SetOfSequences) or isinstance(obj, Sequence):
            if hasattr(obj, '_aligned') and getattr(obj, '_aligned'):
                seqFiles += [os.path.abspath(obj.getAlignmentFileName())]
            else:
                seqFiles += [self.showDefView(obj, outDir)]

        elif isinstance(obj, SequenceVariants):
            seqFiles += [self.showDefView(obj, outDir)]

        elif isinstance(obj, SetOfSequenceROIs) or isinstance(obj, MultiEpitope):
            outPath = os.path.join(outDir, f'viewSequences_{obj.getSequenceObj().getId()}.fasta')
            obj.exportToFile(outPath)

            seqFiles += [outPath]
        views.append(SequenceAliView(seqFiles, cwd=outDir))

        return views

class SequenceGeneralViewer(BaseInteractionViewer):
  """ Protocol viewer to visualize different type of sequence objects
  """
  _label = 'Sequence viewer'
  _targets = seqTargets + [SetOfSequencesChem]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
      super().__init__(**kwargs)

  def _defineParams(self, form):
    form.addSection(label='Sequence viewer')
    aGroup = form.addGroup('AliView viewer')
    aGroup.addParam('aliLabel', params.LabelParam, label='Display sequences with AliView: ',
                    help='Display the output sequences using AliView')

    if isinstance(self.getOutSequences(), EMSet):
      tGroup = form.addGroup('Table viewer')
      tGroup.addParam('tableLabel', params.LabelParam, label='Display sequences in table format: ',
                      help='Display the output sequences using table format')

    if (isinstance(self.getOutSequences(), SetOfSequencesChem)
            and self.checkIfInteractions()):
        BaseInteractionViewer._defineParams(self, form)

    if self._hasDeepLocResiduePredictions():
        form.addSection(label='DeepLoc residue importance')
        form.addParam(
            'viewLocalization',
            params.LabelParam,
            label='Display localization probabilities: ',
            help='Display the DeepLoc predicted localization probabilities.'
        )
        form.addParam(
            'viewResidueImportance',
            params.LabelParam,
            label='Display residue importance: ',
            help='Display the DeepLoc residue importance over the sequence.'
        )

  def _hasDeepLocResiduePredictions(self):
      seqSet = self.getOutSequences()

      if isinstance(seqSet, SequenceChem):
          sequences = [seqSet]
      elif isinstance(seqSet, SetOfSequencesChem):
          sequences = list(seqSet)
      else:
          return False

      for sequence in sequences:
          try:
              attrDic = sequence.getAttributesDic()
          except (FileNotFoundError, AttributeError):
              continue

          if any(
                  attrName.startswith('DeepLoc_')
                  for attrName in attrDic
          ):
              return True

      return False

  def _getVisualizeDict(self):
      visDict = {
          'aliLabel': self._viewSeqSet,
          'tableLabel': self._viewTable,
      }

      seqSet = self.getOutSequences()

      if isinstance(seqSet, SetOfSequencesChem) and self.checkIfInteractions():
          visDict.update(BaseInteractionViewer._getVisualizeDict(self))

      if self._hasDeepLocResiduePredictions():
          visDict['viewLocalization'] = self._showLocalization
          visDict['viewResidueImportance'] = self._showResidueImportance

      return visDict

  def getOutSequences(self):
      if hasattr(self.protocol, 'outputSequences'):
          return self.protocol.outputSequences

      if isinstance(self.protocol, SetOfSequencesChem):
          return self.protocol

      if hasattr(self.protocol, 'iterOutputAttributes'):
          for oAttr in self.protocol.iterOutputAttributes():
              obj = getattr(self.protocol, oAttr[0])
              if isinstance(obj, tuple(self._targets)):
                  return obj

      return self.protocol

  def _viewSeqSet(self, e=None):
    seqSet = self.getOutSequences()
    setV = SequenceAliViewer(project=self.getProject())
    views = setV._visualize(seqSet)
    return views

  def _showLocalization(self, e=None):
      seqSet = self.getOutSequences()

      localizationPerc = getattr(seqSet, '_localizationPerc', None)

      if localizationPerc is None:
          return

      if hasattr(localizationPerc, 'get'):
          localizationPerc = localizationPerc.get()

      if not os.path.exists(localizationPerc):
          raise FileNotFoundError(
              f"Localization probability file not found: {localizationPerc}"
          )

      plotLocalizationHistogram(localizationPerc)

  def _showResidueImportance(self, paramName=None):
      seqSet = self.getOutSequences()

      if isinstance(seqSet, SequenceChem):
          sequences = [seqSet]
      else:
          sequences = seqSet.iterItems()

      sequenceData = []

      for sequence in sequences:
          attrDic = sequence.getAttributesDic()

          for attrName, values in attrDic.items():
              if attrName.startswith('DeepLoc_'):
                  sequenceData.append({
                      'name': sequence.getSeqName(),
                      'attrName': attrName,
                      'values': list(map(float, values))
                  })

      if sequenceData:
          SequenceAttributeViewer(sequenceData)

          plt.show(block=False)

  def _viewTable(self, e=None):
    seqSet = self.getOutSequences()
    try:
      setV = MDViewer(project=self.getProject())
    except:
      setV = BioinformaticsDataViewer(project=self.getProject())
    views = setV._visualize(seqSet)
    return views

  def checkIfProtocol(self):
    if isinstance(self.protocol, Protocol):
      return True
    else:
      return False

  def getOutputSet(self):
    return self.protocol

  def checkIfInteractions(self):
      seqSet = self.getOutSequences()

      if isinstance(seqSet, SetOfSequencesChem):
          return seqSet.getInteractMols() is not None

      return False

  def getInteractionSet(self):
      return self.getOutSequences()

  def _getEntityNames(self, data):
      seqNames = sorted(data.keys())

      molNames = set()
      scoreTypes = set()

      for seq in data:
          for mol in data[seq]:
              molNames.add(mol)

              for sc in data[seq][mol]:
                  scoreTypes.add(sc)

      return seqNames, sorted(molNames), sorted(scoreTypes)

  def _getLabels(self):
      return "Sequence", "Molecule", "Interaction score"

  def _getMolSet(self):
      return self.getOutSequences().getInteractMols()

  def _generateProts(self, paramName=None):
      project = self.getProject()
      structs = self.getInteractionSet()
      if hasattr(structs, '_getData'):
          data = structs._getData()

          f1 = self.getEnumText('chooseEnt1')
          f2 = self.getEnumText('chooseEnt2')
          fScore = self.getEnumText('chooseScore')

          _, seqNames, _, _ = self._getFilteredData(
              data, f1, f2, fScore
          )

          objIds = []

          for seq in self.getOutSequences():
              if seq.getSeqName() in seqNames:
                  objIds.append(str(seq.getObjId()))

          if not objIds:
              return

          if askokcancel(
                  "Generate sequences subset",
                  f"Generate subset with {len(objIds)} proteins?"):
              prot = project.newProtocol(
                  ProtSubSet,
                  inputFullSet=self.getOutSequences(),
                  selectIds=True,
                  range=','.join(objIds)
              )

              prot.setObjLabel('Filtered sequences')

              project.launchProtocol(prot, wait=True)

from matplotlib.widgets import Button

class SequenceAttributeViewer:
    def __init__(self, sequenceData):
        self.sequenceData = sequenceData
        self.index = 0

        self.fig = plt.figure(figsize=(10, 6))
        self.ax = self.fig.add_axes([0.08, 0.20, 0.88, 0.70])

        axPrev = self.fig.add_axes([0.25, 0.05, 0.18, 0.07])
        axNext = self.fig.add_axes([0.57, 0.05, 0.18, 0.07])

        self.prevButton = Button(axPrev, '◀ Previous')
        self.nextButton = Button(axNext, 'Next ▶')

        self.prevButton.on_clicked(self.previous)
        self.nextButton.on_clicked(self.next)

        self.update()

    def update(self):
        data = self.sequenceData[self.index]

        values = data['values']
        seqName = data['name']
        attrName = data['attrName']

        self.ax.clear()

        xs = np.arange(len(values))
        self.ax.bar(xs, values)

        self.ax.set_xlabel("Sequence position")
        self.ax.set_ylabel("{} value".format(attrName))

        self.ax.set_title(
            "{} — Sequence {}/{} — {}".format(
                attrName,
                self.index + 1,
                len(self.sequenceData),
                seqName
            )
        )

        self.ax.set_xlim(-1, len(values))

        maxY = max(values)
        self.ax.set_ylim(0, maxY + maxY / 10)

        self.fig.canvas.draw_idle()

    def next(self, event):
        if self.index < len(self.sequenceData) - 1:
            self.index += 1
            self.update()

    def previous(self, event):
        if self.index > 0:
            self.index -= 1
            self.update()