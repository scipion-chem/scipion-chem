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
import numpy as np
from subprocess import Popen
import matplotlib.pyplot as plt
import matplotlib
from tkinter.messagebox import askokcancel

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params, Protocol

from pwem.objects import SetOfSequences, Sequence, EMSet
from pwem.protocols import ProtSubSet
from pwem.viewers import EmPlotter
from pwem.viewers.mdviewer.viewer import MDViewer

from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SequenceVariants, SetOfSequenceROIs, SetOfSequencesChem, SequenceChem, MultiEpitope
from pwchem.protocols import ProtExtractInteractingMols
from pwchem.constants import *
from pwchem.utils import getFilteredOutput
from pwchem.viewers import BioinformaticsDataViewer

seqTargets = [SequenceChem, Sequence, SetOfSequences,
              SequenceVariants, SetOfSequenceROIs, MultiEpitope]

def heatmap(data, rowLabels, colLabels, ax=None, cbarKw=None, cbarLabel="", **kwargs):
  """
  Create a heatmap from a numpy array and two lists of labels.

  Parameters
  ----------
  data
      A 2D numpy array of shape (M, N).
  rowLabels
      A list or array of length M with the labels for the rows.
  colLabels
      A list or array of length N with the labels for the columns.
  ax
      A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
      not provided, use current axes or create a new one.  Optional.
  cbarKw
      A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
  cbarLabel
      The label for the colorbar.  Optional.
  **kwargs
      All other arguments are forwarded to `imshow`.
  """
  if ax is None:
    ax = plt.gca()

  if cbarKw is None:
    cbarKw = {}

  # Plot the heatmap
  im = ax.imshow(data, **kwargs)

  # Create colorbar
  cbar = ax.figure.colorbar(im, ax=ax, **cbarKw)
  cbar.ax.set_ylabel(cbarLabel, rotation=-90, va="bottom")

  # Show all ticks and label them with the respective list entries.
  ax.set_xticks(np.arange(data.shape[1]), labels=colLabels)
  ax.set_yticks(np.arange(data.shape[0]), labels=rowLabels)

  # Let the horizontal axes labeling appear on top.
  ax.tick_params(top=True, bottom=False,
                 labeltop=True, labelbottom=False)

  # Rotate the tick labels and set their alignment.
  plt.setp(ax.get_xticklabels(), rotation=-30, ha="right",
           rotation_mode="anchor")

  # Turn spines off and create white grid.
  ax.spines[:].set_visible(False)

  ax.set_xticks(np.arange(data.shape[1] + 1) - .5, minor=True)
  ax.set_yticks(np.arange(data.shape[0] + 1) - .5, minor=True)
  ax.grid(which="minor", color="w", linestyle='-', linewidth=3)
  ax.tick_params(which="minor", bottom=False, left=False)

  return im, cbar

def annotateHeatmap(im, data=None, valfmt="{x:.2f}",
                     textcolors=("black", "white"),
                     threshold=None, **textkw):
  """
  A function to annotate a heatmap.

  Parameters
  ----------
  im
      The AxesImage to be labeled.
  data
      Data used to annotate.  If None, the image's data is used.  Optional.
  valfmt
      The format of the annotations inside the heatmap.  This should either
      use the string format method, e.g. "$ {x:.2f}", or be a
      `matplotlib.ticker.Formatter`.  Optional.
  textcolors
      A pair of colors.  The first is used for values below a threshold,
      the second for those above.  Optional.
  threshold
      Value in data units according to which the colors from textcolors are
      applied.  If None (the default) uses the middle of the colormap as
      separation.  Optional.
  **kwargs
      All other arguments are forwarded to each call to `text` used to create
      the text labels.
  """
  if not isinstance(data, (list, np.ndarray)):
    data = im.get_array()

  # Normalize the threshold to the images color range.
  if threshold is not None:
    threshold = im.norm(threshold)
  else:
    threshold = im.norm(data.max()) / 2.

  # Set default alignment to center, but allow it to be
  # overwritten by textkw.
  kw = dict(horizontalalignment="center",
            verticalalignment="center")
  kw.update(textkw)

  # Get the formatter in case a string is supplied
  if isinstance(valfmt, str):
    valfmt = matplotlib.ticker.StrMethodFormatter(valfmt)

  # Loop over the data and create a `Text` for each "pixel".
  # Change the text's color depending on the data.
  texts = []
  for i in range(data.shape[0]):
    for j in range(data.shape[1]):
      kw.update(color=textcolors[int(im.norm(data[i, j]) > threshold)])
      text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
      texts.append(text)

  return texts

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

class SequenceGeneralViewer(pwviewer.ProtocolViewer):
  """ Protocol viewer to visualize different type of sequence objects
  """
  _label = 'Sequence viewer'
  _targets = seqTargets
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Sequence viewer')
    aGroup = form.addGroup('AliView viewer')
    aGroup.addParam('aliLabel', params.LabelParam, label='Display sequences with AliView: ',
                    help='Display the output sequences using AliView')

    if isinstance(self.getOutSequences(), EMSet):
      tGroup = form.addGroup('Table viewer')
      tGroup.addParam('tableLabel', params.LabelParam, label='Display sequences in table format: ',
                      help='Display the output sequences using table format')

  def _getVisualizeDict(self):
    return {
      # AliView
      'aliLabel': self._viewSeqSet,

      # Table
      'tableLabel': self._viewTable,
    }

  def getOutSequences(self):
    if self.checkIfProtocol():
      for oAttr in self.protocol.iterOutputAttributes():
        for oType in seqTargets:
          if isinstance(getattr(self.protocol, oAttr[0]), oType):
            return getattr(self.protocol, oAttr[0])
    else:
      return self.protocol

  def _viewSeqSet(self, e=None):
    seqSet = self.getOutSequences()
    setV = SequenceAliViewer(project=self.getProject())
    views = setV._visualize(seqSet)
    return views

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
  
class SequenceChemViewer(SequenceGeneralViewer):
    """ Protocol viewer to visualize different type of sequence objects
    """
    _label = 'Sequence chem viewer'
    _targets = [SetOfSequencesChem]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def __init__(self, **kwargs):
        pwviewer.ProtocolViewer.__init__(self, **kwargs)

    def _defineParams(self, form):
        super()._defineParams(form)
        if self.checkIfInteractions():
            self._defineInteractionParams(form)

    def _defineInteractionParams(self, form):
        outSeqs = self.getOutSequences()
        seqNames = outSeqs.getSequenceNames()
        molNames = outSeqs.getInteractMolNames() if outSeqs.getInteractMolsNumber() < 1000 else []

        form.addSection(label='Interaction viewer')
        sGroup = form.addGroup('Score filter')
        sGroup.addParam('chooseSeq', params.EnumParam, label='Display results for protein: ',
                        choices=['All'] + seqNames, default=0,
                        help='Display the selected results only for the specific protein')
        sGroup.addParam('chooseMol', params.EnumParam, label='Display results for molecules: ',
                        choices=['All'] + molNames, default=0,
                        help='Display the selected results only for the specific molecule')
        sGroup.addParam('scThres', params.FloatParam, label='Score threshold: ', default=0,
                        help='Display the interaction results over the selected score threshold')

        hGroup = form.addGroup('Scores View')
        hGroup.addParam('displayHeatMap', params.LabelParam, label='Display DTI heatmap: ',
                        help='Display a heatmap showing the scores for each protein-molecule pair')
        hGroup.addParam('displayHistogram', params.LabelParam, label='Display scores histogram: ',
                        help='Display a histogram showing the score distribution for the selected proteins '
                             'and molecules')
        hGroup.addParam('intervals', params.IntParam, default=11, label="Intervals",
                        help="Number of labels of the scale", expertLevel=params.LEVEL_ADVANCED)

        oGroup = form.addGroup('Generate output')
        oGroup.addParam('genProts', params.LabelParam, label='Generate protein sequences: ',
                        help='Generate an output of the filtered sequences of proteins'
                             '\ne.g: proteins passing the score filter specified above')
        oGroup.addParam('genMols', params.LabelParam, label='Generate small molecules: ',
                        help='Generate an output of the filtered molecules. '
                             '\ne.g: molecules passing the score filter specified above')
        return form

    def _getVisualizeDict(self):
        vDic = super()._getVisualizeDict()
        vDic.update({
            # Interaction views
            'displayHeatMap': self._viewHeatMap,
            'displayHistogram': self._viewHistogram,
            'genProts': self._generateProts,
            'genMols': self._generateMols,
        })
        return vDic

    def _viewHeatMap(self, paramName=None):
        outSeqs = self.getOutSequences()
        seqNames, nMols = outSeqs.getSequenceNames(), outSeqs.getInteractMolsNumber()
        if len(seqNames) * nMols > 1000:
            try:
                msg = "A heatmap with the scores of the selected sequence-molecule pairs will be displayed.\n" \
                      "Big sets of sequences/molecules might take a while to process and display.\n" \
                      "Do you want to continue?"
                answer = askokcancel("Display heatmap", msg, parent=None)
            except:
                answer = True
        else:
            answer = True

        if answer:
            filtSeqNames, filtMolNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
            intAr, seqNames, molNames = getFilteredOutput(outSeqs, filtSeqNames, filtMolNames, self.scThres.get())

            fig, ax = plt.subplots()
            im, _ = heatmap(intAr, seqNames, molNames, ax=ax,
                               cmap="YlGn", cbarLabel="ConPLex interaction score")
            annotateHeatmap(im, valfmt="{x:.2f}")
            fig.tight_layout()
            plt.show()

    def _viewHistogram(self, paramName=None):
      outSeqs = self.getOutSequences()
      filtSeqNames, filtMolNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
      if filtMolNames == 'All':
        intAr, seqNames, molNames = getFilteredOutput(outSeqs, filtSeqNames, filtMolNames, self.scThres.get())
        scoreList = intAr.flatten()

        self.plotter = EmPlotter(x=1, y=1, windowTitle='Score distribution')
        outName = f'{filtSeqNames}'
        a = self.plotter.createSubPlot(f"{outName} distribution", f'{outName} scores', f"Score counts")
        low, high = min(scoreList), max(scoreList)
        a.set_xlim([low, high])

        n = 5
        mult = 10 ** n
        stepSize = int(round((high - low) / self.intervals.get(), n) * mult)
        bins = [i / mult for i in range(int(low * mult), int(high * mult), stepSize)]
        _, _, bars = a.hist(scoreList, bins=bins, linewidth=1, label="Map", rwidth=0.9)

        a.grid(True)
        return [self.plotter]

    def _generateProts(self, paramName=None):
        project = self.getProject()
        outSeqs = self.getOutSequences()
        filtSeqNames, filtMolNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
        _, seqNames, _ = getFilteredOutput(outSeqs, filtSeqNames, filtMolNames, self.scThres.get())

        objIds = []
        for seq in self.getOutSequences():
            if seq.getSeqName() in seqNames:
                objIds.append(str(seq.getObjId()))

        try:
            msg = f"An output of the following {len(seqNames)} sequences will be generated:\n{','.join(seqNames)}"
            answer = askokcancel("Generate sequences output", msg, parent=None)
        except:
            answer = True

        if answer:
            protFilter = project.newProtocol(
                ProtSubSet, inputFullSet=self.getOutSequences(),
                selectIds=True, range=','.join(objIds)
            )

            protFilter.setObjLabel('Filtered sequences')
            project.launchProtocol(protFilter, wait=True)

    def _generateMols(self, paramName=None):
        project = self.getProject()
        outSeqs = self.getOutSequences()
        filtSeqNames, filtMolNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
        _, _, filtMolNames = getFilteredOutput(outSeqs, filtSeqNames, filtMolNames, self.scThres.get())

        try:
            msg = f"An output of the following {len(filtMolNames)} molecules will be generated:\n{','.join(filtMolNames)}"
            answer = askokcancel("Generate sequences output", msg, parent=None)
        except:
            answer = True

        if answer:
            seqNames, molNames = self.getEnumText('chooseSeq'), self.getEnumText('chooseMol')
            scThres = self.scThres.get()
            protFilter = project.newProtocol(
                ProtExtractInteractingMols, inputSequences=self.getOutSequences(),
                chooseSeq=seqNames, chooseMol=molNames,
                scThres=scThres
            )
            project.launchProtocol(protFilter, wait=True)

    ############ UTILS ##############

    def checkIfInteractions(self):
        if not self.checkIfProtocol():
            seqSet = self.protocol
        else:
            for oAttr in self.protocol.iterOutputAttributes():
                if type(getattr(self.protocol, oAttr[0])) == SetOfSequencesChem:
                    seqSet = getattr(self.protocol, oAttr[0])
        return seqSet.getInteractMols() is not None
