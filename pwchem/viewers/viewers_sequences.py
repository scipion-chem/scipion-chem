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
from pwchem.constants import *
from pwchem.viewers.viewers_data import BioinformaticsDataViewer
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
  
class SequenceChemViewer(BaseInteractionViewer):
    """ Protocol viewer to visualize different type of sequence objects
    """
    _label = 'Sequence chem viewer'
    _targets = [SetOfSequencesChem]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def __init__(self, **kwargs):
        super().__init__(**kwargs)


    def _getData(self):
        outSeqs = self.getOutSequences()

        data = {}

        for seq in outSeqs:
            seqName = seq.getSeqName()

            if seqName not in data:
                data[seqName] = {}

            interactMols = seq.getInteractMols()

            if interactMols is None:
                continue

            for molName, scores in interactMols.items():
                data[seqName][molName] = scores

        return data

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
        data = self._getData()

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

    # ---------------------------------------------------

    def checkIfInteractions(self):
        if not self.checkIfProtocol():
            seqSet = self.protocol
        else:
            for oAttr in self.protocol.iterOutputAttributes():
                obj = getattr(self.protocol, oAttr[0])

                if isinstance(obj, SetOfSequencesChem):
                    seqSet = obj

        return seqSet.getInteractMols() is not None
