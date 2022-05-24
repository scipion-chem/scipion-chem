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
from subprocess import Popen

import pyworkflow.viewer as pwviewer
from pwem.objects import SetOfSequences, Sequence

from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SequencesAlignment, SequenceVariants, SetOfSequenceROIs
from pwchem.constants import *


class SequenceAliView(pwviewer.CommandView):
    """ View for calling an external command. """

    def __init__(self, seqFiles, cwd, **kwargs):
        pwviewer.CommandView.__init__(self, '{}/aliview/aliview {}'.
                                      format(pwchem_plugin.getProgramHome(ALIVIEW_DIC), ' '.join(seqFiles)),
                                      cwd=cwd, **kwargs)

    def show(self):
        Popen(self._cmd, cwd=self._cwd, shell=True)


class SequenceAliViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of sequence objects
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [Sequence, SetOfSequences, SequencesAlignment, SequenceVariants, SetOfSequenceROIs]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        seqFiles = []
        if issubclass(cls, SetOfSequences) or issubclass(cls, Sequence):
            outPath = os.path.abspath(self.protocol._getExtraPath('viewSequences.fasta'))
            seqFiles += [outPath]
            if os.path.exists(outPath):
                os.remove(outPath)
            obj.exportToFile(outPath)

        elif issubclass(cls, SequencesAlignment):
            seqFiles += [obj.getAlignmentFileName()]

        elif issubclass(cls, SequenceVariants):
            outPath = os.path.abspath(self.protocol._getExtraPath('viewSequences.fasta'))
            seqFiles += [outPath]
            if os.path.exists(outPath):
                os.remove(outPath)
            obj._sequence.exportToFile(outPath)

        elif issubclass(cls, SetOfSequenceROIs):
            outPath = os.path.abspath(self.protocol._getExtraPath('viewSequences_{}.fasta'.
                                                                  format(obj.getSequenceObj().getId())))
            seqFiles += [outPath]
            if os.path.exists(outPath):
                os.remove(outPath)
            wholeSeq = obj.getFirstItem().getSequence()
            obj.getFirstItem()._sequence.exportToFile(outPath)
            for roi in obj:
                roiSeq, roiIdx = roi.getROISequence(), roi.getROIIdx()
                tmpSeq = (roiIdx-1) * '-' + roiSeq + (len(wholeSeq) - len(roiSeq) - roiIdx + 1) * '-'
                tmpSeqObj = Sequence(sequence=tmpSeq, id=roi._ROISequence.getId())
                tmpSeqObj.appendToFile(outPath)

        views.append(SequenceAliView(seqFiles, cwd=self.protocol._getExtraPath()))

        return views

