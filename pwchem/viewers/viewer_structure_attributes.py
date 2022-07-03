# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
from pwem.viewers import ChimeraAttributeViewer
from pwem.objects import Sequence
from pwchem.protocols import ProtExtractSeqsROI, ProtCalculateSASA
from pwchem.viewers.viewers_sequences import SequenceAliView

class ConservationViewer(ChimeraAttributeViewer):
    """ Viewer for attribute conservation of an AtomStruct.
      Includes visualization in chimera and in histograms"""
    _targets = [ProtExtractSeqsROI]
    _label = 'Sequence conservation analysis viewer'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization of conservation regions')
        form.addParam('viewROIs', params.LabelParam, label='View sequence conservation ROIs: ',
                       help='View sequence conservation ROIs extracted in the protocol')
        if getattr(self.protocol, 'inputAS'):
            super()._defineParams(form)
            # Overwrite defaults
            from pwem.wizards.wizard import ColorScaleWizardBase
            group = form.addGroup('Color settings')
            ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=1, defaultIntervals=21,
                                                        defaultColorMap='RdBu_r')

    def _getVisualizeDict(self):
        visDic = {'viewROIs': self._showROIs}
        if getattr(self.protocol, 'inputAS'):
            visDic.update(super()._getVisualizeDict())
        return visDic

    def _showROIs(self, paramName=None):
        obj = self.protocol.outputROIs
        outPath = os.path.abspath(self.protocol._getExtraPath('viewSequences_{}.fasta'.
                                                              format(obj.getSequenceObj().getId())))
        obj.exportToFile(outPath)
        return [SequenceAliView([outPath], cwd=self.protocol._getExtraPath())]


class SASAStructureViewer(ChimeraAttributeViewer):
    """ Viewer for attribute SASA of an AtomStruct.
      Includes visualization in chimera and in histograms"""
    _targets = [ProtCalculateSASA]
    _label = 'Atomic structure attributes viewer'

    def __init__(self, **kwargs):
      super().__init__(**kwargs)

    def _defineParams(self, form):
      super()._defineParams(form)
      # Overwrite defaults
      from pwem.wizards.wizard import ColorScaleWizardBase
      group = form.addGroup('Color settings')
      ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=200, defaultIntervals=21,
                                                  defaultColorMap='RdBu')


