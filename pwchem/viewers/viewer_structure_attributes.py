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
import numpy as np
import matplotlib.pyplot as plt

from pyworkflow.protocol import params
from pwem.viewers import ChimeraAttributeViewer
from pwchem.protocols import ProtCalculateSASA, ProtSeqCalculateConservation
from deeploc.protocols import ProtDeepLoc
from pwchem.viewers.viewers_sequences import SequenceAliView

def plotSequenceAttribute(attrValues, attrName='Attribute', thres=None):
    attrValues = list(map(float, attrValues))
    maxY = max(attrValues)
    xs = np.arange(len(attrValues))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(xs, attrValues)
    yloc = plt.MaxNLocator(10)
    ax.yaxis.set_major_locator(yloc)
    ax.set_ylim(0, maxY + maxY / 10)
    if thres:
        ax.axhline(y=thres, color='r', linestyle='-', linewidth=1)

    plt.xlabel("Sequence position")
    plt.ylabel("{} value".format(attrName))
    plt.title('{} values along sequence'.format(attrName))

    plt.show()

class ConservationViewer(ChimeraAttributeViewer):
    """ Viewer for attribute conservation of an AtomStruct.
      Includes visualization in chimera and in histograms"""
    _targets = [ProtSeqCalculateConservation]
    _label = 'Sequence conservation analysis viewer'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self, form):
        form.addSection(label='Visualization of sequence conservation')
        form.addParam('viewSequence', params.LabelParam, label='View sequence: ',
                       help='View output sequence')
        form.addParam('viewConservation', params.LabelParam,
                      label='Display conservation over sequence: ',
                      help='Display a graph witht the values of the selected attribute over the sequence.')

        if hasattr(self.protocol, 'inputAS') and getattr(self.protocol, 'inputAS').get():
            super()._defineParams(form)
            # Overwrite defaults
            from pwem.wizards.wizard import ColorScaleWizardBase
            group = form.addGroup('Color settings')
            ColorScaleWizardBase.defineColorScaleParams(group, defaultLowest=0, defaultHighest=1, defaultIntervals=21,
                                                        defaultColorMap='RdBu_r')

    def _getVisualizeDict(self):
        visDic = {'viewSequence': self._showSequence, 'viewConservation': self._showConservation}
        if hasattr(self.protocol, 'inputAS') and getattr(self.protocol, 'inputAS').get():
            visDic.update(super()._getVisualizeDict())
        return visDic

    def _showSequence(self, paramName=None):
        obj = self.protocol.outputSequence
        outPath = os.path.abspath(self.protocol._getExtraPath('viewSequences_{}.fasta'.
                                                              format(obj.getId())))
        obj.exportToFile(outPath)
        return [SequenceAliView([outPath], cwd=self.protocol._getExtraPath())]

    def _showConservation(self, paramName=None):
        prot = self.protocol
        attrValues = list(prot.getConsDic().values())
        plotSequenceAttribute(attrValues, attrName=prot.getEnumText('method'))



class SASAStructureViewer(ChimeraAttributeViewer):
    """Viewer for residue attributes stored in an AtomStruct.
    Includes structure visualization in ChimeraX and attribute
    visualization as histograms or accessibility sequence regions.
    """

    _targets = [ProtCalculateSASA, ProtDeepLoc]
    _label = 'Accessibility viewer'

    def __init__(self, **kwargs):
        super().__init__(**kwargs)

    def _defineParams(self, form):
        if hasattr(self.protocol, 'outputSequence'):
            if isinstance(self.protocol, ProtDeepLoc):
                sectionLabel = 'Visualization of DeepLoc residue importance'
                sequenceLabel = 'View sequence: '
                attributeLabel = 'Display DeepLoc residue importance: '
            else:
                sectionLabel = 'Visualization of sequence SASA'
                sequenceLabel = 'View sequence: '
                attributeLabel = 'Display SASA over sequence: '

            form.addSection(label=sectionLabel)

            form.addParam(
                'viewSequence',
                params.LabelParam,
                label=sequenceLabel,
                help='View output sequence'
            )

            form.addParam(
                'viewSequenceAttribute',
                params.LabelParam,
                label=attributeLabel,
                help='Display the residue attribute over the sequence.'
            )
        if hasattr(self.protocol, 'outputSequences'):
            if isinstance(self.protocol, ProtDeepLoc):
                sectionLabel = 'Visualization of DeepLoc residue importance'
            else:
                sectionLabel = 'Visualization of sequence SASA'

            form.addSection(label=sectionLabel)

            form.addParam(
                'viewSequences',
                params.LabelParam,
                label='View all sequences: ',
                help='View all output sequences.'
            )
            form.addParam(
                'viewSequencesAttribute',
                params.LabelParam,
                label='Display residue attribute for all sequences: ',
                help='Display the residue attribute for all output sequences.'
            )

        if hasattr(self.protocol, 'outputAtomStruct'):
            super()._defineParams(form)
            from pwem.wizards.wizard import ColorScaleWizardBase
            group = form.addGroup('Color settings')

            if isinstance(self.protocol, ProtDeepLoc):
                lowest = 0
                highest = 1
            else:
                lowest = 0
                highest = 200

            ColorScaleWizardBase.defineColorScaleParams(
                group,
                defaultLowest=lowest,
                defaultHighest=highest,
                defaultIntervals=21,
                defaultColorMap='RdBu'
            )
        if hasattr(self.protocol, 'outputAtomStructs'): #todo

            form.addSection(
                label='Visualization of structure set'
            )

            form.addParam(
                'viewAtomStructs',
                params.LabelParam,
                label='View all structures: ',
                help='View all output structures in ChimeraX.'
            )

    def _getVisualizeDict(self):
        visDic = {}
        if hasattr(self.protocol, 'outputSequence'):

            visDic.update({
                'viewSequence':
                    self._showSequenceAttrs,

                'viewSequenceAttribute':
                    self._showSequenceAttribute
            })
        if hasattr(self.protocol, 'outputSequences'):

            visDic.update({
                'viewSequences':
                    self._showSequences,

                'viewSequencesAttribute':
                    self._showSequencesAttribute
            })
        if hasattr(self.protocol, 'outputAtomStruct'):

            visDic.update(
                super()._getVisualizeDict()
            )
        if hasattr(self.protocol, 'outputAtomStructs'):

            visDic.update({
                'viewAtomStructs':
                    self._showAtomStructs
            })

        return visDic

    def _showSequenceAttrs(self, paramName=None):

        obj = self.protocol.outputSequence

        outPath = os.path.abspath(
            self.protocol._getExtraPath(
                'viewSequences_{}.fasta'.format(obj.getId())
            )
        )

        obj.exportToFile(outPath)

        return [
            SequenceAliView(
                [outPath],
                cwd=self.protocol._getExtraPath()
            )
        ]

    def _showSequenceAttribute(self, paramName=None):

        prot = self.protocol
        attrName = prot._ATTRNAME

        attrDic = prot.outputSequence.getAttributesDic()

        plotSequenceAttribute(
            attrDic[attrName],
            attrName=attrName
        )

        # ================================================================
        # Set of Sequences
        # ================================================================

    def _showSequences(self, paramName=None):

        outPath = os.path.abspath(
            self.protocol._getExtraPath('viewSequences.fasta')
        )

        if os.path.exists(outPath):
            os.remove(outPath)

        self.protocol.outputSequences.exportToFile(outPath)

        return [
            SequenceAliView(
                [outPath],
                cwd=self.protocol._getExtraPath()
            )
        ]

    def _showSequencesAttribute(self, paramName=None):

        attrName = self.protocol._ATTRNAME

        for sequence in self.protocol.outputSequences:
            attrDic = sequence.getAttributesDic()

            if attrName not in attrDic:
                continue

            plotSequenceAttribute(
                attrDic[attrName],
                attrName='{} - {}'.format(
                    sequence.getSeqName(),
                    attrName
                )
            )

        # ================================================================
        # Set of AtomStructs
        # ================================================================

    def _showAtomStructs(self, paramName=None):

        modelFiles = []

        for atomStruct in self.protocol.outputAtomStructs:
            modelFiles.append(
                os.path.abspath(
                    atomStruct.getFileName()
                )
            )

        return [
            ChimeraView(
                modelFiles,
                cwd=self.protocol._getExtraPath()
            )
        ]

