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
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
from pwem.wizards.wizard import ColorScaleWizardBase

def plotSequenceAttribute(attrValues, attrName='Attribute', thres=None):
    attrValues = list(map(float, attrValues))
    maxY = max(attrValues)
    xs = np.arange(len(attrValues))

    _, ax = plt.subplots()
    ax.bar(xs, attrValues)
    ax.yaxis.set_major_locator(plt.MaxNLocator(10))
    ax.set_ylim(0, maxY + maxY / 10)
    if thres:
        ax.axhline(
            y=thres,
            color='r',
            linestyle='-',
            linewidth=1
        )

    ax.set_xlabel('Sequence position')
    ax.set_ylabel('{} value'.format(attrName))
    ax.set_title('{} values along sequence'.format(attrName))

    plt.show()


def plotInteractive(data, histogram=False):
    if not data:
        return

    currentIndex = [0]

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.20)

    def drawPlot():
        ax.clear()

        name, attrName, attrValues = data[currentIndex[0]]
        attrValues = list(map(float, attrValues))

        if histogram:
            ax.hist(attrValues, bins=20)
            ax.set_xlabel(attrName)
            ax.set_ylabel('Frequency')
        else:
            xs = np.arange(1, len(attrValues) + 1)
            ax.bar(xs, attrValues)
            ax.set_xlabel('Residue')
            ax.set_ylabel(attrName)
            ax.set_xlim(
                0.5,
                len(attrValues) + 0.5
            )

        ax.set_title(
            '{} - {}'.format(name, attrName)
        )

        fig.canvas.draw_idle()

    previousAx = plt.axes([0.25, 0.05, 0.20, 0.075])
    nextAx = plt.axes([0.55, 0.05, 0.20, 0.075])

    previousButton = Button(previousAx, 'Previous')
    nextButton = Button(nextAx, 'Next')

    def previous(event):
        if currentIndex[0] > 0:
            currentIndex[0] -= 1
            drawPlot()

    def next(event):
        if currentIndex[0] < len(data) - 1:
            currentIndex[0] += 1
            drawPlot()

    previousButton.on_clicked(previous)
    nextButton.on_clicked(next)

    drawPlot()
    plt.show()


def plotSequenceAttributesInteractive(sequenceData):
    plotInteractive(sequenceData)


def plotAtomStructAttributesInteractive(structureData):
    plotInteractive(structureData, histogram=True)


def plotAtomStructSequenceAttributesInteractive(structureData):
    plotInteractive(structureData)

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
    _targets = [ProtCalculateSASA, ProtDeepLoc]
    _label = 'Accessibility viewer'

    def getAtomStructObject(self):
        if hasattr(self.protocol, 'outputAtomStructs'):
            return self._getSelectedAtomStruct()
        return super().getAtomStructObject()

    def getEnumText(self, paramName):
        if (
            paramName == 'attrName'
            and hasattr(self.protocol, 'outputAtomStructs')
        ):
            return self.protocol._ATTRNAME
        return super().getEnumText(paramName)

    def _isDeepLoc(self):
        return isinstance(self.protocol, ProtDeepLoc)

    def _getColorScale(self):
        return (0, 1) if self._isDeepLoc() else (0, 200)

    def _defineColorScale(self, form):
        lowest, highest = self._getColorScale()

        group = form.addGroup('Color settings')

        ColorScaleWizardBase.defineColorScaleParams(
            group,
            defaultLowest=lowest,
            defaultHighest=highest,
            defaultIntervals=21,
            defaultColorMap='RdBu'
        )

    def _defineParams(self, form):
        if hasattr(self.protocol, 'outputSequence'):
            if self._isDeepLoc():
                sectionLabel = 'Visualization of DeepLoc residue importance'
                attributeLabel = 'Display DeepLoc residue importance: '
            else:
                sectionLabel = 'Visualization of sequence SASA'
                attributeLabel = 'Display SASA over sequence: '

            form.addSection(label=sectionLabel)

            form.addParam(
                'viewSequence',
                params.LabelParam,
                label='View sequence: ',
                help='View output sequence'
            )

            form.addParam(
                'viewSequenceAttribute',
                params.LabelParam,
                label=attributeLabel,
                help='Display the residue attribute over the sequence.'
            )

        if hasattr(self.protocol, 'outputSequences'):
            sectionLabel = (
                'Visualization of DeepLoc residue importance'
                if self._isDeepLoc()
                else 'Visualization of sequence SASA'
            )

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
            self._defineColorScale(form)

        if hasattr(self.protocol, 'outputAtomStructs'):
            self._defineAtomStructSetParams(form)

    def _defineAtomStructSetParams(self, form):
        form.addSection(label='Visualization of structure set')

        self._structureNames = [
            os.path.splitext(
                os.path.basename(atomStruct.getFileName())
            )[0]
            for atomStruct in self.protocol.outputAtomStructs
        ]

        form.addParam(
            'atomStruct',
            params.EnumParam,
            choices=self._structureNames,
            default=0,
            label='Structure: ',
            help='Select the structure to visualize.'
        )

        form.addParam(
            'viewAtomStruct',
            params.LabelParam,
            label='View structure: ',
            help='View the selected structure in ChimeraX.'
        )

        form.addParam(
            'viewAtomStructAttribute',
            params.LabelParam,
            label='Display attribute histogram: ',
            help='Display the residue attribute distribution '
                 'for the selected structure.'
        )

        form.addParam(
            'chain_name',
            params.StringParam,
            default='A',
            allowsNull=True,
            label='Chain of interest: ',
            help='Specify the chain of interest (e.g. A).'
        )

        form.addParam(
            'viewAtomStructSequenceAttribute',
            params.LabelParam,
            label='Display attribute over sequence: ',
            help='Display the residue attribute over the sequence '
                 'for the selected structure.'
        )

        self._defineColorScale(form)

    def _getVisualizeDict(self):
        visDic = {}

        if hasattr(self.protocol, 'outputSequence'):
            visDic.update({
                'viewSequence': self._showSequenceAttrs,
                'viewSequenceAttribute': self._showSequenceAttribute
            })

        if hasattr(self.protocol, 'outputSequences'):
            visDic.update({
                'viewSequences': self._showSequences,
                'viewSequencesAttribute': self._showSequencesAttribute
            })

        if hasattr(self.protocol, 'outputAtomStruct'):
            visDic.update(
                super()._getVisualizeDict()
            )

        if hasattr(self.protocol, 'outputAtomStructs'):
            visDic.update({
                'viewAtomStruct': self._showAtomStruct,
                'viewAtomStructAttribute': self._showAtomStructAttribute,
                'viewAtomStructSequenceAttribute':
                    self._showAtomStructSequenceAttribute
            })

        return visDic

    # ------------------------------------------------------------------
    # Single sequence
    # ------------------------------------------------------------------

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
        attrName = self.protocol._ATTRNAME
        attrDic = self.protocol.outputSequence.getAttributesDic()

        plotSequenceAttribute(
            attrDic[attrName],
            attrName=attrName
        )

    # ------------------------------------------------------------------
    # Set of sequences
    # ------------------------------------------------------------------

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

        sequenceData = [
            (
                sequence.getSeqName(),
                attrName,
                sequence.getAttributesDic()[attrName]
            )
            for sequence in self.protocol.outputSequences
            if attrName in sequence.getAttributesDic()
        ]

        plotSequenceAttributesInteractive(sequenceData)

    # ------------------------------------------------------------------
    # Set of AtomStructs
    # ------------------------------------------------------------------

    def _getSelectedAtomStruct(self):
        selectedName = self._structureNames[self.atomStruct.get()]

        for atomStruct in self.protocol.outputAtomStructs:
            name = os.path.splitext(
                os.path.basename(atomStruct.getFileName())
            )[0]

            if name == selectedName:
                return atomStruct

        return None

    def _showAtomStruct(self, paramName=None):
        atomStruct = self._getSelectedAtomStruct()

        print(
            'SELECTED STRUCTURE:',
            atomStruct.getFileName()
        )

        self._atomStruct = atomStruct

        intermediateFile = self.protocol._getExtraPath(
            'chimeraAttribute_{}.cif'.format(
                self.protocol._ATTRNAME
            )
        )

        if os.path.exists(intermediateFile):
            os.remove(intermediateFile)

        return self._showChimera(paramName)

    def _showAtomStructAttribute(self, paramName=None):
        return self._showHistogram(paramName)

    def _showAtomStructSequenceAttribute(self, paramName=None):
        return self._showSequence(paramName)
