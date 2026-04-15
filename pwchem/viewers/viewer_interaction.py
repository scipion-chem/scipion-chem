# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import numpy as np
import matplotlib.pyplot as plt
from tkinter.messagebox import askokcancel

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.protocol import params

from pwem.objects import SetOfAtomStructs
from pwem.viewers import EmPlotter
from pwchem.viewers.viewers_sequences import heatmap, annotateHeatmap
from pwem.protocols import ProtSubSet


class BaseInteractionViewer(ProtocolViewer):
    """Generic viewer for interaction matrices"""

    def _getData(self):
        """Return the full interaction dictionary"""
        raise NotImplementedError

    def _getEntityNames(self, data):
        """Return (entities1, entities2, scoreTypes)"""
        raise NotImplementedError

    def _getLabels(self):
        """Return labels for UI (entity1, entity2, title)"""
        return "Entity1", "Entity2", "Score"

    def _defineParams(self, form):
        data = self._getData()

        if data:
            self._defineInteractionParams(form, data)
        else:
            form.addSection(label='Viewer')
            form.addInfo('No interaction data found.')

    def _defineInteractionParams(self, form, data):
        ent1, ent2, scoreTypes = self._getEntityNames(data)
        label1, label2, _ = self._getLabels()

        form.addSection(label='Interaction viewer')
        sGroup = form.addGroup('Filters')

        sGroup.addParam('chooseEnt1', params.EnumParam,
                        label=f'{label1}: ',
                        choices=['All'] + ent1, default=0)

        sGroup.addParam('chooseEnt2', params.EnumParam,
                        label=f'{label2}: ',
                        choices=['All'] + ent2, default=0)

        sGroup.addParam('chooseScore', params.EnumParam,
                        label='Score type: ',
                        choices=scoreTypes, default=0)

        sGroup.addParam('scThres', params.FloatParam,
                        label='Score threshold: ', default=0.0)

        hGroup = form.addGroup('Visualizations')
        hGroup.addParam('displayHeatMap', params.LabelParam, label='Heatmap')
        hGroup.addParam('displayHistogram', params.LabelParam, label='Histogram')
        hGroup.addParam('intervals', params.IntParam, default=10, label='Number of bins')

    def _getVisualizeDict(self):
        return {
            'displayHeatMap': self._viewHeatMap,
            'displayHistogram': self._viewHistogram
        }

    def _getFilteredData(self, filt1, filt2, filtScore):
        data = self._getData()

        ent1 = [filt1] if filt1 != 'All' else sorted(data.keys())

        allEnt2 = set()
        for e in data:
            allEnt2.update(data[e].keys())

        ent2 = [filt2] if filt2 != 'All' else sorted(allEnt2)

        matrix = np.zeros((len(ent1), len(ent2)))
        thres = self.scThres.get()

        for i, e1 in enumerate(ent1):
            for j, e2 in enumerate(ent2):
                val = data.get(e1, {}).get(e2, {}).get(filtScore, np.nan)
                matrix[i, j] = val if val >= thres else np.nan

        return matrix, ent1, ent2, filtScore

    def _viewHeatMap(self, paramName=None):
        f1 = self.getEnumText('chooseEnt1')
        f2 = self.getEnumText('chooseEnt2')
        fScore = self.getEnumText('chooseScore')

        matrix, e1, e2, label = self._getFilteredData(f1, f2, fScore)

        if len(e1) * len(e2) > 500:
            if not askokcancel("Large dataset", "Heatmap may be large. Continue?"):
                return

        fig, ax = plt.subplots(figsize=(10, 8))
        im, _ = heatmap(matrix, e1, e2, ax=ax,
                        cmap="YlOrRd",
                        cbarLabel=f"{label}")

        if len(e1) * len(e2) < 150:
            annotateHeatmap(im, valfmt="{x:.2f}")

        plt.show()

    def _viewHistogram(self, paramName=None):
        f1 = self.getEnumText('chooseEnt1')
        f2 = self.getEnumText('chooseEnt2')
        fScore = self.getEnumText('chooseScore')

        matrix, _, _, label = self._getFilteredData(f1, f2, fScore)
        scores = matrix.flatten()
        scores = scores[~np.isnan(scores)]

        plotter = EmPlotter(x=1, y=1, windowTitle='Score Distribution')
        ax = plotter.createSubPlot(f"{label} distribution", label, "Count")

        ax.hist(scores, bins=self.intervals.get())
        ax.grid(True)

        plotter.show()