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


class OmniBindPredictionViewer(ProtocolViewer):
    _label = 'Interactions Viewer'
    _targets = [SetOfAtomStructs]
    _environments = [DESKTOP_TKINTER]

    def _defineParams(self, form):
        self.structSet = self._getStructSet()

        if self.structSet and self.checkIfInteractions():
            self._defineInteractionParams(form)
        else:
            form.addSection(label='Viewer')
            form.addInfo('No interaction data found for these structures.')

    def _defineInteractionParams(self, form):
        scoresFile = self.structSet._interactScoresFile.get()
        with open(scoresFile, 'r') as f:
            data = json.load(f)

        protNames = sorted(list(data.keys()))
        molNames = sorted(list(data[protNames[0]].keys())) if protNames else []
        scoreTypes = sorted(list(data[protNames[0]][molNames[0]].keys())) if molNames else []

        form.addSection(label='Interaction Analysis')

        sGroup = form.addGroup('Filters')
        sGroup.addParam('chooseStruct', params.EnumParam, label='Protein structure: ',
                        choices=['All'] + protNames, default=0,
                        help='Filter results by a specific protein structure')

        sGroup.addParam('chooseMol', params.EnumParam, label='Small molecule: ',
                        choices=['All'] + molNames, default=0,
                        help='Filter results by a specific molecule')

        sGroup.addParam('chooseScore', params.EnumParam, label='Score type: ',
                        choices=scoreTypes, default=0)

        sGroup.addParam('scThres', params.FloatParam, label='Score threshold: ', default=0.0,
                        help='Display results above this value')

        hGroup = form.addGroup('Visualizations')
        hGroup.addParam('displayHeatMap', params.LabelParam, label='Show Heatmap',
                        help='Matrix of protein-ligand affinities')
        hGroup.addParam('displayHistogram', params.LabelParam, label='Show Histogram',
                        help='Distribution of predicted scores')
        hGroup.addParam('intervals', params.IntParam, default=10, label="Bins",
                        expertLevel=params.LEVEL_ADVANCED)

    def _getVisualizeDict(self):
        return {
            'displayHeatMap': self._viewHeatMap,
            'displayHistogram': self._viewHistogram
        }

    def _viewHeatMap(self, paramName=None):
        filtProt = self.getEnumText('chooseStruct')
        filtMol = self.getEnumText('chooseMol')
        filtScore = self.getEnumText('chooseScore')

        intAr, prots, mols, scoreLabel = self._getFilteredData(filtProt, filtMol, filtScore)

        if len(prots) * len(mols) > 500:
            if not askokcancel("Large dataset", "The heatmap might be too large. Continue?"):
                return

        fig, ax = plt.subplots(figsize=(10, 8))
        im, _ = heatmap(intAr, prots, mols, ax=ax, cmap="YlOrRd",
                        cbarLabel=f"OmniBind {scoreLabel}")

        if len(prots) * len(mols) < 150:
            annotateHeatmap(im, valfmt="{x:.2f}")

        fig.tight_layout()
        plt.show()

    def _viewHistogram(self, paramName=None):
        filtProt = self.getEnumText('chooseStruct')
        filtMol = self.getEnumText('chooseMol')
        filtScore = self.getEnumText('chooseScore')

        intAr, _, _, scoreLabel = self._getFilteredData(filtProt, filtMol, filtScore)
        scoreList = intAr.flatten()
        scoreList = scoreList[~np.isnan(scoreList)]

        plotter = EmPlotter(x=1, y=1, windowTitle='OmniBind Score Distribution')
        a = plotter.createSubPlot(f"Distribution of {scoreLabel}", 'Predicted Affinity', 'Frequency')

        a.hist(scoreList, bins=self.intervals.get(), color='skyblue', edgecolor='black', alpha=0.7)
        a.grid(axis='y', linestyle='--', alpha=0.6)
        plotter.show()

    # ----------- UTILS -----------

    def _getStructSet(self):
        if hasattr(self.protocol, 'outputAtomStructs'):
            return self.protocol.outputAtomStructs
        return self.protocol

    def checkIfInteractions(self):
        return hasattr(self.structSet, '_interactScoresFile') and self.structSet._interactScoresFile.get() is not None

    def _getFilteredData(self, filtProt, filtMol, filtScore):
        with open(self.structSet._interactScoresFile.get(), 'r') as f:
            fullData = json.load(f)

        prots = [filtProt] if filtProt != 'All' else sorted(fullData.keys())
        allMols = set()
        for p in fullData:
            allMols.update(fullData[p].keys())
        mols = [filtMol] if filtMol != 'All' else sorted(list(allMols))

        matrix = np.zeros((len(prots), len(mols)))
        thres = self.scThres.get()

        for i, p in enumerate(prots):
            for j, m in enumerate(mols):
                val = fullData.get(p, {}).get(m, {}).get(filtScore, np.nan)
                matrix[i, j] = val if val >= thres else np.nan

        return matrix, prots, mols, filtScore