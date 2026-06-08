# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os
import webbrowser
import numpy as np
import matplotlib.pyplot as plt

import pyworkflow.viewer as pwviewer
import pyworkflow.protocol.params as params
from pwchem.viewers.viewers_data import PyMolViewer
from pwchem.objects import SetOfStructROIs
from pwchem.protocols.VirtualDrugScreening.protocol_structROI_voting import ProtROIVoting


def plot_residue_attribute(attrValues, attrName='Attribute'):
    attrValues = list(map(float, attrValues))
    maxY = max(attrValues)
    xs = np.arange(len(attrValues))

    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.bar(xs, attrValues, color='steelblue')
    yloc = plt.MaxNLocator(10)
    ax.yaxis.set_major_locator(yloc)
    ax.set_ylim(0, maxY + maxY / 10)

    plt.xlabel("Sequence position")
    plt.ylabel(f"{attrName} value")
    plt.title(f'{attrName} values along sequence')
    plt.show()


class ViewerROIVoting(pwviewer.ProtocolViewer):
    _label = 'Viewer ROI Voting'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [ProtROIVoting]

    def _defineParams(self, form):
        form.addSection(label='3D Visualization')
        form.addParam('displayPyMol', params.LabelParam,
                      label='Display residues colored by frequency (white to red)',
                      help='Colors each residue by its voting frequency percentage in PyMol.')

        if hasattr(self.protocol, 'outputSequence'):
            form.addSection(label='Sequence')
            form.addParam('viewFrequency', params.LabelParam,
                          label='Display frequency over sequence: ',
                          help='Display a bar chart with frequency values per residue.')
            form.addParam('viewSequenceColored', params.LabelParam,
                          label='Display sequence colored by frequency: ',
                          help='Shows the sequence with each AA colored white to red by frequency.')

    def _getVisualizeDict(self):
        visDic = {
            'displayPyMol': self._showColorMap,
        }
        if hasattr(self.protocol, 'outputSequence'):
            visDic['viewFrequency'] = self._showFrequency
            visDic['viewSequenceColored'] = self._showSequenceColored
        return visDic

    def _validate(self):
        return []

    def _showColorMap(self, paramName=None):
        if isinstance(self.protocol, SetOfStructROIs):
            roiSet = self.protocol
        elif hasattr(self.protocol, 'outputStructROIs'):
            roiSet = self.protocol.outputStructROIs
        elif hasattr(self.protocol, 'outputFilteredROIs'):
            roiSet = self.protocol.outputFilteredROIs
        else:
            print("NO ROISET FOUND")
            return []

        projectPath = self.getProject().getPath()
        proteinFile = os.path.join(projectPath, roiSet.getProteinFile())
        outDir = roiSet.getSetDir()
        pmlFile = os.path.join(outDir, 'colormap.pml')
        pmlFile = self._generatePmlScript(roiSet, proteinFile, pmlFile)
        pymolV = PyMolViewer(project=self.getProject())
        return pymolV._visualize(pmlFile, cwd=os.path.dirname(pmlFile))

    def _generatePmlScript(self, roiSet, proteinFile, pmlFile):
        with open(pmlFile, 'w') as f:
            f.write(f'load {proteinFile}\n')
            f.write('hide everything\n')
            f.write('show cartoon\n')
            f.write('color white\n')

            for roi in roiSet:
                residue = roi._contactResidues.get()
                perc = roi._percentage.get() if hasattr(roi, '_percentage') else 0.0

                try:
                    chain, resnum = residue.split('_')
                    intensity = max(0.0, min(1.0, perc / 100.0))
                    r = 1.0
                    g = 1.0 - intensity
                    b = 1.0 - intensity
                    f.write(f'set_color col_{resnum}, [{r:.3f}, {g:.3f}, {b:.3f}]\n')
                    f.write(f'color col_{resnum}, chain {chain} and resi {resnum}\n')
                except Exception:
                    continue

        return pmlFile

    def _showFrequency(self, paramName=None):
        attrDic = self.protocol.outputSequence.getAttributesDic()
        plot_residue_attribute(attrDic['frequency'], attrName='ROI Voting Frequency')

    def _showSequenceColored(self, paramName=None):
        attrDic = self.protocol.outputSequence.getAttributesDic()
        freqValues = list(map(float, attrDic['frequency']))
        seqStr = self.protocol.outputSequence.getSequence()
        maxFreq = max(freqValues) if max(freqValues) > 0 else 1

        htmlContent = """
        <html><head><title>ROI Voting Sequence</title>
        <style>
            body { font-family: monospace; padding: 20px; background: #1e1e1e; color: white; }
            h2 { color: #ccc; }
            .seq { display: flex; flex-wrap: wrap; gap: 2px; margin-top: 20px; }
            .aa { width: 28px; height: 36px; display: flex; flex-direction: column;
                  align-items: center; justify-content: center; font-size: 13px;
                  font-weight: bold; border-radius: 3px; color: black; }
            .pos { font-size: 8px; color: #333; }
            .legend { margin-top: 20px; display: flex; align-items: center; gap: 10px; }
            .grad { width: 200px; height: 20px;
                    background: linear-gradient(to right, white, red); border: 1px solid #555; }
        </style></head><body>
        <h2>ROI Voting Frequency along sequence</h2>
        <div class='seq'>
        """

        for i, (aa, freq) in enumerate(zip(seqStr, freqValues)):
            intensity = freq / maxFreq
            r = 255
            g = int(255 * (1 - intensity))
            b = int(255 * (1 - intensity))
            htmlContent += f"""
            <div class='aa' style='background: rgb({r},{g},{b});' title='Position {i+1}, freq={int(freq)}'>
                {aa}<span class='pos'>{i+1}</span>
            </div>"""

        htmlContent += """
        </div>
        <div class='legend'>
            <span>0</span>
            <div class='grad'></div>
            <span>max</span>
        </div>
        </body></html>"""

        outPath = os.path.abspath(self.protocol._getExtraPath('sequence_colored.html'))
        with open(outPath, 'w') as f:
            f.write(htmlContent)
        webbrowser.open(f"file://{outPath}")
        return []
      