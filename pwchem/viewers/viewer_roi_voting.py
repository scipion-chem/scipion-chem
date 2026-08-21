# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
# *          Judith Maestro (judith.maestro@cnb.csic.es)
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

import pyworkflow.protocol.params as params
from pwem.viewers import ChimeraAttributeViewer
from pwchem.protocols.VirtualDrugScreening.protocol_structROI_voting import ProtROIVoting
from pwchem.viewers.viewer_structure_attributes import plotSequenceAttribute


class ViewerROIVoting(ChimeraAttributeViewer):
    _label = 'Viewer ROI Voting'
    _targets = [ProtROIVoting]

    def _getSequenceOutputsDic(self):
        '''{chainId: SequenceChem output} for every per-chain sequence output of the protocol'''
        seqOutputs = {}
        for attrName in dir(self.protocol):
            if attrName.startswith('outputSequence_'):
                seqOutputs[attrName[len('outputSequence_'):]] = getattr(self.protocol, attrName)
        return seqOutputs

    def _defineParams(self, form):
        seqOutputs = self._getSequenceOutputsDic()
        if seqOutputs:
            form.addSection(label='Sequence')
            if len(seqOutputs) > 1:
                form.addParam('seqChain', params.EnumParam, choices=sorted(seqOutputs.keys()), default=0,
                              label='Chain: ', help='Chain whose sequence frequency will be displayed.')
            form.addParam('viewFrequency', params.LabelParam,
                          label='Display frequency over sequence: ',
                          help='Display a bar chart with frequency values per residue.')
            form.addParam('viewSequenceColored', params.LabelParam,
                          label='Display sequence colored by frequency: ',
                          help='Shows the sequence with each AA colored white to red by frequency.')

        if hasattr(self.protocol, self.protocol._OUTNAME):
            super()._defineParams(form)

    def _getVisualizeDict(self):
        visDic = {}
        if self._getSequenceOutputsDic():
            visDic['viewFrequency'] = self._showFrequency
            visDic['viewSequenceColored'] = self._showSequenceColored
        if hasattr(self.protocol, self.protocol._OUTNAME):
            visDic.update(super()._getVisualizeDict())
        return visDic

    def _getSelectedSequenceOutput(self):
        seqOutputs = self._getSequenceOutputsDic()
        if len(seqOutputs) > 1:
            chainId = sorted(seqOutputs.keys())[self.seqChain.get()]
        else:
            chainId = next(iter(seqOutputs.keys()))
        return seqOutputs[chainId]

    def _showFrequency(self, paramName=None):
        attrDic = self._getSelectedSequenceOutput().getAttributesDic()
        plotSequenceAttribute(attrDic['frequency'], attrName='ROI Voting Frequency')

    def _showSequenceColored(self, paramName=None):
        seqOutput = self._getSelectedSequenceOutput()
        attrDic = seqOutput.getAttributesDic()
        freqValues = list(map(float, attrDic['frequency']))
        seqStr = seqOutput.getSequence()
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
      