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
# * e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

"""
This protocol visualizes the results from 'ProtROIVoting'
on a protein sequence using a color map depending on the
frequency/percentage of each residue position.
"""

import pyworkflow.object as pwobj
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfROIColorMap, ROIColorItem


class ProtROIColorMap(EMProtocol):
    _label = 'ROI Color Map (visual sequence intensity)'

    # ------------- PARAMS -------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVotes', params.PointerParam,
                      label='ROI Voting Output',
                      pointerClass='SetOfROIVotes',
                      help='Select the ROI voting result to visualize.')
        form.addParam('proteinSequence', params.StringParam,
                      label='Protein sequence (AA)',
                      help='Paste the full amino acid sequence here.')
        form.addParam('chainPrefix', params.StringParam, default='A_',
                      label='Residue prefix (e.g., A_ for chain A)',
                      help='Prefix used in ROI residue naming (e.g., A_123).')

    # ------------- STEPS -------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createColorMapStep')

    # ------------- MAIN STEP -------------
    def createColorMapStep(self):
        votesSet = self.inputVotes.get()
        sequence = self.proteinSequence.get().strip()
        prefix = self.chainPrefix.get()

        # Map each residue position to its voting percentage
        res_dict = {}
        for roi in votesSet:
            try:
                resName = roi._residue.get()
                if resName.startswith(prefix):
                    pos = int(resName.split('_')[1])
                    res_dict[pos] = roi._percentage.get()
            except Exception:
                continue

        # Color scaling
        max_color = 100.0
        outSet = SetOfROIColorMap(filename=self._getPath('ColorMap.sqlite'))

        # Fill sequence items
        for i, aa in enumerate(sequence, start=1):
            perc = res_dict.get(i, 0.0)
            color = self._mapColor(perc / max_color)
            item = ROIColorItem()
            item._residue = pwobj.String(f"{prefix}{i}")
            item._frequency = pwobj.Integer(0)
            item._percentage = pwobj.Float(round(perc, 2))
            item._color = pwobj.String(color)
            outSet.append(item)

        # Register Scipion output (visible in GUI)
        if len(outSet) > 0:
            self._defineOutputs(outputColorMap=outSet)
            self._defineSourceRelation(self.inputVotes, self.outputColorMap)
            self.info(f"Color map created for {len(sequence)} residues.")
        else:
            self.warning("No residues were mapped to the sequence.")

    # ------------- COLOR UTILITY -------------
    def _mapColor(self, intensity):
        """
        Map a normalized value in [0,1] to a color (blue→red gradient).
        """
        r = int(255 * intensity)
        g = int(255 * (1 - intensity))
        b = 255 - r
        return f'#{r:02x}{g:02x}{b:02x}'

    # ------------- INFO -------------
    def _summary(self):
        seq = self.proteinSequence.get()
        n = len(seq) if seq else 0
        return [f"Generated color map for {n} residues.",
                "Color intensity corresponds to residue voting percentage "
                "(blue→red gradient)."]

    def _methods(self):
        return ["Each residue receives a color (blue→red) proportional to its voting percentage."]
