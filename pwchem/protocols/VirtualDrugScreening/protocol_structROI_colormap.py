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
This protocol visualizes the results from "protocol_structROI_voting" on a protein sequence
using a color map depending on the frequency/percentage of each residue position.
"""

import pyworkflow.object as pwobj
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import EMSet
from pwchem.objects import SetOfROIColorMap, ROIColorItem


class ProtROIColorMap(EMProtocol):
    _label = 'ROI Color Map (visual sequence intensity)'

    # ---------------- DEFINE PARAMS ----------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputVotes', params.PointerParam,
                      label='ROI Voting Output',
                      pointerClass='EMSet',
                      help='Select the ROI voting result to visualize (currently stored as EMSet).')
        form.addParam('proteinSequence', params.StringParam,
                      label='Protein sequence (AA)',
                      help='Paste the full amino acid sequence here.')
        form.addParam('chainPrefix', params.StringParam, default='A_',
                      label='Residue prefix (e.g., A_ for chain A)',
                      help='Prefix used in ROI residue naming (e.g., A_123).')

    # ---------------- INSERT STEPS ----------------
    def _insertAllSteps(self):
        self._insertFunctionStep('createColorMapStep')

    # ---------------- MAIN STEP ----------------
    def createColorMapStep(self):
        votesSet = self.inputVotes.get()
        sequence = (self.proteinSequence.get() or '').strip()
        prefix = (self.chainPrefix.get() or 'A_').strip()

        if not sequence:
            self.warning("Protein sequence is empty. Please provide a valid amino acid sequence.")
            return

        # ---------------- Read voting set into dictionaries ----------------
        # position (1-based) -> frequency / percentage
        res_freq = {}
        res_perc = {}

        for roi in votesSet:
            # Residue name can be stored as _residue or _sequence (depending on upstream object type)
            resObj = getattr(roi, '_residue', None) or getattr(roi, '_sequence', None)
            resName = resObj.get() if resObj is not None else ''
            if not resName:
                continue

            # Determine position:
            # Prefer _roiIdx if present, otherwise parse from "A_123"
            idxObj = getattr(roi, '_roiIdx', None)
            if idxObj is not None:
                try:
                    pos = int(idxObj.get())
                except Exception:
                    continue
            else:
                # Only take residues from the requested chain prefix
                if not resName.startswith(prefix):
                    continue
                try:
                    pos = int(resName.split('_')[1])
                except Exception:
                    continue

            # Read frequency and percentage (pwobj -> python types)
            freqObj = getattr(roi, '_frequency', None)
            percObj = getattr(roi, '_percentage', None)
            freqVal = freqObj.get() if freqObj is not None else 0
            percVal = percObj.get() if percObj is not None else 0.0

            res_freq[pos] = int(freqVal)
            res_perc[pos] = float(percVal)

        # ---------------- Build output colormap set ----------------
        outSet = SetOfROIColorMap(filename=self._getPath('ColorMap.sqlite'))
        max_color = 100.0  # percentages are expected in [0, 100]

        for i, aa in enumerate(sequence, start=1):
            freq = res_freq.get(i, 0)
            perc = res_perc.get(i, 0.0)

            # Normalize to [0, 1] safely
            intensity = max(0.0, min(1.0, (perc / max_color) if max_color else 0.0))
            color = self._mapColor(intensity)

            item = ROIColorItem()
            item._residue = pwobj.String(f"{prefix}{i}")
            item._frequency = pwobj.Integer(freq)
            item._percentage = pwobj.Float(round(perc, 2))
            item._color = pwobj.String(color)
            outSet.append(item)

        # ---------------- Register output ----------------
        if len(outSet) > 0:
            outSet.setStore(True)
            outSet.write()
            self._defineOutputs(outputColorMap=outSet)
            self._defineSourceRelation(self.inputVotes, self.outputColorMap)
            self.info(f"Color map created for {len(sequence)} residues.")
        else:
            self.warning("No residues were mapped to the sequence (output set is empty).")

    # ---------------- COLOR UTILITY ----------------
    def _mapColor(self, intensity):
        """
        Map a normalized value in [0,1] to an RGB color string.
        Blue (low) → Red (high).
        """
        # Clamp just in case
        intensity = max(0.0, min(1.0, float(intensity)))
        r = int(255 * intensity)
        g = int(255 * (1 - intensity))
        b = 255 - r
        return f'#{r:02x}{g:02x}{b:02x}'

    # ---------------- INFO ----------------
    def _summary(self):
        seq = self.proteinSequence.get()
        n = len(seq) if seq else 0
        return [f"Generated color map for {n} residues.",
                "Each residue color intensity (blue→red) reflects its voting percentage.",
                "Frequencies are read directly from the ROI Voting results."]

    def _methods(self):
        return ["Each residue receives a color (blue to red) proportional to its voting percentage "
                "and frequency is preserved from the ROI Voting results."]