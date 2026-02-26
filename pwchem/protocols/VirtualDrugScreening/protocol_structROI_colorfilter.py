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
This protocol filters ROI voting results by frequency threshold.
Only residues with frequency >= threshold are kept.
The output contains only residue, frequency, and percentage columns.
"""

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwem.objects import EMSet, EMObject

class ROIFilterItem(EMObject):
    """Single filtered residue with frequency and percentage."""
    _possibleAttributes = ['_residue', '_frequency', '_percentage']

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residue = pwobj.String()
        self._frequency = pwobj.Integer()
        self._percentage = pwobj.Float()


class SetOfROIFiltered(EMSet):
    """Clean set with only 3 columns (residue, frequency, percentage)."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._itemType = ROIFilterItem


class ProtROIFrequencyFilter(EMProtocol):
    _label = 'ROI Frequency Filter (clean output)'

    # -------- Params --------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam,
                      pointerClass='EMSet',
                      label='Input ROI Voting Results',
                      help='Select the ROI voting output to filter.')
        form.addParam('minFrequency', params.IntParam, default=10,
                      label='Minimum frequency threshold',
                      help='Only residues with frequency ≥ this value are kept.')

    # -------- Steps ---------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')

    # -------- Main logic --------
    def filterStep(self):
        input_set = self.inputROIs.get()
        threshold = self.minFrequency.get()

        outSet = SetOfROIFiltered(filename=self._getPath('FilteredROIs.sqlite'))
        kept = 0

        for roi in input_set:
            # ---- frequency (pwobj.Integer -> int) ----
            freqObj = getattr(roi, '_frequency', None)
            freqVal = freqObj.get() if freqObj is not None else None
            if freqVal is None or freqVal < threshold:
                continue

            # ---- residue name can be _residue or _sequence ----
            resObj = getattr(roi, '_residue', None) or getattr(roi, '_sequence', None)
            resVal = resObj.get() if resObj is not None else ''

            # ---- percentage (pwobj.Float -> float) ----
            percObj = getattr(roi, '_percentage', None)
            percVal = percObj.get() if percObj is not None else 0.0

            new_roi = ROIFilterItem()
            new_roi._residue = pwobj.String(resVal)
            new_roi._frequency = pwobj.Integer(freqVal)
            new_roi._percentage = pwobj.Float(percVal)

            outSet.append(new_roi)
            kept += 1

        if kept > 0:
            outSet.setStore(True)
            outSet.write()
            self._defineOutputs(outputFilteredROIs=outSet)
            self._defineSourceRelation(self.inputROIs, self.outputFilteredROIs)
            self.info(f"Filter applied successfully — kept {kept} residues (freq ≥ {threshold}).")
        else:
            self.warning(f"No residues met the threshold ≥ {threshold}.")

    # -------- Summary --------
    def _summary(self):
        s = []
        if hasattr(self, 'outputFilteredROIs'):
            s.append(f"Residues filtered with frequency ≥ {self.minFrequency.get()}.")
        else:
            s.append("Protocol not executed or no residues met the threshold.")
        return s
