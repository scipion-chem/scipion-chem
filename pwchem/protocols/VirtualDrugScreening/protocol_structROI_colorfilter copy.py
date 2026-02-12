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
This protocol applies a frequency threshold to the ROI voting output.
Residues with frequency >= threshold are labeled 'YES'; otherwise 'NO'.
"""

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwchem.objects import SetOfStructROIs, StructROI


class ProtROIFrequencyFilter(EMProtocol):
    _label = 'ROI Frequency Filter (threshold YES/NO)'

    # ---------------- PARAMETERS ----------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputROIs', params.PointerParam,
                      pointerClass='SetOfROIVotes',
                      label='Input ROI Voting Results',
                      help='Select the ROI voting output to filter.')
        form.addParam('minFrequency', params.IntParam, default=10,
                      label='Minimum frequency threshold',
                      help='Only residues with frequency >= this value are kept.')

    # ---------------- EXECUTION -----------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')

    def filterStep(self):
        input_set = self.inputROIs.get()
        output_set = SetOfStructROIs(filename=self._getPath('FilteredROIs.sqlite'))

        threshold = self.minFrequency.get()
        kept = 0

        for roi in input_set:
            freq = getattr(roi, '_frequency', None)
            if freq is not None and freq >= threshold:
                new_roi = StructROI()
                new_roi._residue = pwobj.String(getattr(roi, '_residue', ''))
                new_roi._frequency = pwobj.Integer(freq)
                new_roi._percentage = pwobj.Float(getattr(roi, '_percentage', 0.0))
                output_set.append(new_roi)
                kept += 1

        if kept > 0:
            self._defineOutputs(outputFilteredROIs=output_set)
            self._defineSourceRelation(self.inputROIs, self.outputFilteredROIs)
            self.info(f"Filter applied successfully — kept {kept} residues (freq ≥ {threshold}).")
        else:
            self.warning(f"No residues met the threshold ≥ {threshold}.")

    # ---------------- DISPLAY -------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'outputFilteredROIs'):
            summary.append(f"Residues filtered with frequency ≥ {self.minFrequency.get()}.")
        else:
            summary.append("Protocol not executed or no residues met the threshold.")
        return summary
