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

# *****************************************************************************
# * Scipion-CHEM - Compute Interaction Intensity Protocol
# *****************************************************************************

"""
This protocol combines the results from several models and counts how many times each residue appears.
"""

import sqlite3
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwchem.objects import SetOfStructROIs, StructROI


class ProtROIVoting(EMProtocol):
    _label = 'ROI voting (normalized)'

    # --------------------- PARAMS ---------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Structural ROIs (runs)',
                      pointerClass='ProtDefineStructROIs',
                      help='Select all Define Structural ROIs runs to combine for voting.')

    # --------------------- STEPS ----------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeVotingStep')

    # --------------------- MAIN LOGIC -----------------
    def computeVotingStep(self):
        """
        Reads each input StructROIs.sqlite file, aggregates residue frequencies (only chain A),
        and saves a SetOfStructROIs with:
            _residue (String)
            _frequency (Integer)
            _percentage (Float)
        The percentage is normalized so that the residue with the highest frequency = 100%.
        """
        residue_counts = {}

        # Iterate through all selected DefineStructuralROIs runs
        for p in self.roisList:
            try:
                prot = p.get()
                db_path = prot.outputStructROIs.getFileName()
                conn = sqlite3.connect(db_path)
                cur = conn.cursor()

                # Identify the correct column for residues
                cols = [c[1] for c in cur.execute("PRAGMA table_info(Objects);")]
                if "c06" in cols:
                    cur.execute("SELECT c06 FROM Objects;")
                elif "_contactResidues" in cols:
                    cur.execute("SELECT _contactResidues FROM Objects;")
                else:
                    self.warning(f"No residue column found in {db_path}")
                    conn.close()
                    continue

                # Count each residue occurrence
                for (val,) in cur.fetchall():
                    if not val:
                        continue
                    for r in (x.strip() for x in val.split('-') if x.strip()):
                        if r.startswith("A_"):  # keep only chain A residues
                            residue_counts[r] = residue_counts.get(r, 0) + 1

                conn.close()
            except Exception as e:
                self.warning(f"Error reading {p}: {e}")

        if not residue_counts:
            self.warning("No residues found for chain A in any ROI input.")
            return

        # Normalize frequencies to the maximum
        max_count = max(residue_counts.values())
        self.top_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)

        # Create output SetOfStructROIs
        outSet = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for residue, count in self.top_residues:
            roi = StructROI()
            roi._residue = pwobj.String(residue)
            roi._frequency = pwobj.Integer(int(count))
            roi._percentage = pwobj.Float(round((count / max_count) * 100.0, 2))
            outSet.append(roi)

        # Register output in Scipion environment
        if len(outSet) > 0:
            self._defineOutputs(outputStructROIs=outSet)
            self._defineSourceRelation(self.roisList, self.outputStructROIs)
            self.info(f"Voting results saved: {len(outSet)} residues processed "
                      f"(max votes = {max_count}).")
        else:
            self.warning("No ROIs were generated after voting.")

    # --------------------- METHODS INFO ---------------------
    def _methods(self):
        return ["Residues are scored by frequency across input ROI sets, "
                "with percentages normalized to the highest count (100%)."]
