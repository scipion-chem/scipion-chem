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
This protocol combines the results from several ROI models and counts how many times each residue appears.
It outputs a clean Scipion Set (SetOfROIVotes) visible in the environment.
"""



import sqlite3
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwchem.objects.base import ROISequence   # importamos la clase base


class ProtROIVoting(EMProtocol):
    _label = 'ROI voting (multi-model integration)'

    # ------------------ PARAMETERS ------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Structural ROI runs',
                      pointerClass='ProtDefineStructROIs',
                      help='Select all Define Structural ROIs runs to combine.')

    # ------------------ STEPS ------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeVotingStep')

    # ------------------ MAIN LOGIC ------------------
    def computeVotingStep(self):
        residue_counts = {}

        # Iterate through all input runs
        for p in self.roisList:
            try:
                prot = p.get()
                db_path = prot.outputStructROIs.getFileName()
                conn = sqlite3.connect(db_path)
                cur = conn.cursor()

                # Detect residue column
                cols = [c[1] for c in cur.execute("PRAGMA table_info(Objects);")]
                if "c06" in cols:
                    cur.execute("SELECT c06 FROM Objects;")
                elif "_contactResidues" in cols:
                    cur.execute("SELECT _contactResidues FROM Objects;")
                else:
                    self.warning(f"No residue column found in {db_path}")
                    conn.close()
                    continue

                # Count residues
                for (val,) in cur.fetchall():
                    if not val:
                        continue
                    for r in (x.strip() for x in val.split('-') if x.strip()):
                        if r.startswith("A_"):  # only chain A
                            residue_counts[r] = residue_counts.get(r, 0) + 1
                conn.close()

            except Exception as e:
                self.warning(f"Error reading {p}: {e}")

        if not residue_counts:
            self.warning("No residues found for chain A.")
            return

        # Normalize frequencies so max = 100%
        max_count = max(residue_counts.values())
        self.top_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)

        # ✅ Create unified Scipion Set (using ROISequence)
        outSet = self._createSetOfBaseObjects(self._getPath('ROIVoting.sqlite'))

        for residue, count in self.top_residues:
            item = ROISequence()
            item._sequence = pwobj.String(residue)
            item._roiIdx = pwobj.Integer(int(residue.split('_')[1]) if '_' in residue else None)
            item._frequency = pwobj.Integer(count)
            item._percentage = pwobj.Float(round((count / max_count) * 100.0, 2))
            outSet.append(item)

        # Register output in Scipion (official blue output)
        if len(outSet) > 0:
            self._defineOutputs(outputSet=outSet)
            self._defineSourceRelation(self.roisList, self.outputSet)
            self.info(f"Voting results saved ({len(outSet)} residues).")
        else:
            self.warning("No output generated after voting.")

    # ------------------ SUMMARY ------------------
    def _summary(self):
        summary = []
        if hasattr(self, 'top_residues') and self.top_residues:
            summary.append("Top-5 most frequent residues:")
            for r, n in self.top_residues[:5]:
                summary.append(f"  {r}: {n} votes")
        else:
            summary.append("No residues detected or protocol not executed yet.")
        return summary

    def _methods(self):
        return ["Residues are scored by frequency across ROI runs, "
                "and normalized so that the highest count equals 100%."]

    # ------------------ INTERNAL HELPERS ------------------
    def _createSetOfBaseObjects(self, filename):
        """Creates a new Set based on the ROISequence base class."""
        from pwem.objects import EMSet
        outSet = EMSet(filename=filename)
        outSet._itemType = ROISequence
        return outSet
