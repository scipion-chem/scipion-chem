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
This protocol combines the results from several models and counts how many times each residue appears.
"""

import sqlite3
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pwem.objects import EMSet, EMObject


class ROIVote(EMObject):
    """Objeto individual: residuo + frecuencia + porcentaje."""
    _possibleAttributes = ['_residue', '_frequency', '_percentage']  # fuerza schema mínimo

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._residue = pwobj.String()
        self._frequency = pwobj.Integer()
        self._percentage = pwobj.Float()


class SetOfROIVotes(EMSet):
    """Set limpio sin columnas adicionales."""
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._itemType = ROIVote


class ProtROIVoting(EMProtocol):
    _label = 'ROI voting (clean schema)'

    # ----------------- PARAMS -----------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Structural ROIs (runs)',
                      pointerClass='ProtDefineStructROIs',
                      help='Select all Define Structural ROIs runs to combine for voting.')

    # ----------------- STEPS -----------------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeVotingStep')

    # ----------------- MAIN LOGIC -----------------
    def computeVotingStep(self):
        residue_counts = {}

        # Leer cada run
        for p in self.roisList:
            try:
                prot = p.get()
                db_path = prot.outputStructROIs.getFileName()
                conn = sqlite3.connect(db_path)
                cur = conn.cursor()

                # Detectar columna con residuos
                cols = [c[1] for c in cur.execute("PRAGMA table_info(Objects);")]
                if "c06" in cols:
                    cur.execute("SELECT c06 FROM Objects;")
                elif "_contactResidues" in cols:
                    cur.execute("SELECT _contactResidues FROM Objects;")
                else:
                    self.warning(f"No residue column in {db_path}")
                    conn.close()
                    continue

                # Contar residuos
                for (val,) in cur.fetchall():
                    if not val:
                        continue
                    for r in (x.strip() for x in val.split('-') if x.strip()):
                        if r.startswith("A_"):
                            residue_counts[r] = residue_counts.get(r, 0) + 1
                conn.close()
            except Exception as e:
                self.warning(f"Error reading {p}: {e}")

        if not residue_counts:
            self.warning("No residues found for chain A.")
            return

        # Normalizar a 100%
        max_count = max(residue_counts.values())
        self.top_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)

        # Crear set limpio
        outSet = SetOfROIVotes(filename=self._getPath('ROIVoting.sqlite'))
        for residue, count in self.top_residues:
            item = ROIVote()
            item._residue = pwobj.String(residue)
            item._frequency = pwobj.Integer(count)
            item._percentage = pwobj.Float(round((count / max_count) * 100.0, 2))
            outSet.append(item)

        # Registrar output
        if len(outSet) > 0:
            self._defineOutputs(outputROIVotes=outSet)
            self._defineSourceRelation(self.roisList, self.outputROIVotes)
            self.info(f"Voting results saved ({len(outSet)} residues).")
        else:
            self.warning("No output generated after voting.")

    def _methods(self):
        return ["Residues are scored by frequency across ROI sets; "
                "percentages normalized to 100% max residue."]
    

