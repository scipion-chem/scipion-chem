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
from pwchem.objects import SetOfStructROIs, StructROI


class ProtROIVoting(EMProtocol):
    _label = 'ROI voting (multi-run integration)'

    # --------- Params ---------
    def _defineParams(self, form):
        form.addSection(label='Input')
        # IMPORTANT: seleccionas los *protocols* DefineStructuralROIs,
        # porque de ahí leemos outputStructROIs.getFileName()
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Structural ROIs (runs)',
                      pointerClass='ProtDefineStructROIs',
                      help='Selecciona todos los runs de "Define structural ROIs" que quieres combinar.')

    # --------- Steps ----------
    def _insertAllSteps(self):
        self._insertFunctionStep('computeVotingStep')

    def computeVotingStep(self):
        """
        Lee cada StructROIs.sqlite de los runs de entrada,
        acumula frecuencia por residuo (A_###, B_###, ...),
        y guarda un SetOfStructROIs con _residue, _frequency, _percentage.
        """
        residue_counts = {}
        outSet = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))

        # Recorrer todos los runs seleccionados
        for p in self.roisList:
            try:
                prot = p.get()
                db_path = prot.outputStructROIs.getFileName()  # <- de cada run cogemos su sqlite
                conn = sqlite3.connect(db_path)
                cur = conn.cursor()

                # Detectar columna con residuos
                cols = [c[1] for c in cur.execute("PRAGMA table_info(Objects);")]
                if "c06" in cols:
                    cur.execute("SELECT c06 FROM Objects;")
                elif "_contactResidues" in cols:
                    cur.execute("SELECT _contactResidues FROM Objects;")
                else:
                    self.warning(f"No contact residue column in {db_path}")
                    conn.close()
                    continue

                for (val,) in cur.fetchall():
                    if not val:
                        continue
                    for r in [x.strip() for x in val.split('-') if x.strip()]:
                        residue_counts[r] = residue_counts.get(r, 0) + 1

                conn.close()
            except Exception as e:
                self.warning(f"Error reading run {p}: {e}")

        if not residue_counts:
            self.warning("No residues found in any ROI input.")
            return

        total_models = max(len(self.roisList), 1)
        self.top_residues = sorted(residue_counts.items(), key=lambda x: x[1], reverse=True)

        # Añadir items con atributos pwobj.* para que se persistan como columnas
        for residue, count in self.top_residues:
            # Solo incluir residuos de la cadena A
            if not residue.startswith("A_"):
                continue

            roi = StructROI()
            roi._residue    = pwobj.String(residue)
            roi._frequency  = pwobj.Integer(int(count))
            roi._percentage = pwobj.Float(round((count / total_models) * 100.0, 2))
            outSet.append(roi)  # <- este sí existe en SetOfStructROIs

        if len(outSet) > 0:
            self._defineOutputs(outputStructROIs=outSet)
            self._defineSourceRelation(self.roisList, self.outputStructROIs)
            self.info(f"Voting results saved: {len(outSet)} residues total.")
        else:
            self.warning("No ROIs were generated after voting.")

    # --------- UI info ---------
    def _summary(self):
        s = []
        if hasattr(self, 'top_residues') and self.top_residues:
            s.append("Top-5 most frequent residues:")
            for r, n in self.top_residues[:5]:
                s.append(f"  {r}: {n} votes")
        else:
            s.append("No residues detected or protocol not executed yet.")
        return s

    def _methods(self):
        return ["Residues are scored by occurrence across input ROI sets."]
