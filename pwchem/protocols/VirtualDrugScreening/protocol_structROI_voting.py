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

"""
This protocol combines the results from several ROI sets and counts how many
times each residue appears across models (one count per model, not per pocket).
It outputs a SetOfStructROIs and a SequenceChem with per-residue frequency attributes.
"""

import os
import sqlite3
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler
from pyworkflow.protocol import params
from pwchem.objects import SetOfStructROIs, StructROI, SequenceChem
from pwchem.utils import getBaseName


class ProtROIVoting(EMProtocol):
    _label = 'ROI voting (multi-model integration)'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Sets of Structural ROIs',
                      pointerClass='SetOfStructROIs',
                      help='Select all SetOfStructROIs to combine.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeVotingStep)

    def computeVotingStep(self):
        residueCounts = {}
        proteinFile = None

        for p in self.roisList:
            try:
                roiSet = p.get()
                if proteinFile is None:
                    proteinFile = roiSet.getProteinFile()

                dbPath = roiSet.getFileName()
                conn = sqlite3.connect(dbPath)
                cur = conn.cursor()

                cols = [c[1] for c in cur.execute("PRAGMA table_info(Objects);")]
                if "c06" in cols:
                    cur.execute("SELECT c06 FROM Objects;")
                elif "_contactResidues" in cols:
                    cur.execute("SELECT _contactResidues FROM Objects;")
                else:
                    self.warning(f"No residue column found in {dbPath}")
                    conn.close()
                    continue

                modelResidues = set()
                for (val,) in cur.fetchall():
                    if not val:
                        continue
                    for r in (x.strip() for x in val.split('-') if x.strip()):
                        if r.startswith("A_"):
                            modelResidues.add(r)
                conn.close()

                for r in modelResidues:
                    residueCounts[r] = residueCounts.get(r, 0) + 1

            except Exception as e:
                self.warning(f"Error reading {p}: {e}")

        if not residueCounts:
            self.warning("No residues found for chain A.")
            return

        maxCount = max(residueCounts.values())
        topResidues = sorted(residueCounts.items(), key=lambda x: x[1], reverse=True)

        outSet = SetOfStructROIs(filename=self._getPath('ROIVoting.sqlite'))

        for residue, count in topResidues:
            chain, resnum = residue.split('_')
            pocketCoords = []
            with open(proteinFile) as f:
                for line in f:
                    if line.startswith('ATOM'):
                        parts = line.split()
                        if parts[4] == chain and parts[5] == resnum:
                            pocketCoords.append(line)

            pocketFile = self._getExtraPath(f'roi_{residue}.pdb')
            with open(pocketFile, 'w') as f:
                for i, line in enumerate(pocketCoords):
                    parts = line.split()
                    x, y, z = float(parts[6]), float(parts[7]), float(parts[8])
                    f.write(f"HETATM{i:5d} APOL STP C   1    {x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          Ve\n")

            roi = StructROI(pocketFile, proteinFile=None)
            roi._proteinFile = pwobj.String(proteinFile)
            roi._contactResidues = pwobj.String(residue)
            roi._frequency = pwobj.Integer(count)
            roi._percentage = pwobj.Float(round((count / maxCount) * 100.0, 2))
            outSet.append(roi)

        if len(outSet) > 0:
            self._defineOutputs(outputStructROIs=outSet)
            self._defineSourceRelation(self.roisList, self.outputStructROIs)

        if proteinFile and os.path.exists(proteinFile):
            ash = AtomicStructHandler()
            ash.read(proteinFile)
            seqStr = str(ash.getSequenceFromChain(modelID=0, chainID='A'))
            nRes = len(seqStr)

            freqValues = [residueCounts.get(f'A_{i}', 0) for i in range(1, nRes + 1)]
            base = getBaseName(proteinFile)
            outSeq = SequenceChem(
                name=f'{base}_A',
                sequence=seqStr,
                id=f'{base}_A',
                attributesFile=self._getExtraPath('sequenceAttributes.txt')
            )
            outSeq.addAttributes({'frequency': freqValues})
            self._defineOutputs(outputSequence=outSeq)
        else:
            self.warning("Could not find protein file for sequence output.")

    def _summary(self):
        summary = []
        if hasattr(self, 'outputStructROIs'):
            summary.append(f"Voting results: {self.outputStructROIs.getSize()} residues.")
        else:
            summary.append("No residues detected or protocol not executed yet.")
        return summary

    def _methods(self):
        return ["Residues are scored by frequency across ROI runs (one vote per model), "
                "and normalized so that the highest count equals 100%."]