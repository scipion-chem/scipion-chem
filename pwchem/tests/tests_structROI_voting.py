# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
# *              Irene Sánchez Martín (100495638@alumnos.uc3m.es)
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
# ***************************************************************************

from pwchem.protocols import ProtROIVoting
from pwchem.utils import assertHandle

from .tests_structROIs import TestDefineStructROIs

defROIsStr = '''1) Coordinate: {"X": 45, "Y": 65, "Z": 60}
2) Residues: {"model": 0, "chain": "A", "index": "1-4", "residues": "VLSP"}'''

defROIsStr2 = '''1) Residues: {"model": 0, "chain": "A", "index": "61-63", "residues": "KVA"}
2) Residues: {"model": 0, "chain": "A", "index": "80-84", "residues": "LSALS"}'''


class TestROIVoting(TestDefineStructROIs):
    @classmethod
    def _runROIVoting(cls, protDefs):
        prot = cls.newProtocol(ProtROIVoting)
        for protDef in protDefs:
            prot.roisList.append(protDef)
            prot.roisList[-1].setExtended('outputStructROIs')
        cls.proj.launchProtocol(prot)
        return prot

    def test(self):
        pDef1 = self._runDefStructROIs(defROIsStr)
        pDef2 = self._runDefStructROIs(defROIsStr2)
        self._waitOutput(pDef1, 'outputStructROIs', sleepTime=10)
        self._waitOutput(pDef2, 'outputStructROIs', sleepTime=10)

        pVoting = self._runROIVoting([pDef1, pDef2])
        self._waitOutput(pVoting, 'outputStructROIs', sleepTime=10)
        assertHandle(self.assertIsNotNone, getattr(pVoting, 'outputStructROIs', None),
                     cwd=pVoting.getWorkingDir())