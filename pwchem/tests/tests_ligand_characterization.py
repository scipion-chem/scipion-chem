# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
# *				 Irene Sánchez Martín (100495638@alumnos.uc3m.es)
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
# ***************************************************************************

import os
import unittest

from pwchem.protocols.VirtualDrugScreening.protocol_ligand_characterization import ProtocolLigandCharacterization
from pwchem.tests.tests_imports import TestImportBase
from pwchem.utils import assertHandle
import pwchem

# RDKit y script externo 
try:
    from rdkit import Chem  
    RDKit_OK = True
except Exception:
    RDKit_OK = False

SCRIPTS_DIR = os.path.join(os.path.dirname(pwchem.__file__), "scripts")
DESCR_SCRIPT = os.path.join(SCRIPTS_DIR, "ligand_descriptor_calc.py")
SCRIPT_EXISTS = os.path.exists(DESCR_SCRIPT)


class TestLigandCharacterization(TestImportBase):
    """Tests for ProtocolLigandCharacterization."""

    #Test 1: Always run
    def test_protocol_instantiates(self):
        """El protocolo se instancia y expone sus parámetros sin crashear."""
        prot = self.newProtocol(ProtocolLigandCharacterization)
        self.assertEqual(prot.getLabel(), "Ligand characterization")
        form = prot.getForm()
        self.assertIsNotNone(form)

    # To launch protocol
    @classmethod
    def _runLigandCharacterization(cls, inProt):
        protChar = cls.newProtocol(ProtocolLigandCharacterization)

        protChar.inputSmallMolecules.set(inProt)
        protChar.inputSmallMolecules.setExtended("outputSmallMolecules")
        cls.proj.launchProtocol(protChar, wait=False)
        return protChar
    
    # Only if RDKIT + script
    @unittest.skipIf(not RDKit_OK, "RDKit no disponible en el entorno de tests")
    @unittest.skipIf(not SCRIPT_EXISTS, "Falta pwchem/scripts/ligand_descriptor_calc.py")
    def test_ligand_characterization_runs(self):
        """Ejecuta el protocolo y verifica que genera outputSmallMolecules con tamaño > 0."""
        protChar = self._runLigandCharacterization(inProt=self.protImportSmallMols)

        # 0) Wait
        self._waitOutput(protChar, "outputSmallMolecules", sleepTime=10)

        # 1) It exists
        out = getattr(protChar, "outputSmallMolecules", None)
        assertHandle(self.assertIsNotNone, out, cwd=protChar.getWorkingDir())

        # 2) Has elements
        assertHandle(self.assertGreater, out.getSize(), 0, cwd=protChar.getWorkingDir())

