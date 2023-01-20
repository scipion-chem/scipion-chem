# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_EXPORT_CSV.PY
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import os, csv
from pyworkflow.tests import *

from pwchem.tests.tests_imports import TestImportBase
from pwchem.protocols import ProtChemExportCSV


class TestExportcsv(TestImportBase):
    @classmethod
    def _runExport(cls, inProt):
        protExport = cls.newProtocol(
            ProtChemExportCSV,
        )
        protExport.inputSet.set(inProt)
        protExport.inputSet.setExtended('outputSmallMolecules')

        cls.proj.launchProtocol(protExport, wait=True)
        return protExport

    def test(self):
        pExp = self._runExport(inProt=self.protImportSmallMols)

        csvFile = os.path.abspath(pExp._getPath('output.csv'))
        self.assertTrue(os.path.exists(csvFile), "CSV file was not create. Check if its location changed")

