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
from bioinformatics.protocols import ProtBioinformaticsImportSmallMolecules as PISmallM
from bioinformatics.protocols import ProtBioinformaticsExportCSV as PEcsv


class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet("ligandLibraries")


class TestImportData(TestImportBase):
    """ Import a set of small molecules, a set of database IDs and a set of binding sites (last 2 coming soon)
    """

    def _importSetMolecules(self, path, pattern):
        MULTIPLE_T = True
        args = {'multiple': MULTIPLE_T,
                'filesPath': path,
                'filesPattern': pattern
                }

        prot1 = self.newProtocol(PISmallM, **args)
        self.launchProtocol(prot1)
        small_1 = prot1.outputSmallMols
        return small_1

class TestExportcsv(TestImportData):

    def testCSVMolecules(self):

        print("\nCSV Export Experiment:  4 small molecules different formats")

        #Import Set of Small Molecules
        path = self.dsModBuild.getFile(os.path.join("mix"))
        pattern = "*"
        smallM = self._importSetMolecules(path,pattern)

        #Export CSV
        args = {'inputSet': smallM}
        pcsv = self.newProtocol(PEcsv, **args)
        self.launchProtocol(pcsv)

        #Search csv output in the directory tree of Runs
        for root, directories, files in os.walk(os.path.join(os.getcwd(), "Runs")):
            for name in files:
                if name=="output.csv":
                    pathcsv = os.path.join(root,name)

        #Check if the file exits, if is not none and if it has got a correct number of entries
        self.assertTrue(os.path.exists(pathcsv), "CSV file was not create. Check if its location changed")

        csvtest = open(pathcsv)
        self.assertIsNotNone(csvtest, "Created CSV is empty")

        reader = csv.reader(csvtest)
        lines = len(list(reader)) #Number of lines in the csv
        self.assertTrue(lines==(4+1), "Created CSV is incomplete (missing entries (rows))")
        #47 molecules and 1 header row

        csvtest.close()

