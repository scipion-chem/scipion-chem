# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_LIST_OPERATE.PY
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

import os
from pathlib import Path
from pyworkflow.tests import *
from pyworkflow.protocol import *
from pwem.protocols.protocol_import import ProtImportPdb
from pwchem.protocols import ProtChemListOperate as LOperate
from pwchem.protocols import ProtChemDali as DALI


class TestImportBase(BaseTest):

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet("ligandLibraries")


class TestCreateOutput(TestImportBase):
    """ Create a realistic output of Dali using an old job (pdb id used:6zow)
    """

    def _importPDB(self, path):
        inputPdbData = 1 #file
        args = {'inputPdbData': inputPdbData,
                'pdbFile': path
                }

        prot1 = self.newProtocol(ProtImportPdb, **args)
        self.launchProtocol(prot1)
        global pdb
        pdb = prot1.outputPdb
        return pdb

    def getDalioutputs(self, protocol):
        fnBaseDir = self.dsModBuild.getFile(os.path.join("dali_output"))
        for fn in Path(fnBaseDir).rglob('*.txt'):
            protocol.constructOutput(str(fn), protocol)



class TestListOperate(TestCreateOutput):

    def test_1filter(self):
        """1. Filter a column using different ways
        """

        print("\n Filtering a SetOfDatabaseID object by a column")

        # Import a pdb file
        path = self.dsModBuild.getFile(os.path.join("proteinpdb", "5xjh.pdb"))
        pdb = self._importPDB(path)

        #Dali
        args = {'inputStructure': pdb,
                'method': 0, #search
                'title': 'test_dali_5xjh',
                'email': 'amparraperez@gmail.com'}

        prot1 = self.newProtocol(DALI, **args)
        prot1._store()
        self.getDalioutputs(prot1)
        prot1.setStatus(STATUS_FINISHED)

        global outputDali; global outputDali2
        outputDali = prot1.outputDatabaseIds90
        outputDali2 = prot1.outputDatabaseIds50

        self.assertIsNotNone(outputDali, "Error in creation of SetOfDatabaseID by DALI - It is NONE")
        self.assertTrue(outputDali.getSize() == 432, "The size of the SetOfDatabaseID is wrong. It should be 432 ")
        n_columns = len(list(outputDali.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 11)  # 2 fixed columns + 9 given


        # SET filtering :  _DaliZscore column ( >= 40)
        args = {'operation': 0,  # Filter
                'inputSet':outputDali,
                'filterColumn': '_DaliZscore',
                'filterOp': 2,  # >=
                'filterValue': 40}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 7, "Error in filtering >= 40 of _DaliZscore column")


        # SET filtering :  _DaliZscore column ( 25 >= value <= 52 )
        args = {'operation': 0,  # Filter
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore',
                'filterOp': 6,  # between
                'filterValue': 52,
                'filterValue2': 25}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize()==8, "Error in filtering 25 >= value <= 52 of _DaliZscore column")

        # SET filtering :  _DaliZscore column ( >= 40)
        args = {'operation': 0,  # Filter
                'inputSet': outputDali,
                'filterColumn': '_pdbId',
                'filterOp': 7,  # startwith
                'filterValue': "6"}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 57, "Error in filtering >_pdbId column. Startwith does not work")

        # SET filtering :  _DaliZscore column ( >= 40)
        args = {'operation': 0,  # Filter
                'inputSet': outputDali,
                'filterColumn': '_DaliDescription',
                'filterOp': 9,  # contains
                'filterValue': "ESTERASE"}

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new filtered SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 104, "Error in filtering _DaliDescription column searching the ESTERASE word")


    def test_2keepcolumns(self):
        """2. Keep only 2 columns and all entries
        """
        print("\n Keeping 2 interesting columns in a new SetOfDatabaseID")

        # Keep_column :  _pdbId column
        args = {'operation': 1,
                'inputSet': outputDali,
                'keepColumns': '_pdbId ; _DaliZscore'
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output


        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432, "Error in creation of a new SetOfDatabaseID. It is different from the original regarding the number of entries")
        n_columns = len(list(setf.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 2+2) #2 fixed columns + 2 given


    def test_3unique(self):
        """3. Keep unique entries in 1 column
        """
        print("\n Keeping unique entries in 1 column in a new SetOfDatabaseID. No biological sense")


        args = {'operation': 2,  # Unique
                'inputSet': outputDali,
                'filterColumn': '_DaliDescription'
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 302, "Error in creation of a new SetOfDatabaseID. There is not unique entries")


    def test_4top(self):
        """4. Keep 10 top and 10 top 5% entries regarding 1 column
        """
        print("\n Keeping a number (10 or 5%) of entries regarding if it is in the top  regarding 1 column")

        # Keep 10 top entries :
        args = {'operation': 3,
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore',
                'N': 10
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output


        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 10, "Error in creation of a new SetOfDatabaseID. The number of entries is not the indicated one")

        for entry in setf:
            first_value = entry.getAttributeValue('_DaliZscore')
            break

        self.assertTrue(first_value == 51.4, "Failed to extract all top 10 entries regarding _DaliZscorecolumn")

        # Keep 5 % top entries :
        args = {'operation': 5,
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore',
                'N': 10
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 23, "Error in creation of a new SetOfDatabaseID. There is not unique entries")

        for entry in setf:
            first_value = entry.getAttributeValue('_DaliZscore')
            break

        self.assertTrue(first_value == 51.4, "Failed to extract all 5% top 10 entries regarding _DaliZscorecolumn")


    def test_5count(self):
        """5. Count the number of same entries
        """
        print("\n Count the number of same entries regarding 1 column")


        args = {'operation': 7,
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore'
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432, "Error in creation of a new SetOfDatabaseID.")
        n_columns = len(list(setf.getFirstItem().getAttributes()))
        self.assertTrue(n_columns == 2+9+1, "Column count was not created") #2 fixed columns + 9 given + 1 count


    def test_6intersect(self):
        """6. Intersect 2 SetOfDatabaseID using the column called _pdbId
        """
        print("\n Intersect 2 SetOfDatabaseID using 1 column")

        args = {'operation': 8,
                'inputSet': outputDali,
                'secondSet': outputDali2,
                'filterColumn': '_pdbId' # Without chain information. With that info there is only 1 difference (311 in intersection 312 in 50%)
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 314, "Error in creation of a new SetOfDatabaseID. The intersection was wrong")


    def test_7sort(self):
        """6. Sort (Ascending or Descending way) a SetOfDatabaseID regarding 1 column (_DaliZscore)
        """
        print("\n Sort a SetOfDatabaseID regarding 1 column (Ascending or Descending way)")

        args = {'operation': 9,
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore',# Without chain information. With that info there is only 1 difference (311 in intersection 312 in 50%)
                'direction': 0 #ascending
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432, "Error in creation of a new SetOfDatabaseID")

        for entry in setf:
            first_value = entry.getAttributeValue('_DaliZscore')
            break

        self.assertTrue(first_value == 2.1, "Failed to sort (ascending) the SetDatabaseID regarding _DaliZscorecolumn")


        args = {'operation': 9,
                'inputSet': outputDali,
                'filterColumn': '_DaliZscore',# Without chain information. With that info there is only 1 difference (311 in intersection 312 in 50%)
                'direction': 1 #descending
                }

        setf = self.newProtocol(LOperate, **args)
        self.launchProtocol(setf)
        setf = setf.output

        self.assertIsNotNone(setf, "Error in creation of a new SetOfDatabaseID - It is NONE")
        self.assertTrue(setf.getSize() == 432, "Error in creation of a new SetOfDatabaseID")

        for entry in setf:
            first_value = entry.getAttributeValue('_DaliZscore')
            break

        self.assertTrue(first_value == 51.4, "Failed to sort (descending) the SetDatabaseID regarding _DaliZscorecolumn")

