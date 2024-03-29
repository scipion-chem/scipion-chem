# **************************************************************************
# *
# * Name:     TEST OF PROTOCOL_IMPORT_SMALLMOLECULES.PY
# *
# * Authors:    Alberto M. Parra Pérez (amparraperez@gmail.com)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

# General imports
import os

# Scipion em imports
from pyworkflow.tests import BaseTest, DataSet, setupTestProject

# Scipion chem imports
from pwchem.protocols import *
from pwchem.constants import *
from pwchem.utils import assertHandle

# Global variables
iedb_csv_str = '''Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Epitope,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object,Related Object
Epitope ID,Object Type,Description,Epitope Modified Residue(s),Epitope Modification(s),Starting Position,Ending Position,Non-peptidic epitope Accession,Epitope Synonyms,Antigen Name,Antigen Accession,Parent Protein,Parent Protein Accession,Organism Name,Parent Organism,Parent Organism ID,Epitope Comments,Epitope Relationship,Object Type,Description,Starting Position,Ending Position,Non-peptidic object Accession,Synonyms,Antigen Name,Parent Protein,Organism Name,Parent Organism
"234","Linear peptide","AAISDYDYY","","","4840","4848","","","orf1ab polyprotein [Severe acute respiratory syndrome coronavirus 2]","YP_009724389.1","Replicase polyprotein 1ab","P0DTD1","Severe acute respiratory syndrome coronavirus 2","Severe acute respiratory syndrome coronavirus 2","2697049","","","","","","","","","","","",""
"956","Linear peptide","AEGSRGGSQA","","","173","182","","","nucleocapsid phosphoprotein [Severe acute respiratory syndrome coronavirus 2]","YP_009724397.2","Nucleoprotein","P0DTC9","Severe acute respiratory syndrome coronavirus 2","Severe acute respiratory syndrome coronavirus 2","2697049","","","","","","","","","","","",""
"1220","Linear peptide","AEVQIDRLI","","","989","997","","","surface glycoprotein [Severe acute respiratory syndrome coronavirus 2]","YP_009724390.1","Spike glycoprotein","P0DTC2","Severe acute respiratory syndrome coronavirus 2","Severe acute respiratory syndrome coronavirus 2","2697049","","","","","","","","","","","",""'''

importErrorMessage = "There was a problem with the import"
emptySetOfMoleculesErrorMessage = "There was a problem with the import and the SetOfSmallMolecules is empty"

class TestImportBase(BaseTest):
	@classmethod
	def setUpClass(cls):
		cls.dsLig = DataSet.getDataSet("smallMolecules")
		setupTestProject(cls)

		cls._runImportSmallMols()
		cls._waitOutput(cls.protImportSmallMols, 'outputSmallMolecules', sleepTime=5)

	@classmethod
	def _runImportSmallMols(cls):
		cls.protImportSmallMols = cls.newProtocol(
			ProtChemImportSmallMolecules,
			filesPath=cls.dsLig.getFile('mol2'))
		cls.proj.launchProtocol(cls.protImportSmallMols, wait=False)

class TestImportSmallMolecules(TestImportBase):
	MULTIPLE_T = True
	MULTIPLE_F = False
	filesPattern = '*'

	def testImport_one_mol_smi(self):
		""" Import a single file of a small molecule provided by the user in smi format. """
		print("\nImport Experiment: 1 molecule smi format")

		kwargs = {
			'multipleFiles': self.MULTIPLE_F,
			'filePath': self.dsLig.getFile(os.path.join("mix", "2000_smi.smi"))
		}

		prot1 = self.newProtocol(ProtChemImportSmallMolecules, **kwargs)
		self.launchProtocol(prot1)
		small1 = getattr(prot1, 'outputSmallMolecules', None)
		assertHandle(self.assertIsNotNone, small1, message=importErrorMessage, cwd=prot1.getWorkingDir())
		assertHandle(self.assertTrue, small1.getSize()==1,
								 message=emptySetOfMoleculesErrorMessage, cwd=prot1.getWorkingDir())

	def testImport_one_mol_pdb(self):
			""" Import a single file of a small molecule provided by the user in pdb format
			"""
			print("\nImport Experiment: 1 molecule pdb format")

			kwargs = {
				'multipleFiles': self.MULTIPLE_F,
				'filePath': self.dsLig.getFile(os.path.join("mix", "2000_pdb.pdb"))
			}

			prot1 = self.newProtocol(ProtChemImportSmallMolecules, **kwargs)
			self.launchProtocol(prot1)
			small1 = getattr(prot1, 'outputSmallMolecules', None)
			assertHandle(self.assertIsNotNone, small1, message=importErrorMessage, cwd=prot1.getWorkingDir())
			assertHandle(self.assertTrue, small1.getSize()==1,
								 message=emptySetOfMoleculesErrorMessage, cwd=prot1.getWorkingDir())

	def testImport_one_mol_sdf(self):
			""" Import a single file of a small molecule provided by the user in sdf format
			"""
			print("\nImport Experiment: 1 molecule sdf format")

			kwargs = {
				'multipleFiles': self.MULTIPLE_F,
				'filePath': self.dsLig.getFile(os.path.join("mix", "2000_sdf.sdf"))
			}

			prot1 = self.newProtocol(ProtChemImportSmallMolecules, **kwargs)
			self.launchProtocol(prot1)
			small1 = getattr(prot1, 'outputSmallMolecules', None)
			assertHandle(self.assertIsNotNone, small1, message=importErrorMessage, cwd=prot1.getWorkingDir())
			assertHandle(self.assertTrue, small1.getSize()==1,
								 message=emptySetOfMoleculesErrorMessage, cwd=prot1.getWorkingDir())

	def testImport_mols_mix(self):
			""" Import several files of a small molecule provided by the user in different formats. """
			print("\nImport Experiment: 4 molecules in several formats")

			kwargs = {
				'multipleFiles': self.MULTIPLE_T,
				'filesPath': self.dsLig.getFile(os.path.join("mix")),
				'filesPattern': self.filesPattern
			}

			prot1 = self.newProtocol(ProtChemImportSmallMolecules, **kwargs)
			self.launchProtocol(prot1)
			small1 = getattr(prot1, 'outputSmallMolecules', None)
			assertHandle(self.assertIsNotNone, small1, message=importErrorMessage, cwd=prot1.getWorkingDir())
			assertHandle(self.assertTrue, small1.getSize()==4,
								 message=emptySetOfMoleculesErrorMessage, cwd=prot1.getWorkingDir())

class TestImportSequences(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.ds = DataSet.getDataSet('model_building_tutorial')
		cls.pImpDBs = cls._runImportDBIds()
		cls._waitOutput(cls.pImpDBs, 'outputDatabaseIDs', sleepTime=5)

	@classmethod
	def getIdsFilePath(cls, inType=0, dbName='UniProt'):
		return os.path.abspath(cls.proj.getPath('inputIDs_{}_{}.txt'.format(inType, dbName)))

	@classmethod
	def _runImportDBIds(cls, inType=0, dbName='UniProt'):
		idsFile = cls.getIdsFilePath(inType, dbName)
		with open(idsFile, 'w') as f:
			f.write('\n'.join(DB_IDS[inType][dbName]))

		protImport = cls.newProtocol(
			ProtChemImportSetOfDatabaseIDs,
			filePath=idsFile, databaseName=dbName
		)
		cls.proj.launchProtocol(protImport, wait=False)
		return protImport

	@classmethod
	def _runImportSeqs(cls):
		protImportSeqs = cls.newProtocol(
			ProtChemImportSetOfSequences,
			multiple=False,
			filePath=cls.ds.getFile('Sequences/Several_sequences.fasta'))
		cls.launchProtocol(protImportSeqs)
		cls.protImportSeqs = protImportSeqs

	@classmethod
	def _runImportSeqsFromDB(cls, inProt):
		protImportSeqs = cls.newProtocol(
			ProtChemImportSetOfSequences,
			fromFile=False, database=0)

		protImportSeqs.inputListID.set(inProt)
		protImportSeqs.inputListID.setExtended('outputDatabaseIDs')

		cls.launchProtocol(protImportSeqs)
		cls.protImportSeqsDB = protImportSeqs

	def test(self):
		self._runImportSeqs()
		self._runImportSeqsFromDB(inProt=self.pImpDBs)

		assertHandle(self.assertIsNotNone, getattr(self.protImportSeqs, 'outputSequences', None), cwd=self.protImportSeqs.getWorkingDir())
		assertHandle(self.assertIsNotNone, getattr(self.protImportSeqsDB, 'outputSequences', None), cwd=self.protImportSeqsDB.getWorkingDir())

class TestImportVariants(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)

	@classmethod
	def _runImportVariants(cls):
		protImportVariants = cls.newProtocol(
			ProtChemImportVariants,
			fromID=True,
			inputUniProtKB='P0DTC2')
		cls.proj.launchProtocol(protImportVariants, wait=False)
		cls.protImportVariants = protImportVariants

	def test(self):
		self._runImportVariants()
		self._waitOutput(self.protImportVariants, 'outputVariants')
		assertHandle(self.assertIsNotNone, getattr(self.protImportVariants, 'outputVariants', None), cwd=self.protImportVariants.getWorkingDir())

class TestImportSeqROIs(BaseTest):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		inFile = cls.proj.getTmpPath('epitopes.csv')
		with open(inFile, 'w') as f:
				f.write(iedb_csv_str)

	@classmethod
	def _runImportSeqROIs(cls):
		protImportSeqROIs = cls.newProtocol(
			ProtImportSeqROI,
			inputFile=cls.proj.getTmpPath('epitopes.csv'))
		cls.launchProtocol(protImportSeqROIs)
		return protImportSeqROIs

	def test(self):
		protImportSeqROIs = self._runImportSeqROIs()
		assertHandle(self.assertIsNotNone, getattr(protImportSeqROIs, 'outputROIs_P0DTC2', None), cwd=protImportSeqROIs.getWorkingDir())
