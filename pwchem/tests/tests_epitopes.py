# ***************************************************************************
# *
# * Authors:     Daniel Del Hoyo (daniel.delhoyo.gomez@alumnos.upm.es)
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

# Scipion em imports
from pyworkflow.tests import DataSet, setupTestProject

# Scipion chem imports
from pwchem.tests import TestDefineSequenceROIs
from pwchem.protocols import ProtDefineSeqROI, ProtDefineMultiEpitope, ProtModifyMultiEpitope
from pwchem.utils import assertHandle

defSeqROIsSeq = '''1) Residues: {"index": "1-10", "residues": "MFVFLVLLPL", "desc": "None"}
2) Residues: {"index": "34-40", "residues": "RGVYYPD", "desc": "None"}
3) Residues: {"index": "56-64", "residues": "LPFFSNVTW", "desc": "None"}
4) Residues: {"index": "70-81", "residues": "VSGTNGTKRFDN", "desc": "None"}'''

defMultiEp = '''1) Epitope: SequenceROI (Idx: 34-40, ROI: RGVYYPD) (Set 103, Item 2, Name ROI_34-40)
2) Linker: LLLLLL
3) Epitope: SequenceROI (Idx: 70-81, ROI: VSGTNGTKRFDN) (Set 103, Item 4, Name ROI_70-81)
4) Linker: TTTTTTTT
5) Epitope: SequenceROI (Idx: 1-10, ROI: MFVFLVLLPL) (Set 103, Item 1, Name ROI_1-10)
'''

modMultiEp = '''1) Remove Epitope (ID 4): MFVFLVLLPL
2) Add Epitope : ASDCFERFDGTEF
3) Modify Linker (ID 3): TTTTTTTT to : GGGGGGG'''

class TestDefineMultiEpitope(TestDefineSequenceROIs):
	@classmethod
	def setUpClass(cls):
		setupTestProject(cls)
		cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')
		cls._runImportSequence()
		cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)


	@classmethod
	def _runDefSeqROIs(cls, inProt):
		protDefSeqROIs = cls.newProtocol(
			ProtDefineSeqROI,
			chooseInput=0, inROIs=defSeqROIsSeq
		)
		protDefSeqROIs.inputSequence.set(inProt)
		protDefSeqROIs.inputSequence.setExtended('outputSequence')

		cls.proj.launchProtocol(protDefSeqROIs, wait=False)
		return protDefSeqROIs

	@classmethod
	def _runDefMultiEp(cls, protROIs):
		protDefMultiEp = cls.newProtocol(
			ProtDefineMultiEpitope,
			multiSummary=defMultiEp
		)
		protDefMultiEp.inputROIsSets.set([protROIs])
		protDefMultiEp.inputROIsSets[0].setExtended('outputROIs')

		cls.proj.launchProtocol(protDefMultiEp, wait=False)
		return protDefMultiEp

	def test(self):
		protsROIs = self._runDefSeqROIs(inProt=self.protImportSeq)
		self._waitOutput(protsROIs, 'outputROIs', sleepTime=5)

		protMultiEp = self._runDefMultiEp(protsROIs)
		self._waitOutput(protMultiEp, 'outputROIs', sleepTime=5)
		assertHandle(self.assertIsNotNone, getattr(protMultiEp, 'outputROIs', None),
								 cwd=protMultiEp.getWorkingDir())

class TestModifyMultiEpitope(TestDefineMultiEpitope):
	@classmethod
	def _runModMultiEp(cls, protMultiEp):
		protModMultiEp = cls.newProtocol(
			ProtModifyMultiEpitope,
			modSummary=modMultiEp
		)
		protModMultiEp.inputMultiEpitope.set(protMultiEp)
		protModMultiEp.inputMultiEpitope.setExtended('outputROIs')

		cls.proj.launchProtocol(protModMultiEp, wait=False)
		return protModMultiEp

	def test(self):
		protsROIs = self._runDefSeqROIs(inProt=self.protImportSeq)
		self._waitOutput(protsROIs, 'outputROIs', sleepTime=5)

		protMultiEp = self._runDefMultiEp(protsROIs)
		self._waitOutput(protMultiEp, 'outputROIs', sleepTime=5)

		protMod = self._runModMultiEp(protMultiEp)
		self._waitOutput(protMod, 'outputROIs', sleepTime=5)
		assertHandle(self.assertIsNotNone, getattr(protMod, 'outputROIs', None),
								 cwd=protMod.getWorkingDir())

