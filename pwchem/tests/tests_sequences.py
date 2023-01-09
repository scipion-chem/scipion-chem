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
# ***************************************************************************

from pyworkflow.tests import BaseTest, DataSet, setupTestProject
from pwem.protocols import ProtImportPdb, ProtImportSequence
from pwem.convert.atom_struct import AtomicStructHandler
from pwem.objects import AtomStruct, Sequence, Pointer

from ..utils import getBaseFileName
from ..tests import TestImportVariants, TestImportSequences
from ..protocols import *


CLUSTALO, MUSCLE, MAFFT = 0, 1, 2

defSeqROIsVar = '''1) Residues: {"index": "1-3", "residues": "MFV", "desc": "Residues"}
2) Variant: Alpha
3) Mutations: L18F'''

defSeqROIsSeq = '''1) Residues: {"index": "1-5", "residues": "MFVFL", "desc": "Residues"}
2) Residues: {"index": "42-45", "residues": "VFRS", "desc": "Residues2"}'''

genSeqsStr = '''1) Variant: Original
2) Variant: Alpha
3) Mutations: T478K
'''

defSetASChain, defSetPDBChain = 'A', 'B'
defSetSeqFile = 'Tmp/P0DTC2_FIRST-LAST.fa'
defSetASFile = 'Tmp/4erf_{}_FIRST-LAST.fa'.format(defSetASChain)
defSetPDBFile = 'Tmp/5ni1_{}_FIRST-LAST.fa'.format(defSetPDBChain)

names = ['P0DTC2', '4erf', '5ni1']
defSetChains = [None, defSetASChain, defSetPDBChain]
defSetFiles = [defSetSeqFile, defSetASFile, defSetPDBFile]

defSetSeqs = '''1) {"name": "%s", "index": "FIRST-LAST", "seqFile": "%s"}
2) {"name": "%s", "chain": "%s", "index": "FIRST-LAST", "seqFile": "%s"}
3) {"name": "%s", "chain": "%s", "index": "FIRST-LAST", "seqFile": "%s"}''' % \
             (names[0], defSetSeqFile, names[1], defSetASChain, defSetASFile, names[2], defSetPDBChain, defSetPDBFile)

class TestDefineSequenceROIs(TestImportVariants):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')

        cls._runImportSequence()
        cls._runImportVariants()
        cls._waitOutput(cls.protImportVariants, 'outputVariants', sleepTime=5)
        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)

        cls.inProts = [cls.protImportSeq, cls.protImportVariants]

    @classmethod
    def _runImportSequence(cls):
        args = {'inputSequenceName': 'User_seq',
                'inputProteinSequence': ProtImportSequence.IMPORT_FROM_UNIPROT,
                'uniProtSequence': 'P0DTC2'
                }
        cls.protImportSeq = cls.newProtocol(ProtImportSequence, **args)
        cls.proj.launchProtocol(cls.protImportSeq, wait=False)

    @classmethod
    def _runDefSeqROIs(cls, inProt, mode=0):
        outLab = 'outputSequence' if mode == 0 else 'outputVariants'
        inList = defSeqROIsSeq if mode == 0 else defSeqROIsVar

        protDefSeqROIs = cls.newProtocol(
            ProtDefineSeqROI,
            chooseInput=mode, inROIs=inList
        )
        if mode == 0:
            protDefSeqROIs.inputSequence.set(inProt)
            protDefSeqROIs.inputSequence.setExtended(outLab)
        else:
            protDefSeqROIs.inputSequenceVariants.set(inProt)
            protDefSeqROIs.inputSequenceVariants.setExtended(outLab)

        cls.proj.launchProtocol(protDefSeqROIs, wait=False)
        return protDefSeqROIs

    def test(self):
        prots = []
        for i in range(2):
            prots.append(self._runDefSeqROIs(inProt=self.inProts[i], mode=i))

        for p in prots:
            self._waitOutput(p, 'outputROIs', sleepTime=10)
            self.assertIsNotNone(getattr(p, 'outputROIs', None))


class TestGenerateSequences(TestImportVariants):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls._runImportVariants()
        cls._waitOutput(cls.protImportVariants, 'outputVariants', sleepTime=5)

    @classmethod
    def _runGenerateSequences(cls, inProt):

        protGenSeqs = cls.newProtocol(
            ProtChemGenerateVariants,
            toMutateList=genSeqsStr
        )
        protGenSeqs.inputSequenceVariants.set(inProt)
        protGenSeqs.inputSequenceVariants.setExtended('outputVariants')

        cls.proj.launchProtocol(protGenSeqs, wait=False)
        return protGenSeqs

    def test(self):
        p = self._runGenerateSequences(self.protImportVariants)

        self._waitOutput(p, 'outputSequences', sleepTime=10)
        self.assertIsNotNone(getattr(p, 'outputSequences', None))


class TestExtractSequenceROIs(TestGenerateSequences):
    @classmethod
    def _runExtractROIs(cls, inProt, fThres = 0.4, minSize=3):

        protExtSeqROIs = cls.newProtocol(
            ProtExtractSeqsROI,
            thres=0.2, flexThres=fThres, minSize=minSize
        )
        protExtSeqROIs.inputSequences.set(inProt)
        protExtSeqROIs.inputSequences.setExtended('outputSequences')

        cls.proj.launchProtocol(protExtSeqROIs, wait=False)
        return protExtSeqROIs

    def test(self):
        pGen = self._runGenerateSequences(self.protImportVariants)
        self._waitOutput(pGen, 'outputSequences', sleepTime=10)

        p = self._runExtractROIs(pGen)

        self._waitOutput(p, 'outputROIs', sleepTime=10)
        self.assertIsNotNone(getattr(p, 'outputROIs', None))

class TestOperateSeqROIs(TestExtractSequenceROIs, TestDefineSequenceROIs):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

        cls._runImportVariants()
        cls._waitOutput(cls.protImportVariants, 'outputVariants', sleepTime=5)

    @classmethod
    def _runOperateROIs(cls, inProts):
        protOperateROIs = cls.newProtocol(
            ProtOperateSeqROI,
            operation=1, keepNonOverlaping=False
        )
        protOperateROIs.inputROIsSets.set(inProts)
        for i in range(len(protOperateROIs.inputROIsSets)):
            protOperateROIs.inputROIsSets[i].setExtended('outputROIs')

        cls.proj.launchProtocol(protOperateROIs, wait=False)
        return protOperateROIs

    def test(self):
        pGen = self._runGenerateSequences(self.protImportVariants)
        pDef = self._runDefSeqROIs(inProt=self.protImportVariants, mode=1)

        self._waitOutput(pGen, 'outputSequences', sleepTime=10)
        pExt = self._runExtractROIs(pGen, fThres=0.2, minSize=1)

        self._waitOutput(pExt, 'outputROIs', sleepTime=10)
        self._waitOutput(pDef, 'outputROIs', sleepTime=10)

        pOp = self._runOperateROIs([pExt, pDef])
        self._waitOutput(pOp, 'outputROIs', sleepTime=10)
        self.assertIsNotNone(getattr(pOp, 'outputROIs', None))

class TestDefineSetSequences(TestDefineSequenceROIs):
    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dsModBuild = DataSet.getDataSet('model_building_tutorial')

        cls._runImportSequence()
        cls._runImportPDB()

        cls._waitOutput(cls.protImportSeq, 'outputSequence', sleepTime=5)
        cls._waitOutput(cls.protImportPDB, 'outputPdb', sleepTime=5)

    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='4erf')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _writeSeqFiles(cls, inputObj, seqFile, chainId, name):

        seqFile = cls.proj.getPath(seqFile)

        ASH = AtomicStructHandler()
        if issubclass(type(inputObj), str):
            pdbFile = ASH.readFromPDBDatabase(inputObj, type='mmCif', dir=cls.proj.getTmpPath())
            ASH.read(pdbFile)
            seq = str(ASH.getSequenceFromChain(modelID=0, chainID=chainId))

        elif issubclass(type(inputObj), AtomStruct):
            ASH.read(inputObj.getFileName())
            seq = str(ASH.getSequenceFromChain(modelID=0, chainID=chainId))

        elif issubclass(type(inputObj), Sequence):
            seq = inputObj.getSequence()

        with open(seqFile, 'w') as f:
            f.write('>{}\n{}\n'.format(name, seq))

    @classmethod
    def _runDefineSetSequences(cls):
        protDefSeqs = cls.newProtocol(
            ProtDefineSetOfSequences,
            inputList=defSetSeqs
        )

        inpObjs = [cls.protImportSeq.outputSequence, cls.protImportPDB.outputPdb, '5ni1']
        for i, inObj in enumerate(inpObjs):
            cls._writeSeqFiles(inObj, defSetFiles[i], defSetChains[i], names[i])
            if i < 2:
                protDefSeqs.inputPointers.append(Pointer(inObj))

        cls.proj.launchProtocol(protDefSeqs, wait=False)
        return protDefSeqs

    def test(self):
        pDef = self._runDefineSetSequences()
        self._waitOutput(pDef, 'outputSequences', sleepTime=10)
        self.assertIsNotNone(getattr(pDef, 'outputSequences', None))


class TestMapSeqROIs(TestDefineSetSequences):
    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='6vsb')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runMapROIs(cls, inProts):
        protMapROIs = cls.newProtocol(
            ProtMapSequenceROI,
            chain_name='{"model": 0, "chain": "A", "residues": 92}'
        )
        protMapROIs.inputSequenceROIs.set(inProts[0])
        protMapROIs.inputSequenceROIs.setExtended('outputROIs')

        protMapROIs.inputAtomStruct.set(inProts[1])
        protMapROIs.inputAtomStruct.setExtended('outputPdb')

        cls.proj.launchProtocol(protMapROIs, wait=False)
        return protMapROIs
    
    def test(self):
        pDef = self._runDefSeqROIs(inProt=self.protImportSeq, mode=0)
        self._waitOutput(pDef, 'outputROIs', sleepTime=5)

        pMap = self._runMapROIs(inProts=[pDef, self.protImportPDB])
        self._waitOutput(pMap, 'outputStructROIs', sleepTime=5)
        self.assertIsNotNone(getattr(pMap, 'outputStructROIs', None))


class TestPairwiseAlign(TestDefineSetSequences):
    @classmethod
    def _runImportPDB(cls):
        cls.protImportPDB = cls.newProtocol(
            ProtImportPdb,
            inputPdbData=0, pdbId='6vsb')
        cls.proj.launchProtocol(cls.protImportPDB, wait=False)

    @classmethod
    def _runPairwiseAlign(cls, inProts):
        protAlign = cls.newProtocol(
            ProtChemPairWiseAlignment,
            condAtomStruct1=False, chain_name2='{"model": 0, "chain": "A", "residues": 92}'
        )
        protAlign.inputSequence1.set(inProts[0])
        protAlign.inputSequence1.setExtended('outputSequence')

        protAlign.inputAtomStruct2.set(inProts[1])
        protAlign.inputAtomStruct2.setExtended('outputPdb')

        cls.proj.launchProtocol(protAlign, wait=False)
        return protAlign

    def test(self):
        pAlign = self._runPairwiseAlign(inProts=[self.protImportSeq, self.protImportPDB])
        self._waitOutput(pAlign, 'outputSequences', sleepTime=5)
        self.assertIsNotNone(getattr(pAlign, 'outputSequences', None))


class TestMultipleAlignSequences(TestImportSequences):
    def _runAlignment(self, program=0):
        protAlignSeqs = self.newProtocol(
          ProtChemMultipleSequenceAlignment,
          inputSequences=self.protImportSeqs.outputSequences,
          programList=program)
        self.proj.launchProtocol(protAlignSeqs, wait=False)
        return protAlignSeqs

    def test(self):
        self._runImportSeqs()

        protCLUSTALO = self._runAlignment(CLUSTALO)
        protMUSCLE = self._runAlignment(MUSCLE)
        protMAFFT = self._runAlignment(MAFFT)

        self._waitOutput(protCLUSTALO, 'outputSequences', sleepTime=5)
        self._waitOutput(protMUSCLE, 'outputSequences', sleepTime=5)
        self._waitOutput(protMAFFT, 'outputSequences', sleepTime=5)

        self.assertIsNotNone(protCLUSTALO.outputSequences)
        self.assertIsNotNone(protMUSCLE.outputSequences)
        self.assertIsNotNone(protMAFFT.outputSequences)
