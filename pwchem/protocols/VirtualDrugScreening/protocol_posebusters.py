# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
# *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to perform tests oer predicted molecule poses (https://github.com/maabuu/posebusters).

"""
import os.path

import csv, json, shutil

from pyworkflow.object import Boolean
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC, OPENBABEL_DIC



class ProtocolPoseBusters(EMProtocol):
    """
    AI Generated:

This protocol evaluates the quality and chemical plausibility of docked ligand poses
using the PoseBusters framework (https://posebusters.readthedocs.io/).

PoseBusters applies a series of geometric, chemical, and steric validation tests to
predict whether a docked ligand pose is physically realistic and chemically valid.

Overview
--------
Docking pose validation is a critical step in structure-based drug discovery to ensure
that predicted ligand binding modes are not only energetically favorable but also
chemically and sterically plausible.

This protocol automates the execution of PoseBusters tests on a set of docked ligands
and filters results based on the number of passed validation criteria.

Test Categories
---------------
PoseBusters evaluates ligand poses using multiple types of checks, including:

- Chemical validity:
  * Sanity checks of molecular structure
  * Atom connectivity and valence correctness
  * Radical and charge consistency

- Geometric consistency:
  * Bond lengths and angles
  * Ring planarity
  * Stereochemistry correctness

- Steric and spatial feasibility:
  * Clashes within ligand
  * Protein-ligand steric clashes
  * Distance to receptor atoms or cofactors
  * Volume overlap with protein or water molecules

- Reference-based metrics (optional):
  * RMSD to true (crystal) ligand
  * Centroid displacement

Workflow
--------
1. Input a set of docked small molecules (poses).
2. Optionally provide:
   - A reference (true) ligand structure
   - A receptor (protein structure)
3. Organize and split input ligands into processing batches.
4. Convert ligand structures to a PoseBusters-compatible format (PDB).
5. Run PoseBusters CLI on each batch of poses.
6. Collect and merge CSV results from all batches.
7. Parse per-pose test outcomes.
8. Filter poses based on a minimum number of passed tests.
9. Annotate surviving poses with individual test results.

Filtering Strategy
------------------
Each pose is evaluated against a set of binary pass/fail tests.
A pose is retained if it passes at least:

    testsPassed ≥ N

where N is configurable and depends on whether protein and/or true ligand information
is used.

Output
------
- outputSmallMolecules:
    A filtered SetOfSmallMolecules containing only poses that pass the minimum
    number of PoseBusters validation tests.

- outputSmallMoleculesFailed (optional):
    A SetOfSmallMolecules with the poses that do not pass them, so that nothing
    is lost: the two sets together hold all the input poses. They carry the
    result of each test, to see which ones they fail.

    Each retained molecule includes:
    - Boolean flags for each evaluated test (TRUE/FALSE)
    - Link to the PoseBusters result file
    - Optional reference-based metrics if enabled

Key Features
------------
- Parallel execution over ligand batches
- Optional receptor-based validation, against the receptor each pose was docked into: the one of the
  molecule when it carries its own (flexible docking, cofolding) and the general one of the set otherwise
- Optional ground-truth ligand comparison
- Automatic conversion of input formats (PDBQT → PDB, CIF → PDB)
- Flexible filtering based on test subsets

Use Cases
---------
- Post-docking validation in virtual screening pipelines
- Filtering false-positive docking poses
- Quality control of large-scale docking campaigns
- Benchmarking docking protocols against structural criteria

Notes
-----
- Requires PoseBusters CLI installed and available in the execution environment.
- Results depend strongly on input pose quality and receptor accuracy.
- Protein-based tests are only meaningful if receptor structure is provided.
"""
    _label = 'PoseBusters docking tests'

    NTESTS = {'simple': 12, 'receptor': 22, 'trueMol': 18, 'trueMolAndReceptor': 28}

    DEF_COL_THRES = {'energy_ratio': 2.0, 'number_clashes': 0, 'aromatic_ring_maximum_distance_from_plane': 0.05,
                     'kabsch_rmsd': 2.5, 'centroid_distance': 12.0,
                     'volume_overlap_protein': 0.05, 'smallest_distance_protein': 6.0}
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        iGroup = form.addGroup('Input')
        iGroup.addParam('inputMoleculesSet', params.PointerParam, pointerClass='SetOfSmallMolecules', allowsNull=False,
                        label="Input Docked Small Molecules: ",
                        help='Select the docked molecules to be tested.')

        iGroup.addParam('useTrueMol', params.BooleanParam, default=False, label="Use true molecule: ",
                        help='Choose whether to use the true molecule in the analysis.')
        iGroup.addParam('inputMoleculesRefSet', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label="Input Reference Small Molecules: ", condition='useTrueMol', allowsNull=True,
                        help='Select the reference molecules to be tested against.')
        iGroup.addParam('molTrue', params.StringParam, label='True molecule: ', condition='useTrueMol',
                        help='Choose the ground truth molecule (crystal ligand).')

        iGroup.addParam('molRec', params.BooleanParam, default=True, label="Use receptor to run tests: ",
                        help='Choose whether to output use receptor to run tests. (If True, it will be used to run '
                             'all default tests that require the protein as input.)')
        iGroup.addParam('batchSize', params.IntParam, label='Batch size: ', default=500,
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Maximum batch size to execute')

        fGroup = form.addGroup('Tests')
        fGroup.addParam('testsPassed', params.IntParam, label='Min. tests passed: ', default=self.NTESTS['receptor'],
                        help='Minimum number of tests passed to keep each molecule.')
        fGroup.addParam('outputFailed', params.BooleanParam, label='Create a set with the discarded molecules: ',
                        default=False,
                        help='Create an extra set with the poses that do not pass PoseBusters.')
        fGroup.addParam('doSpecific', params.BooleanParam, label='Check only specific tests: ', default='',
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Check only a subset of specific tests on filtering')
        fGroup.addParam('specificTests', params.TextParam, width=70, height=15, expertLevel=params.LEVEL_ADVANCED,
                        default='', label='List of tests to consider: ', condition='doSpecific',
                        help='List of tests to consider in the filtering.\nRemember to modify the "Min. tests passed" '
                             'param accordingly.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        pSteps = []
        oStep = self._insertFunctionStep(self.organizeInputStep)
        for it in range(self.getNThreads()):
            cStep = self._insertFunctionStep(self.convertInputStep, it, prerequisites=[oStep])
            pSteps += [self._insertFunctionStep(self.poseBustersStep, it, prerequisites=[cStep])]
        self._insertFunctionStep(self.createOutputStep, prerequisites=pSteps)

    def organizeInputStep(self):
        if self.useTrueMol.get():
            molTrue = self.getTrueMol()
            trueFile = self.convertFormat(molTrue, type='crystal')
            shutil.copy(trueFile, self.getTrueMolFile())

        if self.molRec.get():
            self.organizeReceptors()

        nThreads = self.getNThreads()
        molSubSets = makeSubsets(self.inputMoleculesSet.get(), nThreads, True)

        for i, subset in enumerate(molSubSets):
            for mol in subset:
                inDir = self.getInMolsDir(i)
                os.makedirs(inDir, exist_ok=True)
                copyPath = os.path.join(inDir, os.path.basename(mol.getPoseFile()))
                if not os.path.exists(copyPath):
                    os.link(mol.getPoseFile(), copyPath)

    def convertInputStep(self, it):
        inDir, convDir = self.getInMolsDir(it), self.getConvMolsDir(it)
        os.makedirs(convDir, exist_ok=True)
        args = f' --multiFiles -iD "{inDir}" --pattern "*" -of pdb --outputDir "{convDir}"'
        pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=convDir)

    def poseBustersStep(self, it):
        resultsFile = self.getResultsFile(it)

        args = ['bust', f'-t {self.writePosesTable(it)}', '--outfmt csv', f'--output {os.path.abspath(resultsFile)}',
                '--full-report', '--max-workers 1']

        Plugin.runCondaCommand(
            self, args=" ".join(args), condaDic=SCORCH2_DIC, program="",
            cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
        )

    def createOutputStep(self):
        resultsFile = self.getResultsFile()
        if not os.path.exists(resultsFile):
            concatThreadFiles(resultsFile, self._getPath())
        resTestsDic = self.parseResultRows(resultsFile)

        molSet = self.inputMoleculesSet.get()
        newMols = SetOfSmallMolecules.createCopy(molSet, self._getPath(), copyInfo=True)
        failedMols = SetOfSmallMolecules.createCopy(molSet, self._getPath(), suffix='_failed', copyInfo=True) \
            if self.outputFailed.get() else None

        for mol in molSet:
            newMol = mol.clone()
            newMol.PoseBusters_file = String(resultsFile)

            resDic = resTestsDic.get(getBaseName(mol.getPoseFile()))
            if resDic is not None:
                for testName, testVal in resDic.items():
                    newMol.__setattr__(testName, Boolean(testVal.upper() == 'TRUE'))

            if resDic is not None and self.checkPassedTests(resDic) >= self.testsPassed.get():
                newMols.append(newMol)
            elif failedMols is not None:
                failedMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=newMols)
        if failedMols is not None:
            if len(failedMols) > 0:
                self._defineOutputs(outputSmallMoleculesFailed=failedMols)
            else:
                print('Every molecule passed the tests: there is no set of discarded molecules to create')
                failedPath = failedMols.getSetPath()
                failedMols.close()
                os.remove(failedPath)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ['A results file has been created in the extra folder with the results of the tests.']
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputMoleculesSet.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        if self.molRec.get() and not molSet.getProteinFile():
            noRec = [mol.getUniqueName() for mol in molSet if not mol.getProteinFile()]
            if noRec:
                validations += ['The set has no receptor and neither have {} of its molecules (e.g. {}): '
                                'uncheck "Use receptor to run tests" or use docked molecules '
                                'that carry it'.format(len(noRec), ', '.join(noRec[:3]))]

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- WIZARD functions ---------------------------------
    def getNTests(self):
        return len(self.getDefTests().split('\n'))

    def getDefTests(self):
        cols = ["mol_pred_loaded", "sanitization", "inchi_convertible", "all_atoms_connected", "no_radicals",
                "bond_lengths",	"bond_angles", "internal_steric_clash", "aromatic_ring_flatness", "non-aromatic_ring_non-flatness",
                "double_bond_flatness", "internal_energy"]
        if self.useTrueMol.get():
            cols += ['mol_true_loaded', 'molecular_bonds', 'tetrahedral_chirality', 'molecular_formula', 'rmsd_≤_2å',
                     'double_bond_stereochemistry']
        if self.molRec.get():
            cols += ['mol_cond_loaded',
                     'minimum_distance_to_waters', 'volume_overlap_with_waters',
                     'minimum_distance_to_protein', 'volume_overlap_with_protein',
                     'minimum_distance_to_inorganic_cofactors', 'volume_overlap_with_inorganic_cofactors',
                     'minimum_distance_to_organic_cofactors', 'volume_overlap_with_organic_cofactors',
                     'protein-ligand_maximum_distance']
        cols = '\n'.join(cols)
        return cols

    def getColumns(self):
        cols = ['energy_ratio', 'number_clashes', 'aromatic_ring_maximum_distance_from_plane']
        if self.useTrueMol.get():
            cols += ['kabsch_rmsd', 'centroid_distance']
        if self.molRec.get():
            cols += ['volume_overlap_protein', 'smallest_distance_protein']
        return cols

    def getFilterDefValue(self):
        return self.DEF_COL_THRES[self.filterCol.get()]

    # --------------------------- UTILS functions -----------------------------------

    def getPoseReceptor(self, mol):
        '''Receptor a pose was docked into: the one of the molecule when it carries its own (flexible receptors,
        cofolding), else the general one of the set'''
        return mol.getProteinFile() or self.inputMoleculesSet.get().getProteinFile()

    def organizeReceptors(self):
        '''Write the receptors the poses were docked into as PDB, each of them once, and which of them goes with each
        pose: {pose name: receptor file}'''
        poseRecs = {}
        for mol in self.inputMoleculesSet.get():
            recFile = self.getPoseReceptor(mol)
            if not recFile:
                raise ValueError(f'No receptor found for {mol}: neither the molecule nor its set carry one')
            poseRecs[getBaseName(mol.getPoseFile())] = os.path.abspath(recFile)

        pdbFiles = self.convertReceptors(sorted(set(poseRecs.values())))
        with open(self.getPoseRecsFile(), 'w') as f:
            json.dump({pose: pdbFiles[recFile] for pose, recFile in poseRecs.items()}, f, indent=1)

    def convertReceptors(self, recFiles):
        '''Receptors in the format PoseBusters reads: {receptor file: its PDB}. The pdbqt ones are all converted in a single call'''
        recDir, pdbqtDir = self.getRecsDir(), os.path.abspath(self._getTmpPath('receptorsPDBQT'))
        os.makedirs(recDir, exist_ok=True)

        pdbFiles, nPdbqt = {}, 0
        for i, recFile in enumerate(recFiles):
            ext = os.path.splitext(recFile)[1].lower()
            if ext == '.cif':
                pdbFiles[recFile] = os.path.join(recDir, f'rec_{i}.pdb')
                cifToPdb(recFile, pdbFiles[recFile])
            elif ext == '.pdbqt':
                os.makedirs(pdbqtDir, exist_ok=True)
                os.symlink(recFile, os.path.join(pdbqtDir, f'rec_{i}.pdbqt'))
                pdbFiles[recFile] = os.path.join(recDir, f'rec_{i}.pdb')
                nPdbqt += 1
            else:
                pdbFiles[recFile] = recFile

        if nPdbqt:
            args = f' --multiFiles -iD "{pdbqtDir}" --pattern "*.pdbqt" -of pdb --outputDir "{recDir}"'
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=recDir)

        # A receptor that was not written would only show up as every test on it failing
        notConverted = [recFile for recFile, pdbFile in pdbFiles.items() if not os.path.exists(pdbFile)]
        if notConverted:
            raise FileNotFoundError(f'{len(notConverted)} receptors could not be converted to PDB, e.g. '
                                    f'{notConverted[0]}')
        return pdbFiles

    def getPoseReceptors(self):
        '''{pose name: receptor PDB} written when organizing the input'''
        with open(self.getPoseRecsFile()) as f:
            return json.load(f)

    def writePosesTable(self, it):
        '''Table PoseBusters reads for a batch of poses: a row for each of them, with the receptor it was docked into
        and the true molecule, when they are used'''
        useRec, useTrue = self.molRec.get(), self.useTrueMol.get()
        poseRecs = self.getPoseReceptors() if useRec else {}

        convDir = self.getConvMolsDir(it)
        tableFile = os.path.abspath(self._getExtraPath(f'poses_{it}.csv'))
        with open(tableFile, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(['mol_pred'] + (['mol_true'] if useTrue else []) + (['mol_cond'] if useRec else []))
            for poseFile in sorted(os.listdir(convDir)):
                row = [os.path.join(convDir, poseFile)]
                if useTrue:
                    row.append(self.getTrueMolFile())
                if useRec:
                    poseName = getBaseName(poseFile)
                    if poseName not in poseRecs:
                        raise ValueError(f'No receptor was organized for the pose {poseName}')
                    row.append(poseRecs[poseName])
                writer.writerow(row)
        return tableFile

    def getResultsFile(self, it=None):
        if it is None:
            return self._getPath('results.csv')
        return self._getPath(f'results_{it}.csv')
    
    def getTestToPass(self):
        if self.doSpecific.get():
            tests = self.specificTests.get()
        else:
            tests = self.getDefTests()
        return tests.split('\n')
    
    def checkPassedTests(self, resDic):
        n = 0
        for testName, testValue in resDic.items():
            if testValue.upper() == 'TRUE':
                n += 1
        return n
    
    def parseResultRows(self, resultsFile):
        resDic = {}
        csvRows = self.getFileInfo(resultsFile)
        testsToPass = self.getTestToPass()
        for poseName, row in csvRows.items():
            resDic[poseName] = {}
            for test in testsToPass:
                csvVal = row.get(test)
                resDic[poseName][test] = csvVal

        return resDic

    def getFileInfo(self, resultsFile):
        csvRows = {}
        with open(resultsFile) as f:
            reader = csv.DictReader(f, delimiter=',')
            for row in reader:
                poseFile = row['molecule'] if row['molecule'] else row['file']
                csvRows[getBaseName(poseFile)] = row

        return csvRows

    def convertFormat(self, molPred, type=''):
        if type in ('AtomStruct', 'crystal'):
            basename = os.path.basename(molPred.getFileName()).split('.')[0]
            file = molPred.getFileName()
        elif type == 'file':
            file = molPred
            basename = os.path.splitext(os.path.basename(file))[0]
        else:
            basename = os.path.basename(molPred.getPoseFile()).split('.')[0]
            file = molPred.getPoseFile()

        if file.endswith('.cif'):
            inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            cifToPdb(file, inpFile)
        elif (file.endswith('.pdbqt')):
            inpFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            pdbqt2other(self, file, inpFile)
        else:
            inpFile = file
        return inpFile

    def getTrueMol(self):
        myMol = None
        for mol in self.inputMoleculesRefSet.get():
            if mol.__str__() == self.molTrue.get():
                myMol = mol.clone()
                break

        if myMol is None:
            print('The input ligand is not found')
        return myMol

    def getTrueMolFile(self):
        return os.path.abspath(self._getExtraPath('trueMol.pdb'))

    def getRecsDir(self):
        return os.path.abspath(self._getExtraPath('receptors'))

    def getPoseRecsFile(self):
        return os.path.abspath(self._getExtraPath('poseReceptors.json'))

    def getInMolsDir(self, it):
        return os.path.abspath(self._getTmpPath(f'inputMols_{it}'))

    def getConvMolsDir(self, it):
        return os.path.abspath(self._getExtraPath(f'convertedMols_{it}'))

    def getNThreads(self):
        nThreads = self.numberOfThreads.get() - 1
        nThreads = 1 if nThreads < 1 else nThreads

        nMols = len(self.inputMoleculesSet.get())
        nBatches = (nMols // self.batchSize.get()) + 1

        maxThreads = max(nThreads, nBatches)

        return min(maxThreads, nMols)