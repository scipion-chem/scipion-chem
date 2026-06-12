# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
SCORCH2 (SC2) is a machine learning rescoring model designed for interaction-based virtual screening. (SC2, https://github.com/LinCompbio/SCORCH2)

"""
import csv
import logging
from ctypes.wintypes import SMALL_RECT
from pathlib import Path
from pwem.convert import cifToPdb

from pyworkflow.object import Float
from pyworkflow.protocol import params, STEPS_PARALLEL
from pwem.protocols import EMProtocol

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import os, shutil, re, runOpenBabel, makeSubsets, insistentRun
from pwchem import Plugin, SCORCH2_DIC


currentDir = Path(__file__).parent.resolve()



class ProtocolSCORCH2(EMProtocol):
    """
    AI Generated:

This protocol performs molecular rescoring of protein–ligand docking poses
using SCORCH2 (SC2), a machine-learning-based scoring framework designed for
interaction-aware virtual screening.

SCORCH2 combines multiple learned scoring models to estimate the quality of
ligand binding poses, enabling improved ranking of docked compounds beyond
classical docking scores.

Overview
--------
The protocol takes a SetOfSmallMolecules containing docked ligand poses and
computes refined binding scores using SCORCH2 models.

It supports both GPU and CPU execution and can optionally aggregate multiple
poses per compound to produce a single representative score.

SCORCH2 Scoring Components
--------------------------
SCORCH2 evaluates ligand poses using two complementary models:

- SC2-PS (Protein Score):
  Predicts protein–ligand interaction quality based on structural and
  physicochemical features of the binding interface.

- SC2-PB (Pose/Bond Score):
  Evaluates ligand pose plausibility and internal chemical consistency.

These two components are combined using user-defined weights:

    final_score = psWeight * SC2-PS + pbWeight * SC2-PB

Input Handling
--------------
- Input ligands are provided as a SetOfSmallMolecules
- The dataset is split into batches for parallel processing
- Each batch is processed independently to improve scalability

Protein Handling
----------------
- The receptor (protein) is extracted from the input dataset
- It is copied and prepared for compatibility with SCORCH2
- Ligand and protein directories are organized per batch

File Conversion
---------------
- Input ligand and protein files are converted to PDBQT format if needed
- CIF and PDB formats are supported via OpenBabel conversion
- Temporary files are cleaned after conversion to ensure consistency

Workflow
--------
1. Input docked ligand poses and associated receptor
2. Organize protein and ligand structures into batch directories
3. Convert input files to SCORCH2-compatible formats (PDBQT)
4. Run SCORCH2 rescoring script per batch
5. Optionally use GPU acceleration
6. Optionally aggregate multiple poses per compound
7. Collect batch-wise TSV results
8. Parse SC2 scores and map them to ligand identifiers
9. Assign final scores to molecular objects

Output
------
- outputSmallMolecules:
    A copy of the input ligand set enriched with:

    - scorchScore:
        Continuous numeric score representing SCORCH2-predicted binding
        quality. Higher scores indicate more favorable binding poses.

If aggregation is enabled, scores may represent per-compound summaries
rather than per-pose values.

Key Features
------------
- Machine-learning-based rescoring of docking poses
- Combination of protein–ligand and pose-quality models
- GPU acceleration support
- Batch processing for large-scale virtual screening
- Optional aggregation of multiple poses per compound
- Automatic file conversion and preprocessing (PDB/CIF → PDBQT)

Use Cases
---------
- Post-docking rescoring in virtual screening pipelines
- Prioritization of docking hits using ML-based scoring
- Improvement of ranking accuracy over classical docking scores
- Large-scale compound library screening

Notes
-----
- Requires SCORCH2 environment and pretrained models installed
- GPU usage is optional but recommended for large datasets
- Score quality depends on correct ligand naming and preprocessing
- Batch processing may slightly vary results depending on aggregation mode
"""
    _label = 'SCORCH2 rescoring'
    _defaultName = 'prot'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('useGPU', params.BooleanParam, default=True, label="Use GPU: ",
                      help='Whether to use GPU or not. (Unable to choose the GPU id).')
        iGroup = form.addGroup('Input')
        # Pre-extracted features
        iGroup.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input ligand: ',
                        help='Input docked small molecules to rescore')
        iGroup.addParam('batchSize', params.IntParam, default=500, expertLevel=params.LEVEL_ADVANCED,
                        label='Batch size: ', help="Size of the batches send to rescore")
        # Aggregate
        iGroup.addParam('aggregate', params.BooleanParam, default=False,
                        label="Aggregate results: ",
                        help='If Yes, the best pose will be selected per compound with aggregation metadata. If No, all poses will be scored individually and ranked by SC2 score.')

        iGroup.addParam('psWeight', params.FloatParam, default=0.7, expertLevel=params.LEVEL_ADVANCED,
                        label='PS weight: ', help="Weight for SC2-PS predictions")
        iGroup.addParam('pbWeight', params.FloatParam, default=0.3, expertLevel=params.LEVEL_ADVANCED,
                        label='PB weight: ', help="Weight for SC2-PB predictions")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        rStep = self._insertFunctionStep(self.organizeInputStep)

        sSteps = []
        for it in range(self.getNBatches()):
            cStep = self._insertFunctionStep(self.convertInputStep, it, prerequisites=[rStep])
            sSteps += [self._insertFunctionStep(self.scorchStep, it, prerequisites=[cStep])]
        oStep = self._insertFunctionStep(self.createOutputStep, prerequisites=sSteps)

    def organizeInputStep(self):
        # Receptor
        proteinDir = self.getProtDir()
        proteinDir.mkdir(parents=True, exist_ok=True)

        protein = self.inputSmallMolecules.get().getProteinFile()
        proteinPath = Path(protein)

        proteinFile = proteinDir / f"{self._defaultName}_protein{proteinPath.suffix}"
        shutil.copy(proteinPath, proteinFile)

        proteinFiles = list(proteinDir.glob("*"))
        if proteinFiles:
            self.convertFiles(proteinFiles, proteinDir)
            self.removePdbFiles(proteinDir)
        else:
            logging.warning("No protein files found.")

        # Ligands
        nBatches = self.getNBatches()
        molSubSets = makeSubsets(self.inputSmallMolecules.get(), nBatches, True)
        for it, subset in enumerate(molSubSets):
            moleculeDir = self.getMoleculesDir(it)
            moleculeDir.mkdir(parents=True, exist_ok=True)

            ligandOutDir = moleculeDir / self._defaultName
            ligandOutDir.mkdir(parents=True, exist_ok=True)

            for i, ligand in enumerate(subset, start=1):
                ligandPath = Path(ligand.getPoseFile())
                origName = ligandPath.name
                newName = f"{self._defaultName}_{origName}"

                dest = ligandOutDir / newName
                shutil.copy(ligandPath, dest)


    def convertInputStep(self, it):
        ligandDir = self.getMoleculesDir(it)

        ligandFiles = []
        for subfolder in ligandDir.iterdir():
            if subfolder.is_dir():
                ligandFiles.extend(f for f in subfolder.glob("*") if f.is_file())
        if ligandFiles:
            self.convertFiles(ligandFiles, os.path.abspath(subfolder))
            self.removePdbFiles((subfolder))
        else:
            logging.warning("No ligand files found.")


    def scorchStep(self, it):
        oriProteinDir = self.getProtDir()
        proteinDir = self.getProtDir(it)
        shutil.copytree(oriProteinDir, proteinDir)

        ligandDir = self.getMoleculesDir(it)

        scriptRescoringDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'SCORCH2'))
        modelsDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'scorchModels/models'))

        script = os.path.abspath(os.path.join(scriptRescoringDir, 'scorch2_rescoring.py'))
        sc2PsModel = os.path.abspath(os.path.join(modelsDir, 'sc2_ps.xgb'))
        sc2PbModel = os.path.abspath(os.path.join(modelsDir, 'sc2_pb.xgb'))
        psScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_ps_scaler'))
        pbScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_pb_scaler'))

        resDir = self.getResultsDir(it)
        resDir.mkdir(exist_ok=True)
        outputFileFull = self.getBatchDir(it) / 'scorch2_results.tsv'
        args = [
                str(script),
                "--protein-dir", str(proteinDir.resolve()),
                "--ligand-dir", str(ligandDir.resolve()),
                "--sc2_ps_model", str(sc2PsModel),
                "--sc2_pb_model", str(sc2PbModel),
                "--ps_scaler", str(psScaler),
                "--pb_scaler", str(pbScaler),
                "--output", str(outputFileFull.resolve()),
                "--ps_weight", str(self.psWeight.get()),
                "--pb_weight", str(self.pbWeight.get()),
                '--res-dir', str(resDir.resolve()),
                "--num-cores", '1',
            ]

        if self.aggregate.get():
            args.append("--aggregate")

        if self.useGPU.get():
            args.append("--gpu")

        insistentRun(self, 'python', args,
                     envDic=SCORCH2_DIC, nMax=5, cwd=os.path.abspath(Plugin.getVar(SCORCH2_DIC['home'])), sleepTime=5)

    def createOutputStep(self):
        inMols = self.inputSmallMolecules.get()
        newMols = SetOfSmallMolecules.createCopy(inMols, self._getPath(), copyInfo=True)
        scoresDict = self.readScoresTSV()
        for mol in inMols:
            newMol = mol.clone()
            newMol.scorchScore = Float()

            molName = Path(newMol.getPoseFile()).stem
            if molName in scoresDict:
                newMol.setAttributeValue('scorchScore', scoresDict[molName])
            else:
                newMol.setAttributeValue('scorchScore', None)
            newMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=newMols)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["The results have been saved in the extra folder as scorch2.results.tsv.", "Features and normalized features have been saved in the extra/results folder, along with the aggregated results and target results."]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputSmallMolecules.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getBatchDir(self, it):
        return Path(self._getExtraPath()) / f"batch_{it}"

    def getProtDir(self, it=None):
        if it is None:
            return Path(self._getExtraPath()) / "protein"
        return self.getBatchDir(it) / "protein"

    def getMoleculesDir(self, it):
        return self.getBatchDir(it) / "molecule"

    def getResultsDir(self, it):
        return self.getBatchDir(it) / f"results"

    def getPDBId(self):
        protein = self.inputSmallMolecules.get().getProteinFile()
        proteinPath = Path(protein)
        return proteinPath.stem

    def checkPdbqtFiles(self, directory):
        """Check if files are PDBQT"""
        files = list(directory.glob("*"))
        for f in files:
            if f.suffix.lower() != ".pdbqt":
                return False, files
            else:
                return True, files

    def convertFiles(self, fileList, baseDir):
        """Convert PDB or CIF to PDBQT, keeping output in the same folder as the input file"""
        oFiles = []
        for file in fileList:
            suffix = file.suffix.lower()
            basename = file.stem

            if suffix == ".pdbqt":
                continue

            pdbqtFile = Path(baseDir) / f"{basename}.pdbqt"
            outputPath = str(pdbqtFile.resolve())
            oFiles.append(outputPath)
            inputPath = str(file.resolve())

            if suffix == ".cif":
                pdbFile = Path(baseDir) / f"{basename}.pdb"
                cifToPdb(inputPath, str(pdbFile.resolve()))
                inputPath = str(pdbFile.resolve())
                suffix = '.pdb'

            args = f"-i{suffix[1:]} {inputPath} -opdbqt -O {outputPath}"
            runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
        return oFiles


    def removePdbFiles(self, directory):
        """Removes .pdb files only if their corresponding .pdbqt exists in the same directory."""
        for pdbFile in directory.rglob("*.pdb"):
            pdbqtFile = pdbFile.with_suffix(".pdbqt")
            if pdbqtFile.exists():
                try:
                    pdbFile.unlink()
                except Exception as e:
                    logging.warning(f"Could not delete {pdbFile.name}: {e}")

    def readScoresTSV(self):
        scoreDict = {}
        for it in range(self.getNBatches()):
            resultsTsv = self.getBatchDir(it) / "scorch2_results.tsv"
            with open(resultsTsv, "r", encoding="utf-8-sig") as f:
                lines = [line for line in f.readlines() if not line.startswith("#") and line.strip()]

            reader = csv.DictReader(lines, delimiter=",")

            for row in reader:
                compounId = row["compound_id"].strip()
                molName = compounId.split("_", 1)[1].replace(".pdbqt", "")
                score = float(row["sc2_score"])
                scoreDict[molName] = score
        return scoreDict

    def getNBatches(self):
        nThreads = self.numberOfThreads.get() - 1
        nThreads = 1 if nThreads < 1 else nThreads

        nMols = len(self.inputSmallMolecules.get())
        nBatches = (nMols // self.batchSize.get()) + 1

        maxThreads = max(nThreads, nBatches)

        return min(maxThreads, nMols)

