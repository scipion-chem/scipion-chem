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
from pwchem.utils import os, shutil, re, runOpenBabel
from pwchem import Plugin, SCORCH2_DIC


currentDir = Path(__file__).parent.resolve()



class ProtocolSCORCH2(EMProtocol):
    """
    Computes the best poses of protein-ligand interactions using SCORCH2: https://github.com/LinCompbio/SCORCH2
    """
    _label = 'SCORCH2 rescoring'
    _defaultName = 'prot'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        iGroup = form.addGroup('Input')
        # Pre-extracted features
        iGroup.addParam('useFeatures', params.BooleanParam, default=False,
                       label="Use extracted features: ",
                       help='Choose whether to use pre-extracted features directly.')
        iGroup.addParam('inputFeatures', params.PathParam,
                        label='Input features file: ',
                        condition='useFeatures',
                        help='Input file where the pre-extracted features are stored.')

        iGroup.addParam('inputPDBligandFiles', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input ligand: ', allowsNull=True,
                        condition='not useFeatures',
                        help='Input folder with PDB ligand files.')
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
        preExtracted = self.useFeatures.get()
        if not preExtracted:
            self._insertFunctionStep('moveFilesStep')
            self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('scorchStep')
        self._insertFunctionStep('createOutputStep')
        self._insertFunctionStep('renameFilesStep')

    def convertInputStep(self):
        proteinDir = Path(self._getExtraPath()) / "protein"
        ligandDir = Path(self._getExtraPath()) / "molecule"

        proteinFiles = list(proteinDir.glob("*"))
        if proteinFiles:
            self.convertFiles(proteinFiles, proteinDir)
            self.removePdbFiles(proteinDir)
        else:
            logging.warning("No protein files found.")

        ligandFiles = []
        for subfolder in ligandDir.iterdir():
            if subfolder.is_dir():
                ligandFiles.extend(f for f in subfolder.glob("*") if f.is_file())
        if ligandFiles:
            self.convertFiles(ligandFiles, os.path.abspath(subfolder))
            self.removePdbFiles((subfolder))
        else:
            logging.warning("No ligand files found.")

    def scorchStep(self):
        proteinDir = Path(self._getExtraPath()) / "protein"
        ligandDir = Path(self._getExtraPath()) / "molecule"

        projectDir = Path(self.getWorkingDir()).parent.resolve().parent
        proteinFull = projectDir / proteinDir
        ligandFull = projectDir / ligandDir

        outputFile = self._getExtraPath('scorch2_results.tsv')
        outputFileFull = projectDir / outputFile

        scriptRescoringDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'SCORCH2'))
        modelsDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'scorchModels/models'))

        script = os.path.abspath(os.path.join(scriptRescoringDir, 'scorch2_rescoring.py'))
        sc2PsModel = os.path.abspath(os.path.join(modelsDir, 'sc2_ps.xgb'))
        sc2PbModel = os.path.abspath(os.path.join(modelsDir, 'sc2_pb.xgb'))
        psScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_ps_scaler'))
        pbScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_pb_scaler'))

        preExtracted = self.useFeatures.get()
        if not preExtracted:
            args = [
                str(script),
                "--protein-dir", str(proteinFull),
                "--ligand-dir", str(ligandFull),
                "--sc2_ps_model", str(sc2PsModel),
                "--sc2_pb_model", str(sc2PbModel),
                "--ps_scaler", str(psScaler),
                "--pb_scaler", str(pbScaler),
                "--output", str(outputFileFull),
                "--num-cores", str(self.numberOfThreads.get()),
                "--ps_weight", str(self.psWeight.get()),
                "--pb_weight", str(self.pbWeight.get()),
                "--res-dir", str('results')
            ]
        else:
            args = [
                str(script),
                "--features", str(self.inputFeatures.get()),
                "--sc2_ps_model", str(sc2PsModel),
                "--sc2_pb_model", str(sc2PbModel),
                "--ps_scaler", str(psScaler),
                "--pb_scaler", str(pbScaler),
                "--output", str(outputFileFull),
                "--ps_weight", str(self.psWeight.get()),
                "--pb_weight", str(self.pbWeight.get())
            ]

        if self.aggregate.get():
            args.append("--aggregate")

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=SCORCH2_DIC,
            program="python",
            cwd=os.path.abspath(Plugin.getVar(SCORCH2_DIC['home']))
        )
        if not preExtracted:
            resDir = Path(self._getExtraPath()) / "results"
            shutil.move(os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'results')), resDir)

    def createOutputStep(self):
        mols = self.inputPDBligandFiles.get()
        newMols = SetOfSmallMolecules().create(outputPath=self._getPath())
        scoresDict = self.readScoresTSV()
        for mol in mols:
            newMol = SmallMolecule()

            newMol.copy(mol)
            newMol.scorchScore = Float()
            molName = Path(newMol.getPoseFile()).stem
            if molName in scoresDict:
                newMol.setAttributeValue('scorchScore', scoresDict[molName])
            else:

                newMol.setAttributeValue('scorchScore', None)
            newMols.append(newMol)
        newMols.setDocked(True)
        newMols.proteinFile.set(self.inputPDBligandFiles.get().getProteinFile())
        self._defineOutputs(outputSmallMolecules=newMols)

    def moveFilesStep(self):
        extraPath = Path(self._getExtraPath())

        proteinDir = extraPath / "protein"
        moleculeDir = extraPath / "molecule"
        proteinDir.mkdir(parents=True, exist_ok=True)
        moleculeDir.mkdir(parents=True, exist_ok=True)

        protein = self.inputPDBligandFiles.get().getProteinFile()
        proteinPath = Path(protein)
        pdbId = self.getPDBId()
        #change names so it does not crash with _
        proteinFile = proteinDir / f"{self._defaultName}_protein{proteinPath.suffix}"
        shutil.copy(proteinPath, proteinFile)

        ligands = self.inputPDBligandFiles.get()
        ligandOutDir = moleculeDir / self._defaultName
        ligandOutDir.mkdir(parents=True, exist_ok=True)

        for i, ligand in enumerate(ligands, start=1):
            ligandPath = Path(ligand.getPoseFile())
            origName = ligandPath.name
            newName = f"{self._defaultName}_{origName}"

            dest = ligandOutDir / newName
            shutil.copy(ligandPath, dest)

    def renameFilesStep(self):
        extraPath = Path(self._getExtraPath())
        pdbid = self.getPDBId()

        proteinDir = extraPath / "protein"
        moleculeDir = extraPath / "molecule"
        ligandDir = moleculeDir / self._defaultName

        for file in proteinDir.iterdir():
            if file.name.startswith(f"{self._defaultName}_protein"):
                newProteinName = f"{pdbid}_protein{file.suffix}"
                file.rename(proteinDir / newProteinName)

        newLigandFolder = moleculeDir / pdbid
        if ligandDir.exists():
            ligandDir.rename(newLigandFolder)
            ligandDir = newLigandFolder

        for file in ligandDir.iterdir():
            if file.name.startswith(f"{self._defaultName}_"):
                newName = file.name.replace(
                    f"{self._defaultName}_",
                    f"{pdbid}_",
                    1
                )
                file.rename(ligandDir / newName)

        resultsDir = extraPath / "results"
        oldPrefix = f"{self._defaultName}_"
        newPrefix = f"{pdbid}_"

        if resultsDir.exists():
            for path in resultsDir.rglob("*"):
                if path.is_file() and path.name.startswith(oldPrefix):
                    newName = path.name.replace(oldPrefix, newPrefix, 1)
                    path.rename(path.with_name(newName))

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["The results have been saved in the extra folder as scorch2.results.tsv.", "Features and normalized features have been saved in the extra/results folder, along with the aggregated results and target results."]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputPDBligandFiles.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getPDBId(self):
        protein = self.inputPDBligandFiles.get().getProteinFile()
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
        for file in fileList:
            suffix = file.suffix.lower()
            basename = file.stem

            if suffix == ".pdbqt":
                continue

            pdbqtFile = Path(baseDir) / f"{basename}.pdbqt"
            outputPath = str(pdbqtFile.resolve())
            inputPath = str(file.resolve())

            if suffix == ".cif":
                pdbFile = Path(baseDir) / f"{basename}.pdb"
                cifToPdb(inputPath, str(pdbFile.resolve()))
                inputPath = str(pdbFile.resolve())

            args = f"-ipdb {inputPath} -opdbqt -O {outputPath}"
            runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())


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
        extraPath = Path(self._getExtraPath())
        resultsTsv = extraPath / "scorch2_results.tsv"

        scoreDict = {}

        with open(resultsTsv, "r", encoding="utf-8-sig") as f:
            lines = [line for line in f.readlines() if not line.startswith("#") and line.strip()]

        reader = csv.DictReader(lines, delimiter=",")

        for row in reader:
            compounId = row["compound_id"].strip()
            molName = compounId.split("_", 1)[1].replace(".pdbqt", "")
            score = float(row["sc2_score"])
            scoreDict[molName] = score
        return scoreDict

