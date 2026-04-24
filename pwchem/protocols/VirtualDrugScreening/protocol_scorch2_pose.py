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
from Bio.PDB import PDBParser, PDBIO, Select

from pwem.convert import cifToPdb
from pwem.protocols import EMProtocol
from pyworkflow.object import Float
from pyworkflow.protocol import params

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import os, shutil, re, runOpenBabel, makeSubsets, insistentRun, pdbFromASFile
from pwchem import Plugin, SCORCH2_DIC


currentDir = Path(__file__).parent.resolve()


class CoordinateRangeSelect(Select):
    """
    Select residues that have at least one atom within the specified
    coordinate range, extended by N in all directions.
    """

    def __init__(self, coordinateRanges):
        (self.x_min, self.x_max), (self.y_min, self.y_max), (self.z_min, self.z_max) = coordinateRanges

    def accept_residue(self, residue):
        """
        Accept residue if any of its atoms is within the coordinate range.
        """
        # Check all atoms in the residue
        for atom in residue.get_atoms():
            coord = atom.get_coord()
            if (self.x_min <= coord[0] <= self.x_max and
                    self.y_min <= coord[1] <= self.y_max and
                    self.z_min <= coord[2] <= self.z_max):
                return 1  # Accept this residue
        return 0  # Reject this residue

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
        form.addParam('useGPU', params.BooleanParam, default=True, label="Use GPU: ",
                      help='Whether to use GPU or not. (Unable to choose the GPU id).')
        iGroup = form.addGroup('Input')
        # Pre-extracted features
        iGroup.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input ligand: ',
                        help='Input docked small molecules to rescore')
        iGroup.addParam('cropReceptor', params.BooleanParam, default=False, label="Crop receptor: ",
                        help='Crop receptor for each pocket (20A around positions) to accelerate the '
                             'feature calculation.')

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
        molDic = self.getMolsSetsDic()
        rStep = self._insertFunctionStep(self.organizeInputStep, molDic)

        sSteps = []
        for pockId, molList in molDic.items():
            for it in range(self.getNBatches(molList)):
                cStep = self._insertFunctionStep(self.convertInputStep, pockId, it, prerequisites=[rStep])
                sSteps += [self._insertFunctionStep(self.scorchStep, pockId, it, prerequisites=[cStep])]
        oStep = self._insertFunctionStep(self.createOutputStep, prerequisites=sSteps)

    def organizeInputStep(self, molDic):
        inProteinPath = Path(self.inputSmallMolecules.get().getProteinFile())

        # Ligands
        for pockId, molList in molDic.items():
            nBatches = self.getNBatches(molList)
            molSubSets = makeSubsets(molList, nBatches, True)
            for it, subset in enumerate(molSubSets):
                moleculeDir = self.getMoleculesDir(pockId, it)
                moleculeDir.mkdir(parents=True, exist_ok=True)

                ligandOutDir = moleculeDir / self._defaultName
                ligandOutDir.mkdir(parents=True, exist_ok=True)

                for i, ligand in enumerate(subset, start=1):
                    ligandPath = Path(ligand.getPoseFile())
                    origName = ligandPath.name
                    newName = f"{self._defaultName}_{origName}"

                    dest = ligandOutDir / newName
                    shutil.copy(ligandPath, dest)

            # Receptor
            proteinDir = self.getProtDir(pockId)
            proteinDir.mkdir(parents=True, exist_ok=True)

            proteinFile = proteinDir / f"{self._defaultName}_protein{inProteinPath.suffix}"

            tmpFile = os.path.join(os.path.dirname(oFile), 'tempPDB.pdb')
            tmpFile = pdbFromASFile(inProteinPath, tmpFile)

            if self.cropReceptor.get():
                self.cropProteinFile(tmpFile, proteinFile, molList)
            else:
                shutil.copy(tmpFile, proteinFile)

            os.remove(tmpFile)

            proteinFiles = list(proteinDir.glob("*"))
            if proteinFiles:
                self.convertFiles(proteinFiles, proteinDir)
                self.removePdbFiles(proteinDir)
            else:
                logging.warning("No protein files found.")


    def convertInputStep(self, pockId, it):
        ligandDir = self.getMoleculesDir(pockId, it)

        ligandFiles = []
        for subfolder in ligandDir.iterdir():
            if subfolder.is_dir():
                ligandFiles.extend(f for f in subfolder.glob("*") if f.is_file())
        if ligandFiles:
            self.convertFiles(ligandFiles, os.path.abspath(subfolder))
            self.removePdbFiles((subfolder))
        else:
            logging.warning("No ligand files found.")


    def scorchStep(self, pockId, it):
        oriProteinDir = self.getProtDir(pockId)
        proteinDir = self.getProtDir(pockId, it)
        shutil.copytree(oriProteinDir, proteinDir)

        ligandDir = self.getMoleculesDir(pockId, it)

        scriptRescoringDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'SCORCH2'))
        modelsDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'scorchModels/models'))

        script = os.path.abspath(os.path.join(scriptRescoringDir, 'scorch2_rescoring.py'))
        sc2PsModel = os.path.abspath(os.path.join(modelsDir, 'sc2_ps.xgb'))
        sc2PbModel = os.path.abspath(os.path.join(modelsDir, 'sc2_pb.xgb'))
        psScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_ps_scaler'))
        pbScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_pb_scaler'))

        resDir = self.getResultsDir(pockId, it)
        resDir.mkdir(exist_ok=True)
        outputFileFull = self.getBatchDir(pockId, it) / 'scorch2_results.tsv'
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
    def getPocketDir(self, it):
        return Path(self._getExtraPath()) / f"pocket_{it}"

    def getBatchDir(self, pockId, it):
        return self.getPocketDir(pockId) / f"batch_{it}"

    def getProtDir(self, pockId, it=None):
        if it is None:
            return self.getPocketDir(pockId) / "protein"
        return self.getBatchDir(pockId, it) / "protein"

    def getMoleculesDir(self, pockId, it):
        return self.getBatchDir(pockId, it) / "molecule"

    def getResultsDir(self, pockId, it):
        return self.getBatchDir(pockId, it) / f"results"

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

            args = f"-ipdb {inputPath} -opdbqt -O {outputPath}"
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
        molDic = self.getMolsSetsDic()
        for pockId, molList in molDic.items():
            for it in range(self.getNBatches(molList)):
                resultsTsv = self.getBatchDir(pockId, it) / "scorch2_results.tsv"
                with open(resultsTsv, "r", encoding="utf-8-sig") as f:
                    lines = [line for line in f.readlines() if not line.startswith("#") and line.strip()]

                reader = csv.DictReader(lines, delimiter=",")

                for row in reader:
                    compounId = row["compound_id"].strip()
                    molName = compounId.split("_", 1)[1].replace(".pdbqt", "")
                    score = float(row["sc2_score"])
                    scoreDict[molName] = score
        return scoreDict

    def getMolsSetsDic(self):
        '''Return a dictionary {pocketId: [molList]}'''
        inMols = self.inputSmallMolecules.get()
        if self.cropReceptor.get():
            molIdDic = inMols.updateLigandsDic({}, inMols, 'pocket')
            molDic = {pockId: inMols.getMolsFromIds(molIds) for pockId, molIds in molIdDic.items()}
        else:
            molDic = {0: inMols}
        return molDic

    def getNBatches(self, molList):
        nMols = len(molList)
        nBatches = (nMols // self.batchSize.get()) + 1
        return nBatches

    def getExtendedBounds(self, coordsDics, N):
        """
        Returns:
            Tuple of (min_x, max_x, min_y, max_y, min_z, max_z) each extended by N
        """
        # Extract all coordinates
        allCoords = [c for coordsDic in coordsDics for c in list(coordsDic.values())]

        # Get min and max for each dimension
        minX, maxX = min(coord[0] for coord in allCoords) - N, max(coord[0] for coord in allCoords) + N
        minY, maxY = min(coord[1] for coord in allCoords) - N, max(coord[1] for coord in allCoords) + N
        minZ, maxZ = min(coord[2] for coord in allCoords) - N, max(coord[2] for coord in allCoords) + N

        return [(minX, maxX), (minY, maxY), (minZ, maxZ)]

    def cropProteinFile(self, inFile, oFile, molList):
        atomsPosDics = [mol.getAtomsPosDic() for mol in molList]
        limitCoords = self.getExtendedBounds(atomsPosDics, 20)

        selector = CoordinateRangeSelect(limitCoords)
        structModel = PDBParser().get_structure('receptor', inFile)[0]

        # Write the filtered structure
        io = PDBIO()
        io.set_structure(structModel)
        io.save(str(oFile), selector)

        return oFile

