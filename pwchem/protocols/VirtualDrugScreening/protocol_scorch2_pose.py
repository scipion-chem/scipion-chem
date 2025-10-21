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
from pathlib import Path

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.utils import os, shutil, re, runOpenBabel
from pwchem import Plugin, SCORCH2_DIC


currentDir = Path(__file__).parent.resolve()
scriptRescoring =  os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'SCORCH2/scorch2_rescoring.py'))


class ProtocolSCORCH2(EMProtocol):
    """
    Computes the best poses of protein-ligand interactions using SCORCH2: https://github.com/LinCompbio/SCORCH2
    """
    _label = 'SCORCH2 rescoring'

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
        # Ligand and protein file
        iGroup.addParam('inputPDBproteinFile', params.PointerParam, pointerClass='AtomStruct',
                        label='Input protein: ', allowsNull=True,
                        condition='not useFeatures',
                        help='Input PDB protein file.')
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

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        preExtracted = self.useFeatures.get()
        if not preExtracted:
            self._insertFunctionStep('moveFiles')

        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('createOutputStep')

    def convertInputStep(self):
        proteinDir = Path(self._getExtraPath()) / "protein"
        ligandDir = Path(self._getExtraPath()) / "molecule"

        proteinFiles = list(proteinDir.glob("*"))
        if proteinFiles:
            print("convert protein")
            self.convertFiles(proteinFiles, proteinDir)
            self.removePdbFiles(proteinDir)
        else:
            print("[WARNING] No protein files found.")

        ligandFiles = [f for f in ligandDir.rglob("*") if f.is_file()]
        if ligandFiles:
            print("convert ligands")
            self.convertFiles(ligandFiles, ligandDir)
            self.removePdbFiles(ligandDir)
        else:
            print("[WARNING] No ligand files found.")

    def createOutputStep(self):
        proteinDir = Path(self._getExtraPath()) / "protein"
        ligandDir = Path(self._getExtraPath()) / "molecule"

        projectDir = Path(self.getWorkingDir()).parent.resolve().parent
        proteinFull = projectDir / proteinDir
        ligandFull = projectDir / ligandDir

        outputFile = self._getExtraPath('scorch2_results.tsv')
        outputFileFull = projectDir / outputFile

        modelsDir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'scorchModels/models'))

        sc2PsModel = os.path.abspath(os.path.join(modelsDir, 'sc2_ps.xgb'))
        sc2PbModel = os.path.abspath(os.path.join(modelsDir, 'sc2_pb.xgb'))
        psScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_ps_scaler'))
        pbScaler = os.path.abspath(os.path.join(modelsDir, 'sc2_pb_scaler'))

        preExtracted = self.useFeatures.get()
        print(scriptRescoring)
        if not preExtracted:
            args = [
                str(scriptRescoring),
                "--protein-dir", str(proteinFull),
                "--ligand-dir", str(ligandFull),
                "--sc2_ps_model", str(sc2PsModel),
                "--sc2_pb_model", str(sc2PbModel),
                "--ps_scaler", str(psScaler),
                "--pb_scaler", str(pbScaler),
                "--output", str(outputFileFull),
                "--ps_weight", str(self.psWeight.get()),
                "--pb_weight", str(self.pbWeight.get()),
                "--keep-temp",
                "--res-dir", str('results') # why is this created in the repo and not in the gui?
            ]
        else:
            args = [
                str(scriptRescoring),
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
        #move results folder to scipion outputs, idk why it is created in the dir where the script is called
        if not preExtracted:
            resDir = Path(self._getExtraPath()) / "results"
            shutil.move(os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']),'results')), resDir)

    def moveFiles(self):
        extraPath = Path(self._getExtraPath())

        proteinDir = extraPath / "protein"
        moleculeDir = extraPath / "molecule"
        proteinDir.mkdir(parents=True, exist_ok=True)
        moleculeDir.mkdir(parents=True, exist_ok=True)

        protein = self.inputPDBproteinFile.get()
        print(protein)
        proteinPath = Path(protein.getFileName())
        pdbId = proteinPath.stem
        proteinFile = proteinDir / f"{pdbId}_protein{proteinPath.suffix}"
        shutil.copy(proteinPath, proteinFile)

        ligands = self.inputPDBligandFiles.get()
        ligandOutDir = moleculeDir / pdbId
        ligandOutDir.mkdir(parents=True, exist_ok=True)

        # Regex for valid ligand naming pattern
        validPattern = re.compile(rf"^{re.escape(pdbId)}_[A-Za-z0-9]+_pose.*\.pdbqt$")

        for i, ligand in enumerate(ligands, start=1):
            ligandPath = Path(ligand.getFileName())
            origName = ligandPath.name

            if validPattern.match(origName):
                newName = origName
            else:
                # Extract numeric ID from ligand filename
                match = re.match(r"(\d+)", ligandPath.stem)
                number = match.group(1) if match else str(i)

                # Generate new standardized name
                newName = f"{pdbId}_{number}_pose{i}{ligandPath.suffix}"

            dest = ligandOutDir / newName
            shutil.copy(ligandPath, dest)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ["The results have been saved in the extra folder as scorch2.results.tsv.", "Features and normalized features have been saved in the extra/results folder, along with the aggregated results and target results."]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------

    def checkPdbqtFiles(self, directory):
        """Check if files are PDBQT"""
        files = list(directory.glob("*"))
        for f in files:
            if f.suffix.lower() != ".pdbqt":
                return False, files
            else:
                return True, files

    def convertFiles(self, fileList, baseDir):
        """Convert PDB to PDBQT, keeping output in the same folder as the input file"""
        for file in fileList:
            if file.suffix.lower() == ".pdb":
                basename = file.stem
                pdbqtFile = Path(baseDir) / f"{basename}.pdbqt"

                inputAbs = str(file.resolve())
                outputAbs = str(pdbqtFile.resolve())

                args = f"-ipdb {inputAbs} -opdbqt -O {outputAbs}"
                runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
            else:
                pass

    def removePdbFiles(self, directory):
        """Removes .pdb files only if their corresponding .pdbqt exists in the same directory."""
        for pdbFile in directory.rglob("*.pdb"):
            pdbqtFile = pdbFile.with_suffix(".pdbqt")
            if pdbqtFile.exists():
                try:
                    pdbFile.unlink()
                except Exception as e:
                    print(f"[WARNING] Could not delete {pdbFile.name}: {e}")
