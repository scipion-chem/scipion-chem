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

from pwem.objects import AtomStruct
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from Bio.PDB import MMCIFParser, PDBIO
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import subprocess

from pwchem import Plugin, SCORCH2_DIC
from pwchem.utils import fillEmptyAttributes

current_dir = Path(__file__).parent.resolve()
scriptRescoring = current_dir.parent.parent/'scripts'/'scorch2_rescoring.py'


class ProtocolSCORCH2(EMProtocol):
    """
    Computes the best poses of protein-ligand interactions
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
                        label='Input PDB protein file: ', allowsNull=True,
                        condition='not useFeatures',
                        help='Input PDB protein file.')
        iGroup.addParam('inputPDBligandFiles', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input PDB ligand files: ', allowsNull=True,
                        condition='not useFeatures',
                        help='Input folder with PDB ligand files.')
        # Aggregate
        iGroup.addParam('aggregate', params.BooleanParam, default=False,
                        label="Aggregate results: ",
                        #condition='not useFeatures',
                        help='If Yes, the best pose will be selected per compound with aggregation metadata. If No, all poses will be scored individually and ranked by SC2 score.')

        iGroup.addParam('ps_weight', params.FloatParam, default=0.7, expertLevel=params.LEVEL_ADVANCED,
                        label='PS weight: ', help="Weight for SC2-PS predictions")
        iGroup.addParam('pb_weight', params.FloatParam, default=0.3, expertLevel=params.LEVEL_ADVANCED,
                        label='PB weight: ', help="Weight for SC2-PB predictions")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        pre_extracted = self.useFeatures.get()
        if not pre_extracted:
            self._insertFunctionStep('prepareFiles')

        self._insertFunctionStep('createOutputStep')


    def createOutputStep(self):
        protein_dir = Path(self._getExtraPath()) / "protein"
        ligand_dir = Path(self._getExtraPath()) / "molecule"

        project_dir = Path(self.getWorkingDir()).parent.resolve().parent
        protein_full = project_dir / protein_dir
        ligand_full = project_dir / ligand_dir

        output_file = self._getExtraPath('scorch2_results.tsv')
        output_file_full = project_dir / output_file

        models_dir = os.path.abspath(os.path.join(Plugin.getVar(SCORCH2_DIC['home']), 'scorchModels/SCORCH2_models/models'))

        sc2_ps_model = os.path.abspath(os.path.join(models_dir, 'sc2_ps.xgb'))
        sc2_pb_model = os.path.abspath(os.path.join(models_dir, 'sc2_pb.xgb'))
        ps_scaler = os.path.abspath(os.path.join(models_dir, 'sc2_ps_scaler'))
        pb_scaler = os.path.abspath(os.path.join(models_dir, 'sc2_pb_scaler'))

        pre_extracted = self.useFeatures.get()
        if not pre_extracted:
            args = [
                str(scriptRescoring),
                "--protein-dir", str(protein_full),
                "--ligand-dir", str(ligand_full),
                "--sc2_ps_model", str(sc2_ps_model),
                "--sc2_pb_model", str(sc2_pb_model),
                "--ps_scaler", str(ps_scaler),
                "--pb_scaler", str(pb_scaler),
                "--output", str(output_file_full),
                "--ps_weight", str(self.ps_weight.get()),
                "--pb_weight", str(self.pb_weight.get()),
                "--keep-temp",
                "--res-dir", str('results') # why is this created in the repo and not in the gui?
            ]
        else:
            args = [
                str(scriptRescoring),
                "--features", str(self.inputFeatures.get()),
                "--sc2_ps_model", str(sc2_ps_model),
                "--sc2_pb_model", str(sc2_pb_model),
                "--ps_scaler", str(ps_scaler),
                "--pb_scaler", str(pb_scaler),
                "--output", str(output_file_full),
                "--ps_weight", str(self.ps_weight.get()),
                "--pb_weight", str(self.pb_weight.get())
            ]

        if self.aggregate.get():
            args.append("--aggregate")

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=SCORCH2_DIC,
            program="python",
            cwd=str(current_dir.parent.parent)
        )
        #move results folder to scipion outputs, idk why it is created in the dir where the script is called
        if not pre_extracted:
            res_dir = Path(self._getExtraPath()) / "results"
            shutil.move(current_dir.parent.parent/'results', res_dir)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        pass

    def _warnings(self):
        pass

    # --------------------------- UTILS functions -----------------------------------

    def prepareFiles(self):
        extraPath = Path(self._getExtraPath())

        proteinDir = extraPath / "protein"
        moleculeDir = extraPath / "molecule"
        proteinDir.mkdir(parents=True, exist_ok=True)
        moleculeDir.mkdir(parents=True, exist_ok=True)

        protein = self.inputPDBproteinFile.get()
        protein_path = Path(protein.getFileName())
        pdb_id = protein_path.stem
        protein_file = proteinDir / f"{pdb_id}_protein.pdbqt"
        shutil.copy(protein_path, protein_file)

        ligands: SetOfSmallMolecules = self.inputPDBligandFiles.get()
        ligandOutDir = moleculeDir / pdb_id
        ligandOutDir.mkdir(parents=True, exist_ok=True)

        for ligand in ligands:
            ligand_path = Path(ligand.getFileName())
            dest = ligandOutDir / ligand_path.name
            shutil.copy(ligand_path, dest)