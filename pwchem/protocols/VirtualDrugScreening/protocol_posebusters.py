# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
This protocol is used to score docking positions obtained by several software using descriptors and score functions
available in the Open Drug Discovery Toolkit (ODDT, https://github.com/oddt/oddt)

"""
import os, glob

from pyworkflow.protocol import params
from pyworkflow.utils import Message, createLink
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC

from pwem.protocols import EMProtocol
from pyworkflow.object import String

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *
from pwem.convert import cifToPdb


class ProtocolPoseBusters(EMProtocol):
    """
    Performs plausibility checks for generated molecule poses.
    Distance geometry: Use RDKit distance geometry bounds to check the geometry of a molecule.
    Energy ratio: Check whether the internal energy of a molecular conformation is too far from its ground state.
    Flatness: Check whether substructures of molecule are flat.
    Identity: Check two molecules are identical (docking relevant identity).
    Intermolecular distance: Check that predicted molecule is not too close and not too far away from conditioning molecule.
    RMSD: Calculate RMSD and related metrics between predicted molecule and closest ground truth molecule.
    Volume overlap: Check volume overlap between ligand and protein.
    More info about each test: https://posebusters.readthedocs.io/en/latest/api.html
    """
    _label = 'PoseBusters docking tests'
    _tests = [
        'All (run all compatible tests)',
        'Distance geometry (ligand only)',
        'Energy ratio (ligand only)',
        'Flatness (ligand only)',
        'Identity (requires crystal ligand)',
        'Intermolecular distance (requires protein)',
        'RMSD (requires crystal ligand)',
        'Volume overlap (requires protein + crystal ligand)'
    ]

    _scripts = os.path.abspath(os.path.join(POSEB_DIC['home'], 'posebusters/modules'))

    # -------------------------- DEFINE param functions ----------------------
    def _defineDistanceGeom(self, form):

        form.addParam('thrBadB', params.FloatParam, default=0.2, label="Bond length threshold: ",
                       help='Bonds may be up to x% longer than DG bounds.')
        form.addParam('thrClash', params.FloatParam, default=0.2, label="Overlap that constitutes a clash: ",
                      help='If set to 20%, two atoms may be up to 80% of the lower bound apart.')
        form.addParam('thrBadAngle', params.FloatParam, default=0.2, label="Bond angle threshold: ",
                      help='Bonds may be up to x% longer than DG bounds.')
        form.addParam('ignoreH', params.BooleanParam, default=True,
                        label="Ignore hydrogen: ",
                        help='Choose whether to ignore H.')
        form.addParam('sanitize', params.BooleanParam, default=True,
                      label="Sanitize molecule: ",
                      help='Choose whether to sanitize molecule before running.')
        form.addParam('symmetrize', params.BooleanParam, default=True,
                      label="Symmetrize conjugated terminal groups: ",
                      help='Will symmetrize the lower and upper bounds of the terminal conjugated bonds.')
        return form

    def _defineEnergyRatio(self, form):

        form.addParam('thrEnergyRatio', params.FloatParam, default=0.7, label="Energy ratio threshold: ",
                       help='Limit above which the energy ratio is deemed to high.')
        form.addParam('ensNumConf', params.FloatParam, default=100, label="Conformations: ",
                      help='Number of conformations to generate for the ensemble over which to average.')
        return form

    def _defineFlatness(self, form):

        form.addParam('thrFlatness', params.FloatParam, default=0.1, label="Flatness threshold: ",
                       help='Maximum distance from shared plane used as cutoff.')
        form.addParam('nonFlat', params.BooleanParam, default=False,
                        label="Check non-flat: ",
                        help='Whether to check the ring non-flatness instead of flatness. Turns (flatness <= threshold_flatness) to (flatness >= threshold_flatness).')
        return form

    def _defineIntermolDistance(self, form):

        form.addParam('radType', params.EnumParam, default=0, choices=['van der Waals', 'covalent'],
                      label="Type of atomic radius: ",
                       help='Type of atomic radius to use.')
        form.addParam('radScale', params.FloatParam, default=0.8, label="Radius scale: ",
                      help='Scaling factor for the atomic radii.')
        form.addParam('thrClash', params.FloatParam, default=0.05, label="Clash threshold: ",
                      help='Threshold for how much the atoms may overlap before a clash is reported.')
        form.addParam('ignoreTypes', params.StringParam, default='hydrogens',
                      label="Atoms to ignore in protein: ",
                       help="Which types of atoms to ignore in protein. May include: ['hydrogens', 'protein', 'organic cofactors', 'inorganic cofactors', 'waters']")
        form.addParam('maxDist', params.FloatParam, default=5.0, label="Maximum distance: ",
                      help='Maximum distance (in Angstrom) predicted and conditioning molecule may be apart to be considered as valid.')
        return form

    def _defineRMSD(self, form):

        form.addParam('thrRMSD', params.FloatParam, default=2.0, label="Threshold: ",
                      help='Threshold in angstrom for reporting whether RMSD is within threshold.')
        form.addParam('heavyOnly', params.BooleanParam, default=True,
                      label="Consider only heavy atoms: ",
                      help='Whether to only consider heavy atoms for RMSD calculation.')
        return form

    def _defineVolOverlap(self, form):

        form.addParam('thrClash', params.FloatParam, default=0.05, label="Overlap that constitutes a clash: ",
                      help='Threshold for how much volume overlap is allowed. This is the maximum share of volume of predicted molecule allowed to overlap with protein.')
        form.addParam('vdwScale', params.FloatParam, default=0.8, label="Scaling factor: ",
                      help='Scaling factor for the van der Waals radii which define the volume around each atom.')
        form.addParam('ignoreTypes', params.StringParam, default='hydrogens',
                      label="Atoms to ignore in protein: ",
                      help="Which types of atoms to ignore in protein. May include: ['hydrogens', 'protein', 'organic cofactors', 'inorganic cofactors', 'waters']")
        return form

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('tests', params.EnumParam, choices=self._tests, default=0,
                      label="Select test: ",
                      help=(
                          "Select which PoseBusters test to run.\n\n"
                          "? All: runs all tests compatible with the provided inputs\n"
                          "? Ligand-only tests require only a predicted ligand\n"
                          "? Some tests require a crystal ligand (True molecule)\n"
                          "? Some tests require a protein structure (Protein)"
                      ))

        form.addParam('oneFile', params.BooleanParam, default=True,
                        label="Test on one docked molecule: ",
                        help='Choose whether to run the test on one docked molecule or a set.')

        form.addParam('inputMoleculesSets', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Docked Small Molecules: ",
                      help='Select the docked molecules to be tested.')
        form.addParam('molPred', params.StringParam,
                        label='Predicted molecule: ', condition='oneFile',
                        help='Choose the predicted molecule (docked ligand).')
        form.addParam('inputMoleculesRefSets', params.PointerParam, condition='tests in [0, 4, 6]',
                      pointerClass='SetOfSmallMolecules', allowsNull=True,
                      label="Input Reference Small Molecules: ",
                      help='Select the reference molecules to be tested against.')
        form.addParam('molTrue', params.StringParam,
                      label='True molecule: ', condition='tests in [0, 4, 6]',
                      help='Choose the ground truth molecule (crystal ligand).')
        form.addParam('molCond', params.PointerParam, pointerClass="AtomStruct",
                      condition='tests in [0, 5, 7]', allowsNull=True,
                       label='Protein:', help='Choose the conditioning molecule (protein).')

        distGroup = form.addGroup('Distance geometry', condition='tests in [1]')
        self._defineDistanceGeom(distGroup)
        energyRatio = form.addGroup('Energy ratio', condition='tests in [2]')
        self._defineEnergyRatio(energyRatio)
        flatness = form.addGroup('Flatness', condition='tests in [3]')
        self._defineFlatness(flatness)
        intermolDist = form.addGroup('Intermolecular distance', condition='tests in [5]')
        self._defineIntermolDistance(intermolDist)
        rmsd = form.addGroup('RMSD', condition='tests in [6]')
        self._defineRMSD(rmsd)
        volOverlap = form.addGroup('Volume overlap', condition='tests in [7]')
        self._defineVolOverlap(volOverlap)

        group = form.addGroup('Scoring function', condition='tests in [0]')
        group.addParam('fullReport', params.BooleanParam, default=True, condition='tests in [0]',
                      label="Output full report: ",
                      help='Choose whether to output full report or not.')
        group.addParam('outputFormat', params.EnumParam, choices=['short', 'long', 'csv'], default=2,
                       condition='tests in [0]',
                       label="Output format: ",
                       help='Choose whether to output full report or not.')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        if self.tests.get() == 0:
            self._insertFunctionStep(self.allTestsStep)
        else:
            if self.oneFile.get():
                self._insertFunctionStep(self.indivTestsStep)
            else:
                self._insertFunctionStep(self.indivTestsAllMolsStep)
        self._insertFunctionStep(self.createOutputStep)

    def allTestsStep(self):
        args = ['bust ']

        if self.oneFile.get():
            #docked molecule
            molPred = self.getSpecifiedMol('pred')
            inpFile = self.convertFormat(molPred)
            args.append(os.path.abspath(inpFile))

            if self.inputMoleculesRefSets.get() is not None:
                #true molecule
                molTrue = self.getSpecifiedMol('true')
                inpFile = self.convertFormat(molTrue, type='crystal')
                args.append(f'-l {os.path.abspath(inpFile)}')

            if self.molCond.get() is not None:
                #protein
                molCond = self.molCond.get()
                inpFile = self.convertFormat(molCond, type='AtomStruct')
                args.append(f'-p {os.path.abspath(inpFile)}')

        else:
            for dockedMol in self.inputMoleculesSets.get():
                inpFile = self.convertFormat(dockedMol)
                args.append(os.path.abspath(inpFile))

        outputFormat = self.getEnumText('outputFormat')
        args.append(f'--outfmt {outputFormat}')

        if (self.fullReport.get()):
            args.append('--full-report')

        if (self.outputFormat.get() == 2):
            resultsFile = self._getPath('results.csv')
        else:
            resultsFile = self._getPath('results.txt')
        args.append(f'--output {os.path.abspath(resultsFile)}')

        Plugin.runCondaCommand(
            self,
            args=" ".join(args),
            condaDic=SCORCH2_DIC,
            program="",
            cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
        )

    def runPoseBustersForMol(self, molPred, molTrue=None, suffix=""):
        """
        Ejecuta PoseBusters para UNA molécula concreta
        """
        paramsFile = self._getExtraPath(f'testParams_{suffix}.txt')
        outputFolder = self._getPath(f'posebusters_results/{suffix}')

        inpFile = self.convertFormat(molPred)

        with open(paramsFile, 'w') as f:
            f.write(f"output = {os.path.abspath(outputFolder)}\n")
            f.write(f"mol_pred = {os.path.abspath(inpFile)}\n")

            if self.tests.get() == 1:
                f.write(f"test = distance_geometry\n")
                f.write(f"threshold_bad_bond_length = {self.thrBadB.get()}\n")
                f.write(f"threshold_clash = {self.thrClash.get()}\n")
                f.write(f"threshold_bad_angle = {self.thrBadAngle.get()}\n")
                f.write(f"ignore_hydrogens = {self.ignoreH.get()}\n")
                f.write(f"sanitize = {self.sanitize.get()}\n")
                f.write(f"symmetrize_conjugated_terminal_groups = {self.symmetrize.get()}\n")
            elif self.tests.get() == 2:
                f.write(f"test = check_energy_ratio\n")
                f.write(f"threshold_energy_ratio = {self.thrEnergyRatio.get()}\n")
                f.write(f"ensemble_number_conformations = {self.ensNumConf.get()}\n")
            elif self.tests.get() == 3:
                f.write(f"test = check_flatness\n")
                f.write(f"threshold_flatness = {self.thrFlatness.get()}\n")
                f.write(f"check_nonflat = {self.nonFlat.get()}\n")
            elif self.tests.get() == 4:
                f.write(f"test = check_identity\n")
                #molTrue = self.getSpecifiedMol('true')
                inpFile = self.convertFormat(molTrue, type='crystal')
                f.write(f"mol_true = {os.path.abspath(inpFile)}\n")
            elif self.tests.get() == 5:
                f.write(f"test = check_intermolecular_distance\n")
                inpFile = self.convertFormat(self.molCond.get(), type='AtomStruct')
                f.write(f"mol_cond = {inpFile}\n")
                if self.radType.get() == 0:
                    radType = 'vdw'
                else:
                    radType = self.getEnumText(self.radType.get())

                f.write(f"radius_type = {radType}\n")
                f.write(f"radius_scale = {self.radScale.get()}\n")
                f.write(f"clash_cutoff = {self.thrClash.get()}\n")
                types = self.ignoreTypes.get().replace(' ', '_').split(',')
                types = self.ignoreTypes.get().replace('[', '')
                f.write(f"ignore_types = {types}\n")
                f.write(f"max_distance = {self.maxDist.get()}\n")
            elif self.tests.get() == 6:
                f.write(f"test = check_rmsd\n")
                #molTrue = self.getSpecifiedMol('true')
                inpFile = self.convertFormat(molTrue, type='crystal')
                f.write(f"mol_true = {os.path.abspath(inpFile)}\n")
                f.write(f"rmsd_threshold = {self.thrRMSD.get()}\n")
                f.write(f"heavy_only = {self.heavyOnly.get()}\n")
            elif self.tests.get() == 7:
                f.write(f"test = check_volume_overlap\n")
                inpFile = self.convertFormat(self.molCond.get(), type='AtomStruct')
                f.write(f"mol_cond = {inpFile}\n")
                f.write(f"clash_cutoff = {self.thrClash.get()}\n")
                f.write(f"vdw_scale = {self.vdwScale.get()}\n")
                types = self.ignoreTypes.get().replace(' ', '_').split(',')
                types = self.ignoreTypes.get().replace('[', '')
                f.write(f"ignore_types = {types}\n")

        Plugin.runScript(
            self,
            'posebustersTesting.py',
            args=os.path.abspath(paramsFile),
            env=SCORCH2_DIC,
            cwd=self._getPath())


    def indivTestsStep(self):
        molPred = self.getSpecifiedMol('pred')
        molName = os.path.splitext(os.path.basename(molPred.getPoseFile()))[0]
        molTrue = None
        if self.molTrue.get() and self.inputMoleculesRefSets.get():
            for mol in self.inputMoleculesRefSets.get():
                if mol.__str__() == self.molTrue.get():
                    molTrue = mol.clone()
                    break
        self.runPoseBustersForMol(
            molPred,
            molTrue=molTrue,
            suffix=f'{molName}'
        )

    def indivTestsAllMolsStep(self):
        molTrue = None
        if self.molTrue.get() and self.inputMoleculesRefSets.get():
            for mol in self.inputMoleculesRefSets.get():
                if mol.__str__() == self.molTrue.get():
                    molTrue = mol.clone()
                    break

        for mol in self.inputMoleculesSets.get():
            molName = os.path.splitext(os.path.basename(mol.getPoseFile()))[0]
            self.runPoseBustersForMol(
                mol,
                molTrue=molTrue,
                suffix=f'{molName}'
            )


    def createOutputStep(self):
        if self.tests.get() == 0:
            if (self.outputFormat.get() == 2):
                resultsFile = self._getPath('results.csv')
            else:
                resultsFile = self._getPath('results.txt')
        else:
            resultsFileGeneral = self._getPath('posebusters_results')

        newMols = SetOfSmallMolecules.createCopy(self.inputMoleculesSets.get(), self._getPath(), copyInfo=True)

        predPose = self.getSpecifiedMol('pred')
        predPoseFile = predPose.getPoseFile() if predPose else None
        for mol in self.inputMoleculesSets.get():
            newMol = mol.clone()
            molName = os.path.splitext(os.path.basename(newMol.getPoseFile()))[0]
            if self.tests.get() != 0:
                resultsFile = os.path.join(resultsFileGeneral, f'{molName}')
            newMol.PoseBusters_file = String()
            if self.oneFile.get():
                if newMol.getPoseFile() == predPoseFile:
                    newMol.setAttributeValue('PoseBusters_file', resultsFile)
            else:
                newMol.setAttributeValue('PoseBusters_file', resultsFile)
            newMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputMoleculesSets.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def convertFormat(self, molPred, type=''):
        if 'AtomStruct' in type:
            basename = os.path.basename(molPred.getFileName()).split('.')[0]
            file = molPred.getFileName()
        elif 'crystal' in type:
            basename = os.path.basename(molPred.getFileName()).split('.')[0]
            file = molPred.getFileName()
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

    def getSpecifiedMol(self, string, one=False):
        myMol = None
        if string == 'pred':
            for mol in self.inputMoleculesSets.get():
                if mol.__str__() == self.molPred.get():
                    myMol = mol.clone()
                    break
        else :
            for mol in self.inputMoleculesRefSets.get():
                if mol.__str__() == self.molTrue.get():
                    myMol = mol.clone()
                    break

        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            return myMol
