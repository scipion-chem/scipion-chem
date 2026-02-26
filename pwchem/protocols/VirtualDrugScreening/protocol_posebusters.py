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

import csv, re

from pyworkflow.object import Boolean
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *
from pwem.convert import cifToPdb


class ProtocolPoseBusters(EMProtocol):
    """
        Performs plausibility checks for generated molecule poses.
        More info about the tests: https://posebusters.readthedocs.io/en/latest/cli.html
        """
    _label = 'PoseBusters docking tests'

    NTESTS = {'simple': 12, 'receptor': 22, 'trueMol': 18, 'trueMolAndReceptor': 28}

    DEF_COL_THRES = {'energy_ratio': 2.0, 'number_clashes': 0, 'aromatic_ring_maximum_distance_from_plane': 0.05,
                     'kabsch_rmsd': 2.5, 'centroid_distance': 12.0,
                     'volume_overlap_protein': 0.05, 'smallest_distance_protein': 6.0}

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

        fGroup = form.addGroup('Tests')
        fGroup.addParam('testsPassed', params.IntParam, label='Min. tests passed: ', default=self.NTESTS['receptor'],
                        help='Minimum number of tests passed to keep each molecule.')
        fGroup.addParam('doSpecific', params.BooleanParam, label='Check only specific tests: ', default='',
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Check only a subset of specific tests on filtering')
        fGroup.addParam('specificTests', params.TextParam, width=70, height=15, expertLevel=params.LEVEL_ADVANCED,
                        default='', label='List of tests to consider: ', condition='doSpecific',
                        help='List of tests to consider in the filtering.\nRemember to modify the "Min. tests passed" '
                             'param accordingly.')
        
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.poseBustersStep)
        self._insertFunctionStep(self.createOutputStep)

    def poseBustersStep(self):
        resultsFile = self.getResultsFile()

        baseArgs = []
        baseArgs.append('--outfmt csv')
        baseArgs.append(f'--output {os.path.abspath(resultsFile)}')
        baseArgs.append('--full-report')

        args = ['bust']
        for dockedMol in self.inputMoleculesSet.get():
            inpFile = self.convertFormat(dockedMol)
            args.append(os.path.abspath(inpFile))
        
        args.extend(self.getPoseBusterArgs())
        args.extend(baseArgs)

        Plugin.runCondaCommand(
            self, args=" ".join(args), condaDic=SCORCH2_DIC, program="",
            cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
        )

    def createOutputStep(self):
        resultsFile = self.getResultsFile()
        csvRows = self.getFileInfo(resultsFile)
        resTestsDic = self.parseResultRows(csvRows)

        newMols = SetOfSmallMolecules.createCopy(self.inputMoleculesSet.get(), self._getPath(), copyInfo=True)
        for mol in self.inputMoleculesSet.get():
            newMol = mol.clone()
            newMol.PoseBusters_file = String(resultsFile)

            resDic = resTestsDic[os.path.abspath(mol.getPoseFile())]
            if self.checkPassedTests(resDic) >= self.testsPassed.get():
                for testName, testVal in resDic.items():
                    newMol.__setattr__(testName, Boolean(testVal.upper() == 'TRUE'))
                newMols.append(newMol)

        self._defineOutputs(outputSmallMolecules=newMols)

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

    def getPoseBusterArgs(self):
        args = []

        if self.useTrueMol.get():
            molTrue = self.getTrueMol()
            trueFile = self.convertFormat(molTrue, type='crystal')
            args.append(f'-l {os.path.abspath(trueFile)}')

        if self.molRec.get():
            protFile = self.convertFormat(
                self.inputMoleculesSet.get().getProteinFile(), type='file'
            )
            args.append(f'-p {os.path.abspath(protFile)}')

        return args

    def getResultsFile(self):
        return self._getPath('results.csv')
    
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
    
    def parseResultRows(self, csvRows):
        resDic = {}
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
                csvRows[poseFile] = row

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

