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
This protocol is used to perform tests oer predicted molecule poses (https://github.com/maabuu/posebusters).

"""
import os.path

import csv, re

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwchem import Plugin, POSEB_DIC, SCORCH2_DIC

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwem.convert import cifToPdb


class ProtocolPoseBusters(EMProtocol):
    FILTER_COLUMNS = {
        "normal": ['energy_ratio', 'number_clashes', 'aromatic_ring_maximum_distance_from_plane'],
        "true": ['energy_ratio', 'number_clashes', 'aromatic_ring_maximum_distance_from_plane',
                 'kabsch_rmsd', 'centroid_distance'],
        "prot": ['energy_ratio', 'number_clashes', 'aromatic_ring_maximum_distance_from_plane',
                 'volume_overlap_protein', 'smallest_distance_protein'],
        "true_prot": ['energy_ratio', 'number_clashes', 'aromatic_ring_maximum_distance_from_plane',
                      'kabsch_rmsd', 'centroid_distance',
                      'volume_overlap_protein', 'smallest_distance_protein']
    }

    normalTests = 12
    trueMolTests = 18
    trueMolProtTests = 28
    protTests = 22
    """
    Performs plausibility checks for generated molecule poses.
    More info about the tests: https://posebusters.readthedocs.io/en/latest/cli.html
    """
    _label = 'PoseBusters docking tests'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
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
        form.addParam('useTrueMol', params.BooleanParam, default=False,
                      label="Use true molecule: ",
                      help='Choose whether to use the true molecule in the analysis.')
        form.addParam('inputMoleculesRefSets', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=True,
                      label="Input Reference Small Molecules: ", condition='useTrueMol',
                      help='Select the reference molecules to be tested against.')
        form.addParam('molTrue', params.StringParam,
                      label='True molecule: ', condition='useTrueMol',
                      help='Choose the ground truth molecule (crystal ligand).')
        form.addParam('molCond', params.BooleanParam, default=False,
                       label="Use protein to run tests: ",
                       help='Choose whether to output use protein to run tests. (If True, it will be used to run all default tests that require the protein as input.)')


        form.addParam('outputFormat', params.EnumParam, choices=['short', 'long', 'csv'], default=2,
                       label="Output format: ",
                       help='Choose whether to output full report or not.')
        form.addParam('fullReport', params.BooleanParam, default=True, condition='outputFormat != 0',
                      label="Output full report: ",
                      help='Choose whether to output full report or not.')
        form.addParam('filter', params.BooleanParam, default=False,
                      label="Filter molecules: ",
                      help='Choose whether to filter molecules based on tests.')

        form = form.addGroup('Filter', condition='filter')
        form.addParam('chooseFilterLongTxt', params.EnumParam, default=0, condition='outputFormat != 0 and fullReport',
                      label="Choose filtering option: ", choices=['tests passed', 'threshold over numerical result'],
                      help='Choose how to apply the filter over the results.')

        form.addParam('testsPassedNormal', params.IntParam, condition='((chooseFilterLongTxt == 0) or (outputFormat == 0 and not fullReport)) and not useTrueMol and not molCond',
                      label='Tests passed: ', default=self.normalTests,
                      help='Number of tests passed to keep the molecules.')
        form.addParam('testsPassedTrueMol', params.IntParam,
                      condition='((chooseFilterLongTxt == 0) or (outputFormat == 0 and not fullReport)) and useTrueMol and not molCond',
                      label='Tests passed: ', default=self.trueMolTests,
                      help='Number of tests passed to keep the molecules.')
        form.addParam('testsPassedTrueMolProt', params.IntParam,
                      condition='((chooseFilterLongTxt == 0) or (outputFormat == 0 and not fullReport)) and useTrueMol and molCond',
                      label='Tests passed: ', default=self.trueMolProtTests,
                      help='Number of tests passed to keep the molecules.')
        form.addParam('testsPassedProt', params.IntParam,
                      condition='((chooseFilterLongTxt == 0) or (outputFormat == 0 and not fullReport)) and not useTrueMol and molCond',
                      label='Tests passed: ', default=self.protTests,
                      help='Number of tests passed to keep the molecules.')

        form.addParam('filterColLongTxt', params.EnumParam, default=0,
                      condition=' chooseFilterLongTxt == 1 and not molCond and not useTrueMol',
                      label="Column: ", choices=self.FILTER_COLUMNS["normal"],
                      help='Filtering columns options.')
        form.addParam('filterColLongTxtTrueMol', params.EnumParam, default=0,
                      condition='chooseFilterLongTxt == 1 and not molCond and useTrueMol',
                      label="Column: ", choices=self.FILTER_COLUMNS["true"],
                      help='Filtering columns options.')
        form.addParam('filterColLongTxtProt', params.EnumParam, default=0,
                      condition='chooseFilterLongTxt == 1 and molCond and not useTrueMol',
                      label="Column: ", choices=self.FILTER_COLUMNS["prot"],
                      help='Filtering columns options.')
        form.addParam('filterColLongTxtAll', params.EnumParam, default=0,
                      condition='chooseFilterLongTxt == 1 and molCond and useTrueMol',
                      label="Column: ", choices=self.FILTER_COLUMNS["true_prot"],
                      help='Filtering columns options.')

        form.addParam('filterOp', params.EnumParam, label='Filter operation: ',
                       condition='chooseFilterLongTxt==1', default=3,
                       choices=['==', '>', '>=', '<', '<=', 'between'])

        form.addParam('energyRatio', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and ((not molCond and not useTrueMol and filterColLongTxt == 0) or (not molCond and useTrueMol and filterColLongTxtTrueMol == 0) or (molCond and not useTrueMol and filterColLongTxtProt == 0) or \
                              (molCond and useTrueMol and filterColLongTxtAll == 0))',
                      label='Energy ratio threshold: ', default=2.0,
                      help='Threshold to keep molecules with energy ratio below it. Acts as a lower threshold when filtering "between".')
        form.addParam('numClashes', params.IntParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxt == 1 or filterColLongTxtTrueMol==1 or filterColLongTxtProt==1 or filterColLongTxtAll==1)',
                      label='Number of clashes threshold: ', default=0,
                      help='Threshold to keep molecules with number or clashes below it. Acts as a lower threshold when filtering "between".')
        form.addParam('ringPlanarity', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxt == 2 or filterColLongTxtTrueMol==2 or filterColLongTxtProt==2 or filterColLongTxtAll==2)',
                      label='Ring planarity threshold: ', default=0.05,
                      help='Threshold to keep molecules with ring planarity below it. Acts as a lower threshold when filtering "between".')
        form.addParam('rmsd', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxtTrueMol==3) or filterColLongTxtAll==3',
                      label='RMSD threshold: ', default=2.5,
                      help='Threshold to keep molecules with rmsd below it. Acts as a lower threshold when filtering "between".')
        form.addParam('centroidDist', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxtTrueMol==4) or filterColLongTxtAll==4',
                      label='Centroid distance threshold: ', default=12.0,
                      help='Threshold to keep molecules with centroid distance below it. Acts as a lower threshold when filtering "between".')
        form.addParam('volOverlap', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxtProt==3) or filterColLongTxtAll==5',
                      label='Volume overlap threshold: ', default=0.05,
                      help='Threshold to keep molecules with volume overlap below it. Acts as a lower threshold when filtering "between".')
        form.addParam('smallDist', params.FloatParam,
                      condition='chooseFilterLongTxt == 1 and (filterColLongTxtProt==4) or filterColLongTxtAll==5',
                      label='Smallest distance threshold: ', default=6.0,
                      help='Threshold to keep molecules with smallest distance to protein below it. Acts as a lower threshold when filtering "between".')

        form.addParam('upper', params.FloatParam,
                      condition='filterOp==5',
                      label='Upper threshold: ', default='0',
                      help='Upper threshold to keep molecules.')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.poseBustersStep)
        self._insertFunctionStep(self.createOutputStep)

    def poseBustersStep(self):
        outputFormat = self.getEnumText('outputFormat')
        resultsFile = self.getResultsFile()

        baseArgs = []
        baseArgs.append(f'--outfmt {outputFormat}')
        baseArgs.append(f'--output {os.path.abspath(resultsFile)}')

        if self.fullReport.get():
            baseArgs.append('--full-report')

        # Case 1: oneFile - always single run
        if self.oneFile.get():
            args = self.case1Args()

            args.extend(baseArgs)

            Plugin.runCondaCommand(
                self,
                args=" ".join(args),
                condaDic=SCORCH2_DIC,
                program="",
                cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
            )

        # Case 2: multiple molecules
        else:
            if self.useTrueMol.get() or self.molCond.get():

                for dockedMol in self.inputMoleculesSets.get():
                    args , resultFiles = self.case2Args(dockedMol)

                    Plugin.runCondaCommand(
                        self,
                        args=" ".join(args),
                        condaDic=SCORCH2_DIC,
                        program="",
                        cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
                    )

                self.mergeResultFiles(resultFiles)

            # Case 3: no -l, no -p - batch allowed
            else:
                args = ['bust']

                for dockedMol in self.inputMoleculesSets.get():
                    inpFile = self.convertFormat(dockedMol)
                    args.append(os.path.abspath(inpFile))

                args.extend(baseArgs)

                Plugin.runCondaCommand(
                    self,
                    args=" ".join(args),
                    condaDic=SCORCH2_DIC,
                    program="",
                    cwd=os.path.abspath(Plugin.getVar(POSEB_DIC['home']))
                )

    def createOutputStep(self):
        resultsFile = self.getResultsFile()

        newMols = SetOfSmallMolecules.createCopy(self.inputMoleculesSets.get(), self._getPath(), copyInfo=True)

        csvRows, txtRows = self.getFileInfo(resultsFile)

        predPose = self.getSpecifiedMol('pred')
        predPoseFile = predPose.getPoseFile() if predPose else None

        for mol in self.inputMoleculesSets.get():
            keep = True
            add = True

            newMol = mol.clone()
            newMol.PoseBusters_file = String()
            if self.oneFile.get():
                if newMol.getPoseFile() == predPoseFile:
                    newMol.setAttributeValue('PoseBusters_file', resultsFile)
                else:
                    add = False
            else:
                newMol.setAttributeValue('PoseBusters_file', resultsFile)

            #filter
            if self.filter.get():
                keep = self.filterMolecules(resultsFile, csvRows, txtRows, mol)

            if keep and add:
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
        molSet = self.inputMoleculesSets.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        tests, testNum = self.getTestInfo()

        if tests > testNum or tests < 1:
            validations += [f'only {testNum} tests are performed.']

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def case1Args(self):
        args = ['bust']

        molPred = self.getSpecifiedMol('pred')
        inpFile = self.convertFormat(molPred)
        args.append(os.path.abspath(inpFile))

        if self.useTrueMol.get():
            molTrue = self.getSpecifiedMol('true')
            trueFile = self.convertFormat(molTrue, type='crystal')
            args.append(f'-l {os.path.abspath(trueFile)}')

        if self.molCond.get():
            protFile = self.convertFormat(
                self.inputMoleculesSets.get().getProteinFile(),
                type='file'
            )
            args.append(f'-p {os.path.abspath(protFile)}')
        return args

    def case2Args(self, dockedMol):
        args = ['bust']
        resultFiles = []
        outputFormat = self.getEnumText('outputFormat')

        inpFile = self.convertFormat(dockedMol)
        args.append(os.path.abspath(inpFile))

        if self.useTrueMol.get():
            molTrue = self.getSpecifiedMol('true')
            trueFile = self.convertFormat(molTrue, type='crystal')
            args.append(f'-l {os.path.abspath(trueFile)}')

        if self.molCond.get():
            protFile = self.convertFormat(
                self.inputMoleculesSets.get().getProteinFile(),
                type='file'
            )
            args.append(f'-p {os.path.abspath(protFile)}')

        resultsFile = self.getResultFileForMol(dockedMol)
        resultFiles.append(resultsFile)

        args.append(f'--outfmt {outputFormat}')
        args.append(f'--output {os.path.abspath(resultsFile)}')

        if self.fullReport.get():
            args.append('--full-report')
        return args, resultFiles

    def getResultsFile(self):
        return self._getPath(
            'results.csv' if self.outputFormat.get() == 2 else 'results.txt'
        )

    def filterMolecules(self, resultsFile, csvRows, txtRows, mol):
        poseName = os.path.basename(mol.getPoseFile())

        keep = True
        if resultsFile.endswith('.csv'):
            if not self.fullReport.get() or (self.fullReport.get() and self.chooseFilterLongTxt.get() == 0):
                keep = self.rowPassesShortCsvTests(csvRows[poseName])
            else:
                keep = self.rowPassesColValue(csvRows[poseName])
        else:
            if self.outputFormat.get() == 0:  # short txt
                poseLine = None
                poseBase = os.path.splitext(poseName)[0]
                for key, line in txtRows.items():
                    txtBase = os.path.splitext(os.path.basename(line.split()[1]))[0]
                    if poseBase == txtBase:
                        poseLine = line
                        break
                if poseLine is not None:
                    keep = self.molPassesShortTxtTests(poseLine)
            else:  # long txt
                if self.chooseFilterLongTxt.get() == 1:  # filter value w/ threshold
                    poseName = os.path.basename(mol.getPoseFile()).split('.')[0]
                    poseLine = txtRows.get(poseName)
                    if poseLine is not None:
                        keep = self.rowPassesColValueTxt(poseLine)
                else:
                    poseName = os.path.basename(mol.getPoseFile()).split('.')[0]
                    poseLine = txtRows.get(poseName)
                    if poseLine is not None:
                        keep = self.rowPassesColTxt(poseLine)
        return keep

    def getFileInfo(self, resultsFile):
        csvRows = {}
        txtRows = {}
        if resultsFile.endswith('.csv'):
            with open(resultsFile) as f:
                reader = csv.DictReader(f, delimiter=',')
                for row in reader:
                    poseName = os.path.basename(row['molecule'])
                    csvRows[poseName] = row
        else:
            if self.outputFormat.get() == 0:
                with open(resultsFile) as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        poseFile = line.split()[0]
                        poseName = os.path.basename(poseFile)
                        txtRows[poseName] = line
            else:
                with open(resultsFile) as f:
                    currentPose = None
                    buffer = []
                    for line in f:
                        line = line.rstrip()
                        if line.startswith("Long summary for "):
                            # save previous block
                            if currentPose is not None and buffer:
                                txtRows[currentPose] = "\n".join(buffer)
                            # extract pose name correctly
                            posePath = line[len("Long summary for "):].strip()
                            poseFile = posePath.split()[0]
                            currentPose = os.path.basename(poseFile).split('.')[0]
                            buffer = [line]
                        else:
                            if currentPose is not None:
                                buffer.append(line)
                    # save last block
                    if currentPose is not None and buffer:
                        txtRows[currentPose] = "\n".join(buffer)

        return csvRows, txtRows

    def getMode(self):
        if self.useTrueMol.get() and self.molCond.get():
            mode = "true_prot"
        elif self.useTrueMol.get():
            mode = "true"
        elif self.molCond.get():
            mode = "prot"
        else:
            mode = "normal"
        return mode

    def getSelectedColumnAndThreshold(self):
        mode = self.getMode()
        if mode == "normal":
            idx = self.filterColLongTxt.get()
        elif mode == "true":
            idx = self.filterColLongTxtTrueMol.get()
        elif mode == "prot":
            idx = self.filterColLongTxtProt.get()
        else:
            idx = self.filterColLongTxtAll.get()

        colName = self.FILTER_COLUMNS[mode][idx]

        thrMap = {
            'energy_ratio': self.energyRatio.get(),
            'number_clashes': self.numClashes.get(),
            'aromatic_ring_maximum_distance_from_plane': self.ringPlanarity.get(),
            'kabsch_rmsd': self.rmsd.get(),
            'centroid_distance': self.centroidDist.get(),
            'volume_overlap_protein': self.volOverlap.get(),
            'smallest_distance_protein': self.smallDist.get()
        }

        threshold = thrMap.get(colName)
        return colName, threshold

    def valuePasses(self, value):
        _, threshold = self.getSelectedColumnAndThreshold()

        val = float(value)

        op = self.filterOp.get()

        result = False

        if op == 0:
            result = (val == threshold)
        elif op == 1:
            result = (val > threshold)
        elif op == 2:
            result = (val >= threshold)
        elif op == 3:
            result = (val < threshold)
        elif op == 4:
            result = (val <= threshold)
        elif op == 5:
            upper = self.upper.get()
            result = (threshold <= val <= upper)

        return result


    def rowPassesColTxt(self, txtBlock):
        testsToCheck, testsRequired = self.getTestsToPassTxt()
        passed = 0
        for line in txtBlock.splitlines():
            line = line.strip()
            if not line:
                continue
            for test in testsToCheck:
                if line.startswith(test):
                    if '.' in line and 'Fail' not in line:
                        passed += 1
                    break
        return passed >= testsRequired

    def getTestsToPassTxt(self):
        if self.useTrueMol.get() and not self.molCond.get():
            tests = [
                'MOL_PRED loaded',
                'MOL_TRUE loaded',
                'Sanitization',
                'InChI convertible',
                'All atoms connected',
                'No radicals',
                'Molecular formula',
                'Molecular bonds',
                'Double bond stereochemistry',
                'Tetrahedral chirality',
                'Bond lengths',
                'Bond angles',
                'Internal steric clash',
                'Aromatic ring flatness',
                'Non-aromatic ring non-flatness',
                'Double bond flatness',
                'Internal energy',
                'RMSD ≤ 2Å'
            ]
            testNum = self.testsPassedTrueMol.get()
        elif self.useTrueMol.get() and self.molCond.get():
            tests = [
                'MOL_PRED loaded',
                'MOL_COND loaded',
                'Sanitization',
                'InChI convertible',
                'All atoms connected',
                'No radicals',
                'Molecular formula',
                'Molecular bonds',
                'molecular_bonds',
                'Double bond stereochemistry',
                'Tetrahedral chirality',
                'Bond lengths',
                'Bond angles',
                'Internal steric clash',
                'Aromatic ring flatness',
                'Non-aromatic ring non-flatness',
                'Double bond flatness',
                'Internal energy',
                'Protein-ligand maximum distance',
                'Minimum distance to protein',
                'Minimum distance to organic cofactors',
                'Minimum distance to inorganic cofactors',
                'Minimum distance to waters',
                'Volume overlap with protein',
                'Volume overlap with organic cofactors',
                'Volume overlap with inorganic cofactors',
                'Volume overlap with waters',
                'RMSD ≤ 2Å'
            ]
            testNum = self.testsPassedTrueMolProt.get()
        elif self.molCond.get() and not self.useTrueMol.get():
            tests = [
                'MOL_PRED loaded',
                'MOL_COND loaded',
                'Sanitization',
                'InChI convertible',
                'All atoms connected',
                'No radicals',
                'Bond lengths',
                'Bond angles',
                'Internal steric clash',
                'Aromatic ring flatness',
                'Non-aromatic ring non-flatness',
                'Double bond flatness',
                'Internal energy',
                'Protein-ligand maximum distance',
                'Minimum distance to protein',
                'Minimum distance to organic cofactors',
                'Minimum distance to inorganic cofactors',
                'Minimum distance to waters',
                'Volume overlap with protein',
                'Volume overlap with organic cofactors',
                'Volume overlap with inorganic cofactors ',
                'Volume overlap with waters'
            ]
            testNum = self.testsPassedProt.get()
        else:
            tests = [
                'MOL_PRED loaded',
                'Sanitization',
                'InChI convertible',
                'All atoms connected',
                'No radicals',
                'Bond lengths',
                'Bond angles',
                'Internal steric clash',
                'Aromatic ring flatness',
                'Non-aromatic ring non-flatness',
                'Double bond flatness',
                'Internal energy'
            ]
            testNum = self.testsPassedNormal.get()

        return tests, testNum

    def rowPassesColValueTxt(self, txtBlock):
        metrics = {}

        for line in txtBlock.splitlines():
            line = line.strip()
            if not line:
                continue

            parts = line.rsplit(maxsplit=1)
            if len(parts) != 2:
                continue

            key, val = parts
            metrics[key.strip()] = val.strip()

        colName, _ = self.getSelectedColumnAndThreshold()

        val = metrics.get(colName)
        return self.valuePasses(val)

    def getTestInfo(self):
        if self.useTrueMol.get() and not self.molCond.get():
            testNum = self.trueMolTests
            tests = self.testsPassedTrueMol.get()
        elif self.useTrueMol.get() and self.molCond.get():
            testNum = self.trueMolProtTests
            tests = self.testsPassedTrueMolProt.get()
        elif self.molCond.get() and not self.useTrueMol.get():
            testNum = self.protTests
            tests = self.testsPassedProt.get()
        else:
            testNum= self.normalTests
            tests = self.testsPassedNormal.get()
        return tests, testNum

    def molPassesShortTxtTests(self, txtLine):
        txtLine = txtLine.strip()

        match = re.search(r'\((\d+)\s*/\s*(\d+)\)', txtLine)
        if not match:
            print(f"Warning: Could not parse test counts from line:\n{txtLine}")
            return False

        passed = int(match.group(1))
        tests, _ = self.getTestInfo()

        return passed >= tests

    def mergeResultFiles(self, files):
        finalFile = self._getPath(
            'results.csv' if self.outputFormat.get() == 2 else 'results.txt'
        )

        if finalFile.endswith('.csv'):
            self._mergeCsvFiles(files, finalFile)
        else:
            self.mergeTxtFiles(files, finalFile)

    def mergeCsvFiles(self, files, finalFile):
        headerWritten = False

        with open(finalFile, 'w', newline='') as fout:
            writer = None

            for f in files:
                if not os.path.exists(f):
                    continue

                with open(f) as fin:
                    reader = csv.DictReader(fin)

                    if not headerWritten:
                        writer = csv.DictWriter(fout, fieldnames=reader.fieldnames)
                        writer.writeheader()
                        headerWritten = True

                    for row in reader:
                        writer.writerow(row)

    def mergeTxtFiles(self, files, finalFile):
        with open(finalFile, 'w') as fout:
            for f in files:
                if not os.path.exists(f):
                    continue

                with open(f) as fin:
                    fout.write(fin.read())
                    fout.write('\n')

    def getResultFileForMol(self, mol):
        base = os.path.splitext(os.path.basename(mol.getPoseFile()))[0]
        ext = 'csv' if self.outputFormat.get() == 2 else 'txt'
        return self._getExtraPath(f'results_{base}.{ext}')

    def rowPassesShortCsvTests(self, row):
        testsToPass, tests = self.getTestsToPass()

        passed = sum(row.get(test) == 'True' for test in testsToPass)

        return passed >= tests

    def getTestsToPass(self):
        if self.useTrueMol.get() and not self.molCond.get():
            tests = [
                'mol_pred_loaded',
                'mol_true_loaded',
                'sanitization',
                'inchi_convertible',
                'all_atoms_connected',
                'no_radicals',
                'molecular_formula',
                'molecular_bonds',
                'double_bond_stereochemistry',
                'tetrahedral_chirality',
                'bond_lengths',
                'bond_angles',
                'internal_steric_clash',
                'aromatic_ring_flatness',
                'non-aromatic_ring_non-flatness',
                'double_bond_flatness',
                'internal_energy',
                'rmsd_≤_2å'
            ]
            testNum = self.testsPassedTrueMol.get()
        elif self.useTrueMol.get() and self.molCond.get():
            tests = [
                'mol_pred_loaded',
                'mol_true_loaded',
                'mol_cond_loaded',
                'sanitization',
                'inchi_convertible',
                'all_atoms_connected',
                'no_radicals',
                'molecular_formula',
                'molecular_bonds',
                'double_bond_stereochemistry',
                'tetrahedral_chirality',
                'bond_lengths',
                'bond_angles',
                'internal_steric_clash',
                'aromatic_ring_flatness',
                'non-aromatic_ring_non-flatness',
                'double_bond_flatness',
                'internal_energy',
                'protein-ligand_maximum_distance',
                'minimum_distance_to_protein',
                'minimum_distance_to_organic_cofactors',
                'minimum_distance_to_inorganic_cofactors',
                'minimum_distance_to_waters',
                'volume_overlap_with_protein',
                'volume_overlap_with_organic_cofactors',
                'volume_overlap_with_inorganic_cofactors',
                'volume_overlap_with_waters',
                'rmsd_≤_2å'
            ]
            testNum = self.testsPassedTrueMolProt.get()
        elif self.molCond.get() and not self.useTrueMol.get():
            tests = [
                'mol_pred_loaded',
                'mol_cond_loaded',
                'sanitization',
                'inchi_convertible',
                'all_atoms_connected',
                'no_radicals',
                'bond_lengths',
                'bond_angles',
                'internal_steric_clash',
                'aromatic_ring_flatness',
                'non-aromatic_ring_non-flatness',
                'double_bond_flatness',
                'internal_energy',
                'protein-ligand_maximum_distance',
                'minimum_distance_to_protein',
                'minimum_distance_to_organic_cofactors',
                'minimum_distance_to_inorganic_cofactors',
                'minimum_distance_to_waters',
                'volume_overlap_with_protein',
                'volume_overlap_with_organic_cofactors',
                'volume_overlap_with_inorganic_cofactors',
                'volume_overlap_with_waters'
            ]
            testNum = self.testsPassedProt.get()
        else:
            tests = [
                'mol_pred_loaded',
                'sanitization',
                'inchi_convertible',
                'all_atoms_connected',
                'no_radicals',
                'bond_lengths',
                'bond_angles',
                'internal_steric_clash',
                'aromatic_ring_flatness',
                'non-aromatic_ring_non-flatness',
                'double_bond_flatness',
                'internal_energy'
            ]
            testNum = self.testsPassedNormal.get()

        return tests, testNum

    def rowPassesColValue(self, row):
        colName, _ = self.getSelectedColumnAndThreshold()
        csvVal = row.get(colName)

        return self.valuePasses(csvVal)


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

    def getSpecifiedMol(self, string):
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
