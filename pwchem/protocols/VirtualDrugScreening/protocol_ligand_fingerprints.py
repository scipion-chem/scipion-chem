# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *
# # # *
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


from pyworkflow.protocol import params
import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re
from pwchem import Plugin
from pyworkflow.utils.path import copyFile
from pwchem.utils import fillEmptyAttributes
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam
import csv
import glob

scriptName = 'fingerprint_filtering.py'


class ProtocolFingerprintFiltering(EMProtocol):
    """
    Perform 2D molecular descriptors comparison between a set of molecules in smi, sdf, pdb, mol or mol2 format and a
    target molecule.
    Fingerprints are calculated for all compounds (MACCS or Morgan) and
    compared through a similarity coeffciente which ranges from 0 (no similarity) to 1 (top similarity).
    It is possible to choose between Tanimoto and Dice coeffcient as similarity measures

    """
    _label = 'Fingerprint filtering'
    _fpChoices = ['Morgan', 'MACCS']
    _coefChoices = ['Tanimoto', 'Dice']


    ##### -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules: ",
                      help='Select the molecules to be filtered')

        form.addParam('inputRefSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Reference Small Molecules: ",
                      help='Select the molecules where the reference molecule is stored')

        form.addParam('inputReferenceMolecule', params.StringParam,
                      label="Reference molecule: ",
                      help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('fpChoice', params.EnumParam, default=0, label='Fingerprint type: ',
                       choices=self._fpChoices,
                       help="Chosen fingerprint type to perform the filtering")

        group.addParam('coefChoice', params.EnumParam, default=0, label='Similarity coefficient: ',
                        choices=self._coefChoices,
                        help="Chosen fingerprint type to perform the filtering")

        group.addParam('cut', params.FloatParam, default=0, label='Filter cut-off: ',
                       help="Filter cut-off, distances below the threshold will be filtered")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, mols)
        Plugin.runScript(self, scriptName, paramsPath, env='rdkit', cwd=self._getPath())

    # --------------- INFO functions -------------------------

    def _citations(self):
        return []

    def _methods(self):
        return []

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if os.path.exists(self._getPath("all_distances.csv")):
            summary.append("The distance results for each metric and molecule can be found at {}".
                           format(self._getPath("all_distances.csv")))
        else:
            summary.append("The protocol has not finished.")
        return summary

    # --------------------------- UTILS functions -----------------------------------

    def createOutputStep(self):
        filtered_molecules_dict = self.parseResults(self._getPath("results.tsv"))
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        filtered_molecules = list(filtered_molecules_dict.keys())

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            file = os.path.abspath(mol.getFileName())
            if file in filtered_molecules:
                mol._coefficient = pwobj.String(filtered_molecules_dict[file])
                newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    def writeParamsFile(self, paramsFile, molsScipion):
        molFiles = []
        file_mol_dict = {}
        f = open(paramsFile, 'w')
        f.write('outputPath: results.tsv\n')
        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))

        f.write('referenceFile: {}\n'.format(self.getRefFile()))
        f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        f.write('fingerprint: {}\n'.format(self.getEnumText('fpChoice')))
        f.write('coefficient: {}\n'.format(self.getEnumText('coefChoice')))
        f.write('cut-off: {}\n'.format(self.cut.get()))

        return paramsFile


    def getRefFile(self):
        if not hasattr(self, 'refFile'):
            for mol in self.inputRefSmallMolecules.get():
                if mol.__str__() == self.inputReferenceMolecule.get():
                    self.refFile = os.path.abspath(mol.getFileName())
        return self.refFile

    def getCoefficient(self):
        function = self.getEnumText('coefChoice')
        return function

    def parseResults(self, outputFile):
        molecules = {}
        with open(outputFile) as read_tsv:
            for row in read_tsv:
                if row[0] == "#":
                    pass
                else:
                    row1 = row.split(",")
                    molecules[row1[0]] = row1[1]
        return molecules


