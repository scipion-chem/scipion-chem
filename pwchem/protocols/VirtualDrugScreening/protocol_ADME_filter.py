#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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

# General imports
import os

# Scipion em imports
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import RDKIT_DIC
from pwchem.utils import getBaseName, runOpenBabel

scriptName = 'ADME_script.py'

LIPINSKI, RO3 = 'Lipinski Rule of 5', 'Rule of 3'
_ruleChoices = [LIPINSKI, RO3]
RDIC = {LIPINSKI: 'ro5', RO3: 'ro3'}

class ProtocolADMEFiltering(EMProtocol):
    """
    Performs the analysis of ADME (Absortion, Distribution, Metabolism, Excretion) properties.
    This filtering is particularly useful for testing compounds whose pharmacological
    properties/activity have already been proven.

    It is possible to analyze molecules in  mol, mol2, pdb, sdf and smi format.
    """
    _label = 'ADME ligand filtering'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input  Small Molecules: ",
                      help='Select the molecules to be filtered')

        group = form.addGroup('Descriptors')
        group.addParam('ruleChoice', params.EnumParam, default=0, label='Filtering Rule: ',
                       choices=_ruleChoices,
                       help="Chosen rule to perform the ADME (Absortion, Distribution, Metabolism, Excretion) filtering.\nro5: Lipinski's rule of five.\n"
                            "ro3: Rule of three\n"
                            "https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        self.describeFilter(mols)

    def createOutputStep(self):

        filtered_molecules_names = self.parseResults(self._getPath("results.tsv"))
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            molFile = os.path.abspath(mol.getFileName())
            molName = getBaseName(molFile)
            if molName in filtered_molecules_names:
                newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------- INFO functions -------------------------
    def _citations(self):
        return ['Wójcikowski2015']

    # --------------------------- UTILS functions -----------------------------------
    def describeFilter(self, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion)  # , receptorFile)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    def writeParamsFile(self, paramsFile, molsScipion):
        molFiles = []

        for mol in molsScipion:
            molFile = os.path.abspath(mol.getFileName())
            if molFile.endswith('.pdbqt'):
                inpPDBFile = os.path.abspath(self._getTmpPath(os.path.basename(molFile).split('.')[0] + '.pdb'))
                args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(molFile), inpPDBFile)
                runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
                molFile = inpPDBFile
            molFiles.append(molFile)

        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('rule: {}\n'.format(self.getRule()))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile

    def getRule(self):
        return RDIC[self.getEnumText('ruleChoice')]

    def parseResults(self, outputFile):
        molecules = []
        with open(outputFile) as tsv_file:
            for line in tsv_file:
                if line[0] != "#":
                    molName = getBaseName(line)
                    molecules.append(molName)

        return molecules


