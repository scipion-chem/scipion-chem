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

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.utils import *
import os
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
import csv


scriptName = 'filter_oddt.py'


class ProtocolODDTFiltering(EMProtocol):
    """
    Performs the analysis of ADME properties.
    This filtering is particularly useful for testing compounds whose pharmacological
    properties/activity have already been proven.

    It is possible to analyze molecules in  mol, mol2, pdb, sdf and smi format.
    """
    _label = 'LB filtering'
    _ruleChoices = ['ro5', 'l5', 'ro3']

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
                       choices=self._ruleChoices,
                       help="Chosen rule to perform the filtering")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')


    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        self.describeFilter(mols)

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods = "The objective is the qualitative evaluation of the screening compounds pharmacokinetic properties. These properties give a notion of whether the analyzed compounds can fulfil their function as a drug when they are ingested by humans."
        methods1 = "Lipinski rules (Ro5, l5): - Contain no more than 5 hydrogen bond donors (the total number of nitrogen–hydrogen and oxygen–hydrogen bonds).\n - Contain no more than 10 hydrogen bond acceptors (all nitrogen or oxygen atoms).\n  - Have a molecular mass of 160 or less than 500 Da, as large molecules (high molecular weight) might have more difficulties in passing phospholipid membranes.\n - Not exceed 5 in its octanol-water partition coefficient (log P). This coefficient measures the distribution of a compound, usually between a hydrophobic (e.g. 1-octanol) and a hydrophilic (e.g. water) phase."
        methods2 = "Ro3: - Contain no more than 3 hydrogen bond donors (the total number of nitrogen–hydrogen and oxygen–hydrogen bonds).\n - Contain no more than 3 hydrogen bond acceptors (all nitrogen or oxygen atoms).\n  - Have a molecular mass of 160 or less than 300 Da, as large molecules (high molecular weight) might have more difficulties in passing phospholipid membranes.\n - Not exceed 3 in its octanol-water partition coefficient (log P). This coefficient measures the distribution of a compound, usually between a hydrophobic (e.g. 1-octanol) and a hydrophilic (e.g. water) phase."

        return methods, methods1, methods2

    # --------------------------- UTILS functions -----------------------------------

    def createOutputStep(self):

        filtered_molecules = self.parseResults(self._getPath("results.tsv"))
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            file = os.path.abspath(mol.getFileName())
            if file in filtered_molecules:
                newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    def describeFilter(self, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion)  # , receptorFile)
        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())
        #if paramsPath:
            #os.remove(paramsPath)

    def writeParamsFile(self, paramsFile, molsScipion):
        molFiles = []

        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))  # que es getPoseFile ¿???¿????

        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('rule: {}\n'.format(self.getRule()))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile

    def getRule(self):
        function = self.getEnumText('ruleChoice')
        return function

    def parseResults(self, outputFile):
        molecules = []
        with open(outputFile) as tsv_file:
            read_tsv = csv.reader(tsv_file, delimiter="\n")
            for row in read_tsv:
                if row[0] == "#":
                    pass
                else:
                    molecules.append(row)

        molecules = sum(molecules, [])
        print(molecules)
        return molecules

