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
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re
from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes
import csv

scriptName = 'shape_filtering.py'


class ProtocolShapeFiltering(EMProtocol):
    """
    Performs the RMSD calculation between a molecule in smi format and a query.

    """
    _label = 'Shape filtering'
    _RMSDChoices = ['BestRMSD', 'RMSD']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        form.addParam('inputSmallMolecule', params.PointerParam,
                      pointerClass='SmallMolecule', allowsNull=False,
                      label="Input reference: ",
                      help='Select the reference molecule')

        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules: ",
                      help='Select the molecules to be filtered')

        group = form.addGroup('Descriptors')
        group.addParam('rmsdChoice', params.EnumParam, default=0, label='Filter type: ',
                       choices=self._RMSDChoices,
                       help="Chosen shape filter")

        group.addParam('cut', params.FloatParam, default=0, label='Filter cut-off: ',
                       help="Filter cut-off")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        reference_mol = self.inputSmallMolecule.get()
        results = self.describeFilter(reference_mol, mols)

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

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@misc{landrum _2021, title={Rdkit.chem.rdmolalign module¶}, url={https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html}, journal={rdkit.Chem.rdMolAlign module - The RDKit 2022.03.1 documentation}, author={Landrum , Greg}, year={2021}}"]

    def _methods(self):
        methods = "This filter is based on the analysis of the shape of the query molecule/ligand used as a reference in the search for similar structures. A commonly used measure for comparing the structure of two molecules in virtual screening programs is the calculation of the root mean square deviation (RMSD) between the atoms of two molecules. The RMSD is a distance which describes the structural difference between two topologies. The lower the RMSD between two structures, the greater the similarity between them."
        methods2 = "It is possible to estimate the optimally (minimum) RMSD between two molecules by “rdMolAlign” module function “AlignMol()” and it is also possible to calculate the best RMSD (“GetBestRMS()”) by aligning all permutations of matching atoms orders in both molecules"
        methods3 = "The difference between both functions is that “AlignMol()” aligns molecules without changing atom order"
        return methods, methods2, methods3

    # --------------------------- UTILS functions -----------------------------------

    def describeFilter(self, reference_mol, molsScipion):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, reference_mol, molsScipion)
        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())

    def writeParamsFile(self, paramsFile, reference_mol, molsScipion):
        molFiles = []

        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))

        reference_mol = os.path.abspath(reference_mol.getFileName())

        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('filterChoice: {}\n'.format(self.getFilterChoice()))
            f.write('cut-off: {}\n'.format(self.cut.get()))
            f.write('referenceFile: {}\n'.format(reference_mol))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile

    def getFilterChoice(self):
        function = self.getEnumText('rmsdChoice')
        return function

    def parseResults(self, outputFile):
        molecules = []
        with open(outputFile) as tsv_file:
            read_tsv = csv.reader(tsv_file, delimiter="\n")
            for row in read_tsv:
                if row[0] == "#":
                    next
                else:
                    molecules.append(row)

        molecules = sum(molecules, [])
        print(molecules)
        return molecules

