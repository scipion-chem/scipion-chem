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

scriptName = 'shape_script2.py'


class ProtocolShapeDistancesFiltering(EMProtocol):
    """
    Performs the calculation of Tanimoto and Protrude distances between a molecule in smi format and a query.

    """

    _label = 'Shape Distance filtering'
    _distanceType = ['Tanimoto Distance', 'Protude Distance']

    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules: ",
                      help='Select the molecules to be filtered')

        form.addParam('referenceMolecule', BooleanParam, default=True,
                      label='is the reference molecule in the input set?')

        form.addParam('inputReferenceMolecule', params.StringParam, condition='referenceMolecule',
                      label="Reference molecule: ",
                      help='Model molecule to compare others')

        form.addParam('inputReferenceSmile', params.StringParam, condition='not referenceMolecule',
                      label="Reference SMILE: ",
                      help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('distanceType', params.EnumParam, default=0, label='Distance type: ',
                       choices=self._distanceType,
                       help="Chosen distance type to perform the filtering")

        form.addParam('hydrogen', BooleanParam, default=True,
                      label='Do you want to ignore the hydrogens in the filtering?')

        group.addParam('cut', params.FloatParam, default=0, label='Filter cut-off: ',
                       help="Filter cut-off")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        print(self.inputReferenceSmile.get())
        # receptorFile = self.inputSmallMoleculesSets[0].get().getProteinFile() #getProteinFile() es una funcion de objects.py
        self.describeFilter(mols)

    # --------------------------- UTILS functions -----------------------------------

    def createOutputStep(self):
        archivo_output = self._getPath("results.tsv")

        if archivo_output:
            filtered_molecules_dict = self.parseResults(self._getPath("results.tsv"))
            newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
            filtered_molecules = list(filtered_molecules_dict.keys())

            mols = self.inputSmallMolecules.get()
            for mol in mols:
                file = os.path.abspath(mol.getFileName())
                if file in filtered_molecules:
                    mol._shapeDistance = pwobj.String(filtered_molecules_dict[file])
                    newMols.append(mol)

            newMols.updateMolClass()
            self._defineOutputs(outputSmallMolecules=newMols)

        else:
            caution = open("warning.txt", "w")
            caution.write("The ChEMBL query does not exist, try other parameters")
            caution.close()

    def describeFilter(self, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion)  # , receptorFile)

        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())
        # if paramsPath:
        # os.remove(paramsPath)

    def writeParamsFile(self, paramsFile, molsScipion):
        molFiles = []
        file_mol_dict = {}
        f = open(paramsFile, 'w')
        f.write('outputPath: results.tsv\n')
        if self.referenceMolecule.get():
            for mol in molsScipion:
                molFiles.append(os.path.abspath(mol.getFileName()))  # que es getPoseFile ¿???¿????
                name = mol.__str__()
                path = os.path.abspath(mol.getFileName())
                file_mol_dict[name] = path
            string_objetive = self.inputReferenceMolecule.get()
            objective_file = file_mol_dict[string_objetive]

            f.write('referenceFile: {}\n'.format(objective_file))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        else:
            for mol in molsScipion:
                molFiles.append(os.path.abspath(mol.getFileName()))

            objective_smile = self.inputReferenceSmile.get()
            f.write('referenceFile: {}\n'.format(objective_smile))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        f.write('distanceType:{}\n'.format(self.getDistance()))
        f.write('ignoreHydrogen:{}\n'.format(self.hydrogen.get()))
        f.write('cut-off: {}\n'.format(self.cut.get()))

        return paramsFile
    # --------------- INFO functions -------------------------

    def _citations(self):
        return [
            "@misc{landrum _2021, title={Rdkit.chem.rdmolalign module¶}, url={https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html}, journal={rdkit.Chem.rdMolAlign module - The RDKit 2022.03.1 documentation}, author={Landrum , Greg}, year={2021}}"]

    def _methods(self):
        methods = "This filter is based on the analysis of the shape of the query molecule/ligand used as a reference in the search for similar structures. A commonly used measure for comparing the structure of two molecules in virtual screening programs is the calculation of the root mean square deviation (RMSD) between the atoms of two molecules. The RMSD is a distance which describes the structural difference between two topologies. The lower the RMSD between two structures, the greater the similarity between them."
        methods2 = "RDKit current tools for calculating distances between molecules of unequal size are Protrude Distance and Tanimoto distance. Protrude Distance focusses on the volume mismatch and is defined as the percentage of the larger molecule which protrudes/exceeds from the smaller molecule. Tanimoto Distance measures similarity between finite sample sets and is defined as the size of the intersection divided by the size of the union of the sample sets."
        methods3 = "RDKit does not have a tool that performs an alignment prior to the calculation of similarity, so the results provided by these tools must be carefully analyzed by users."
        return methods, methods2, methods3

    # --------------------------- UTILS functions -----------------------------------

    def getDistance(self):
        function = self.getEnumText('distanceType')
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

        print(molecules)
        return molecules

