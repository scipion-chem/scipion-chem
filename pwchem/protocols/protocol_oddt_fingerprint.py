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
    Perform 2D molecular descriptors comparison between a set of molecules in smi, sdf, pdb, mol or mol2 format and a target molecule.
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

        form.addParam('referenceMolecule', BooleanParam, default=True,
                      label='is the reference molecule in the input set?')

        form.addParam('inputReferenceMolecule', params.StringParam, condition='referenceMolecule',
                      label="Reference molecule: ",
                      help='Model molecule to compare others')

        form.addParam('inputReferenceSmile', params.StringParam, condition='not referenceMolecule',
                      label="Reference SMILE: ",
                      help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('fpChoice', params.EnumParam, default=0, label='Fingerprint type: ',
                       choices=self._fpChoices,
                       help="Chosen fingerprint type to perform the filtering")

        group.addParam('coefChoice', params.EnumParam, default=0, label='Similarity coefficient: ',
                        choices=self._coefChoices,
                        help="Chosen fingerprint type to perform the filtering")

        group.addParam('cut', params.FloatParam, default=0, label='Filter cut-off: ',
                       help="Filter cut-off")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        print(self.inputReferenceSmile.get())
        self.describeFilter(mols)  # , receptorFile)

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods = "The protocol analyses the similarity of two molecules using their fingerprints since Molecular fingerprints are essential cheminformatics tools for virtual screening and mapping chemical"
        methods1 = """The protocol offers the user the option of filtering by two types of fingerprints:\n - MACCS (Molecular ACCess System): these fingerprints consist of 166 bits representing predefined structural fragments. At each position, the presence or absence of a fragment is reported. \n - Morgan (circular fingerprints): this type of fingerprint is based on the Morgan algorithm. Each bit corresponds to a circular environment"""
        methods2 = """The similarity between both fingerprints is measured by different similarity coefficients:\n- Tanimoto coefficient: it is the most popular in both chemical informatics and computational medicinal chemistry because of its ease of implementation and quickness.\n - Dice coefficient: also known by other names such as Sørensen&#39;s index, Dice coefficient is a “statistic introduced to compare the similarity of two sample”. """
        methods3 = """because of each coefficient&#39;s calculation, Dice&#39;s coefficient tends to have higher values than those from Tanimoto one."""
        return methods, methods1, methods2, methods3

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

    # Funcion importante
    def describeFilter(self, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion)  # , receptorFile)
        # Aqui mi script: hay que pasarle el nombre del script y el path del documento de parametros que se crea con otra funcion
        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())
        #if paramsPath:
            #os.remove(paramsPath)

    def writeParamsFile(self, paramsFile, molsScipion):  # , recFile):
        # Aqui creo que creamos una lista de los archivos que vamos a usar
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

        else:
            for mol in molsScipion:
                molFiles.append(os.path.abspath(mol.getFileName()))

            objective_smile = self.inputReferenceSmile.get()

            f.write('referenceFile: {}\n'.format(objective_smile))

        # escribir el tipo de descriptor, en nuestro caso el tipo de regla
        f.write('fingerprint: {}\n'.format(self.getFingerprint()))
        f.write('coefficient: {}\n'.format(self.getCoefficient()))
        f.write('cut-off: {}\n'.format(self.cut.get()))
        f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile


    def getFingerprint(self):
        function = self.getEnumText('fpChoice')
        return function

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


        #print(molecules)
        return molecules

