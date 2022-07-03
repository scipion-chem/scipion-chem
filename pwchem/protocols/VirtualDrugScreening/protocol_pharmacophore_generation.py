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

scriptName = 'PBVS_script.py'

class ProtocolPharmacophoreFiltering(EMProtocol):
    """
    Perform the construction of a consensus pharmacophore from a set of PDB ligands and
    the PBVS (using the pharmacophore) against a set of molecules in smi, pdb, mol, mol2 or sdf file format.

    """
    #Nombre del protocolo (aparece en grande arriba)
    _label = 'Pharmacophore filtering'

    ##### -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputLigands', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input reference ligands: ",
                      help='Select the ligands PDB files')

        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules: ",
                      help='Select the molecules to be filtered')

        group = form.addGroup('Descriptors')
        #group.addParam('smileFile', params.PathParam, label='CSV FILE: ',
                       #help="CSV file which contains the info of the ligands")
        group.addParam('numberPharmacophoreFeatures', params.FloatParam, default = 0, label='Number of pharmacophore features: ',
                       help="Which is the number of pharmacophoric features that needs a molecule to pass the filter?")

        # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        ligands = self.inputLigands.get()
        mols = self.inputSmallMolecules.get()
        self.describeFilter(ligands, mols)  # , receptorFile)

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods ="The first part of the process involves the generation of the consensus pharmacophore from the PDB ligands. The second protocol step is the analysis of starting molecules pharmacophore characteristics using the consensus pharmacophore."
        return methods

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

    def describeFilter(self, ligands, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, ligands,molsScipion)  # , receptorFile)

        args = ' {} {}'.format(paramsPath, os.path.abspath(self._getPath()))
        Plugin.runScript(self, scriptName, args, env='rdkit', cwd=self._getPath())

    def writeParamsFile(self, paramsFile, ligands, molsScipion):
        molFiles = []
        ligandsFiles = []
        smiles_ligand = []
        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))  # que es getPoseFile ¿???¿????

        smiles = []
        for ligand in ligands:
            fnLigand = os.path.abspath(ligand.getPoseFile())
            ligandName = os.path.basename(fnLigand).split('.')[0]

            args = ' -i "{}" --outputDir {} -of smi -o {} --overWrite'.\
                format(fnLigand, os.path.abspath(self._getExtraPath()), ligandName)
            Plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=self._getExtraPath())

            with open(self._getExtraPath(ligandName + '.smi')) as f:
                smile = f.read().split()[0].strip()
            if not smile in smiles:
                ligandsFiles.append(fnLigand)
                smiles.append(smile)
                id = ((ligand.getPoseFile()).split("/")[-1]).split('.')[0]
                smiles_ligand.append([id, smile])

        print(smiles_ligand)
        with open(paramsFile, 'w') as f:
            # Esta linea vale --> es la que lleva al archivo output

            f.write('outputPath: results.tsv\n')
            f.write('ligandSmiles:' + str(smiles_ligand) + "\n")
            f.write('ligandFiles: {}\n'.format(' '.join(ligandsFiles)))
            f.write('moleculesFiles: {}\n'.format(' '.join(molFiles)))
            f.write('featuresNumber: {}\n'.format(self.numberPharmacophoreFeatures.get()))
        return paramsFile


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
