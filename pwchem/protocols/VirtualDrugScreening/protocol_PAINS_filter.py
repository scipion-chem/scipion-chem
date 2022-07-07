# Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
import glob
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re
from pwchem import Plugin
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.utils.path import copyFile
from pwchem.utils import fillEmptyAttributes
import csv
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam
import pyworkflow.object as pwobj

scriptName = 'RDKit_PAINS_script.py'


class ProtocolPainsRdkitFiltering(EMProtocol):
    """
    Perform the PAINS search in molecules from smi, sdf, pdb, mol or mol2 file format.
    PAINS file can be uploaded and must have the following format: on each line, the first position is reserved for the PAINS molecule
    in SMART format, and the second position for the PAINS name, both elements must be
    space separated; or it can be used the predefined PAINS catalogue offered by the RDKit library.

    """
    _label = 'PAINS RDKit filtering'
    _possibleOutputs = {"SmallMoleculesCleaned": SetOfSmallMolecules,"SmallMoleculesPains": SetOfSmallMolecules}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input  Small Molecules: ",
                      help='Select the molecules to be filtered')

        form.addParam('referenceMolecule', BooleanParam, default=True,
                      label='have you a file os PAINS to perform the filtering?')

        form.addParam('painsFile', params.PathParam, condition='referenceMolecule',
                      label="Reference file: ",
                      help='Import the path PAINS file')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        results = self.describeFilter(mols)


    def createOutputStep(self):
        list_mols = []
        filtered_molecules = self.parseResults(self._getPath("results.tsv"))
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(),
                                                 suffix='consensusSmallMolecules{}', copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            file = os.path.abspath(mol.getFileName())
            newMols.append(mol)

        newMols.updateMolClass()
        list_mols.append(newMols)
        #self._defineOutputs(outputSmallMolecules=newMols)

        if not self.painsFile.get():
            pains_dict = self.parsePains(self._getPath("with_pains.txt"))
            pains = list(pains_dict.keys())
            newMols2 = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

            mols2 = self.inputSmallMolecules.get()
            for mol2 in mols2:
                file2 = os.path.abspath(mol2.getFileName())
                if file2 in pains:
                    mol2._pain = pwobj.String(pains_dict[file2])
                    newMols2.append(mol2)
            newMols2.updateMolClass()

            list_mols.append(newMols2)
            #self._defineOutputs = ({"without_pains": newMols})

            if newMols:
                self._defineOutputs(SmallMoleculesCleaned=newMols)
            if newMols2:
                self._defineOutputs(SmallMoleculesPains=newMols2)

    # --------------- INFO functions -------------------------

    def _citations(self):
        return [
            "@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods = "PAINS (Pan-Assay INterference compoundS) are compounds that frequently show erroneous results (false positives) in high-throughput screening (HTS). This is a consequence of the functional groups they contain which cause non-specific binding to a wide range of targets. The protocol goal is to detect PAINS structures in a set of molecules. In the protocol, it is possible to choose between the predefined PAINS catalogue offered by the RDKit library or to use a catalogue provided by the user. This user supplied file must be structured as follows: on each line, the first position is reserved for the PAINS molecule in SMART format, and the second position for the PAINS name, both elements must be space separated."
        return methods
    # --------------------------- UTILS functions -----------------------------------

    # Funcion importante
    def describeFilter(self, molsScipion):  # , receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion)
        Plugin.runScript(self, scriptName, paramsPath, env='rdkit', cwd=self._getPath())
        #if paramsPath:
            #os.remove(paramsPath)

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
            # Esta linea vale --> es la que lleva al archivo output

            f.write('outputPath: results.tsv\n')
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

            if self.painsFile.get():
                f.write('painsFile: {}\n'.format(self.painsFile.get()))
            else:
                pass

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

    def parsePains(self, file):
        pains_m = {}
        with open(file) as tsv:
            for row in tsv:
                if row[0] == "#":
                    pass
                else:
                    row2 = row.split(",")

                    pains_m[row2[0]] = row2[1]

        return pains_m

