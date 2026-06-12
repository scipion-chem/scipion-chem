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

# General imports
import os

# Scipion em imports
import pyworkflow.object as pwobj
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem.constants import RDKIT_DIC

scriptName = 'RDKit_PAINS_script.py'

class ProtocolPainsRdkitFiltering(EMProtocol):
    """
    AI Generated:

This protocol performs PAINS (Pan-Assay INterference Compounds) filtering on a set of small molecules.

PAINS are molecular substructures known to frequently generate false positives in high-throughput screening assays
due to non-specific interactions with biological targets. This protocol identifies and separates molecules that
contain PAINS patterns from those that do not.

Core Concepts
-------------
Small Molecules:
    Input dataset consisting of chemical compounds in formats such as SMI, SDF, PDB, MOL, or MOL2.

PAINS Substructures:
    Predefined or user-provided SMARTS patterns corresponding to problematic chemical motifs associated with
    assay interference and non-specific binding.

Filtering Process:
    Each molecule is analyzed to detect the presence of PAINS patterns. Based on the results, molecules are
    classified into clean compounds or PAINS-positive compounds.

PAINS Sources
-------------
RDKit Catalogue:
    Default PAINS substructure definitions provided by the RDKit library.

Custom PAINS File:
    User-defined file containing PAINS patterns. Each line must include:
    - A SMARTS pattern describing the substructure
    - A name/label for the PAINS pattern

Workflow
--------
1. Input a set of small molecules.
2. Convert molecular formats if required for compatibility.
3. Apply PAINS detection using RDKit or a custom PAINS catalogue.
4. Identify molecules containing PAINS substructures.
5. Split molecules into:
   - Clean molecules (no PAINS detected)
   - PAINS-positive molecules (contain at least one PAINS pattern)
6. Annotate PAINS-positive molecules with detected pattern information.
7. Generate separate output datasets.

Output
------
- outputSmallMolecules:
    Set of molecules that passed the PAINS filter (clean compounds).

- outputSmallMoleculesPains:
    Set of molecules flagged as PAINS-positive, annotated with detected PAINS information.

Use Cases
---------
- Pre-filtering compound libraries before virtual screening
- Removal of assay-interfering compounds in drug discovery pipelines
- Quality control of chemical libraries
- Identification of promiscuous binders in screening datasets

Notes
-----
- PAINS filtering is a heuristic method and may produce false positives or negatives.
- It should be used as a prioritization tool rather than an absolute exclusion criterion.
- Custom PAINS definitions allow adaptation to specific screening campaigns or chemical spaces.
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

        form.addParam('referenceMolecule', params.BooleanParam, default=False,
                      label='Custom PAINS file: ',
                      help='Whether to use RDKit default PAINS substructures or provide a file with custom PAINS'
                           ' substructures.')

        form.addParam('painsFile', params.PathParam, condition='referenceMolecule',
                      label="Reference file: ",
                      help='Custom PAINS file to use for filtering. Each line must contain a first column with a '
                           'SMARTS string and a second column with a short description.')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        mols = self.inputSmallMolecules.get()
        results = self.describeFilter(mols)


    def createOutputStep(self):
        list_mols = []
        filtered_molecules_names = self.parseResults(self._getPath("results.tsv"))
        pains_dict = self.parsePains(self._getPath("with_pains.txt"))
        pains_molNames = list(pains_dict.keys())

        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        newMols2 = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(),
                                                 suffix='_pains', copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            file = os.path.abspath(mol.getFileName())
            molName = getBaseName(file)
            if molName in filtered_molecules_names:
                newMols.append(mol)
            elif molName in pains_molNames:
                mol._pain = pwobj.String(pains_dict[molName])
                newMols2.append(mol)

        newMols.updateMolClass()
        newMols2.updateMolClass()

        if newMols:
            self._defineOutputs(outputSmallMolecules=newMols)
        if newMols2:
            self._defineOutputs(outputSmallMoleculesPains=newMols2)

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
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())
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
            for line in tsv_file:
                if not line[0] == "#":
                    molName = getBaseName(line)
                    molecules.append(molName)

        return molecules

    def parsePains(self, file):
        pains_m = {}
        with open(file) as tsv:
            for row in tsv:
                if not row[0] == "#":
                    molFile, pain = row.split(",")
                    molName = getBaseName(molFile)
                    pains_m[molName] = pain

        return pains_m

