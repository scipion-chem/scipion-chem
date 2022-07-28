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

import glob
from pyworkflow.utils.path import copyFile
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re
from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes
import csv
import pyworkflow.object as pwobj

scriptName = 'chembl_accession.py'


class ProtocolChemblAccession(EMProtocol):
    """
    Perform the search of compounds that mimic the biological effects of a drug against a target protein

    """
    _label = 'Chembl accession'
    _idOrigins = ['Chembl', 'Uniprot']
    _moleculetypeChoices = ['Any', 'SINGLE PROTEIN', 'CHIMERIC PROTEIN', 'PROTEIN FAMILY', 'PROTEIN-PROTEIN INTERACTION']
    _assaytypeChoices = ['Any', '(B) Binding', '(F) Functional', '(A) ADMET', '(T) Toxicity', '(P) Physicochemical', '(U) Unclassified']
    _relationChoices = ['=']
    _typeChoices = ['IC50', 'Ki', 'Kd', 'EC50']
    # -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        group = form.addGroup('Input')
        group.addParam('idOrigin', params.EnumParam, default=0, label='ID database: ',
                       choices=self._idOrigins,
                       help="Databse where the ID comes from")
        group.addParam('inputID', params.StringParam,
                      allowsNull=False, label="Uniprot ID: ",
                      help='Select the Uniprot ID of the target molecule')
        group.addParam('targetOrganism', params.StringParam, label="Target organism: ",
                       condition='idOrigin==1', default='',
                       help='Select the organism of the target protein. Empty if any')
        group.addParam('moleculetypeChoice', params.EnumParam, default=0, label='Target type: ',
                       choices=self._moleculetypeChoices, condition='idOrigin==1',
                       help="Target type of the query molecule")

        group = form.addGroup('Descriptors')
        group.addParam('assaytypeChoice', params.EnumParam, default=0, label='Assay type: ',
                       choices=self._assaytypeChoices,
                       help="Classification of the assay searched a. Binding (B): it is based on the binding of ligand molecules to receptors, antibodies, or other macromolecules. This type of assays obtain data for measuring binding of compound to a molecular target, e.g. Ki, IC50, Kd.\n b. Functional (F): a set of systematic in vivo experiments designed to measure the effect or biological role of a compound in a cellular pathway or biological process, e.g., cell death in a cell line.\n c. ADMET (A): assays related to the pharmacokinetic characteristics of a compound e.g., t1/2, oral bioavailability.\n d. Toxicity (T): Data measuring toxicity of a compound, e.g., cytotoxicity.\n e. Physicochemical (P): this type of assay is performed in the absence of biological material and measures physicochemical properties of the compounds e.g., chemical stability, solubility\n f. Unclassified (U): those assays which cannot be included in the rest of type are classified into one of the above categories e.g., ratio of binding vs efficacy.")

        group.addParam('relationChoice', params.EnumParam, default=0, label='Relation type: ',
                       choices=self._relationChoices,
                       help="Relation Choice")

        group.addParam('bioactivityChoice', params.EnumParam, default=0, label='Bioactivity type: ',
                       choices=self._typeChoices,
                       help="Bioactivity measure. a. IC50: (Half maximal inhibitory concentration) indicates the required concentration of a compound to inhibit a biological process by 50%. The lower the IC50 value, the higher the compound bioactivity.\n b. EC50: (Half maximal effective concentration) is the drug concentration required to produce half (50%) of the maximum effect of that compound. The lower the EC50 value, the higher the compound bioactivity. EC50 is equivalent to IC50, but this last is used in the case of working with inhibitors.\n c. Kd: dissociation constant, this term is generic and describes the binding affinity between a molecule and an enzyme or a receptor. \n d. Ki: inhibition constant, the term is equivalent to the dissociation constant (Kd), but it is used when the molecule for which bioactivity is measured is an inhibitor.")

        group.addParam('numberOfStructures', params.IntParam, default=20, label='Number of final molecules : ',
                       help="Number of small molecules chosen after the filtering")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('describeFilter')
        self._insertFunctionStep('importStep')

    def filterStep(self):
        results = self.describeFilter()

    def describeFilter(self):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath)
        Plugin.runRDKitScript(self, scriptName, paramsPath, cwd=self._getPath())
        #if paramsPath:
            #os.remove(paramsPath)

    def importStep(self):
        file_path = self.createOutput()
        fh = open(file_path)
        dic_atrib = {}
        for line in fh.readlines():
            tokens = line.split(',')
            if len(tokens) == 3:
                fhSmile = open(self._getExtraPath(tokens[1].strip() + ".smi"), 'w')
                fhSmile.write(tokens[0].strip() + "\n")
                file_smi = self._getExtraPath(tokens[1].strip() + ".smi")
                atribute = tokens[2]
                dic_atrib[file_smi] = atribute

                fhSmile.close()
        fh.close()

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')

        for fnSmall in glob.glob(self._getExtraPath("*.smi")):
            smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
            smallMolecule._bioactivity = pwobj.String(dic_atrib[fnSmall])
            #añadiratributo
            if len(os.listdir(self._getExtraPath())) <= 100:  # costly
                if not fnSmall.endswith('.mae') or not fnSmall.endswith('.maegz'):
                    fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                    fnOut = self._getExtraPath("%s.png" % fnRoot)
                    args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (fnSmall, fnOut)
                    try:
                        Plugin.runRDKit(self, "python3", args)
                        smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                    except:
                        smallMolecule._PDBLigandImage = pwobj.String("Not available")

            outputSmallMolecules.append(smallMolecule)
        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods ="""Protocol filters the compounds according to the following parameters:\n 
        - Uniprot ID \n
        - Assay type \n
        - Bioactivity measure \n
        - Target type \n
        - Organism of compound origin"""
        return methods

    # --------------------------- UTILS functions -----------------------------------

    def writeParamsFile(self, paramsFile):
        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('idOrigin: {}\n'.format(self.getEnumText('idOrigin')))
            f.write('inputID: {}\n'.format(self.inputID.get()))
            f.write('targetOrganism: {}\n'.format(self.targetOrganism.get()))
            f.write('targetType: {}\n'.format(self.getTargetTypeChoice()))

            f.write('assayType: {}\n'.format(self.getAssayTypeChoice()))
            f.write('relation: {}\n'.format(self.getRelationTypeChoice()))
            f.write('bioactivity: {}\n'.format(self.getBioactivityChoice()))
            f.write('numStructures: {}\n'.format(self.numberOfStructures.get()))

        return paramsFile

    def getTargetTypeChoice(self):
        function = self.getEnumText('moleculetypeChoice')
        return function

    def getAssayTypeChoice(self):
        function = self.getEnumText('assaytypeChoice')
        return function

    def getRelationTypeChoice(self):
        function = self.getEnumText('relationChoice')
        return function

    def getBioactivityChoice(self):
        function = self.getEnumText('bioactivityChoice')
        return function

    def createOutput(self):
        file_path = self._getPath("output.smi")
        return file_path

