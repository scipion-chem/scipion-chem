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
from pwchem.utils import clean_PDB

scriptName = 'ligand_filtering.py'


class ProtocolLigandFiltering(EMProtocol):
    """
    Perform the extraction of ligands that bind to a target protein form PDB database

    """
    # Nombre del protocolo (aparece en grande arriba)
    _label = 'Ligand filtering'
    _experimentalMethodChoices = ['X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY', 'NEUTRON DIFRACTION',
                                  'ELECTRON CRYSTALLOGRAPHY', 'SOLID-STATE NMR', 'SOLUTION SCATTERING',
                                  'THEORETICAL MODEL', 'None']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        form.addParam('inputUniprotID', params.StringParam,
                      allowsNull=False, label="Uniprot ID: ",
                      help='Select the Uniprot ID of the target molecule')

        group = form.addGroup('Descriptors')
        group.addParam('experimentalMethodChoice', params.EnumParam, default=0, label='Experimental method: ',
                       choices=self._experimentalMethodChoices,
                       help="Experimental method used to determine the protein structure")

        group.addParam('molecularWeight', params.FloatParam, default=0, label='Minimum molecular weight: ',
                       help="Ligand minimum molecular weight")

        group.addParam('maxResolution', params.FloatParam, default=0, label='Quality of the structure : ',
                       help="The lower the resolution value, the higher is the quality of the structure")

        group.addParam('polymerCount', params.IntParam, default=1, label='Number of chains : ',
                       help="Number of structure chains")

        group.addParam('date', params.StringParam, default=2020, allowsNull=False,
                       label='Structures deposited before : ',
                       help="Only consider structures that were deposited before given date. The user must choose a year and the date that will be passed to the search criteria is January 1st of the year")

        group.addParam('numberOfStructures', params.IntParam, default=2, label='Number of top structures : ',
                       help="Number of top structures chosen after the pdb filtering")

        group = form.addGroup('Clustering Descriptors')

        group.addParam('eps', params.FloatParam, default=0, label='Maximum distance between instances: ',
                       help="The Maximum distance between two samples for one to be considered as in the neighborhood of the other.")

        group.addParam('min_samples', params.IntParam, default=1, label='Minimum samples of a cluster: ',
                       help="The number of samples in a neighborhood for a point to be considered as a core point.")

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

    def importStep(self):
        files_path = self.parseResults(self._getPath("output.txt"))
        protein = self.obtain_protein_file((self._getPath("protein.txt")))
        smiles_dict = self.parse_csv(self._getPath("PDB_top_ligands.csv"))
        clus = 0
        for cluster in files_path:
            clus += 1
            for filename in cluster:
                fnSmall = self._getExtraPath(os.path.split(filename)[1])

                copyFile(filename, fnSmall)
            suffix = '_{}'.format(clus)
            outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='consensusSmallMolecules{}'.format(suffix))

            for fnSmall in cluster:
                smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
                if len(os.listdir(self._getExtraPath())) <= 100:
                    if not fnSmall.endswith('.mae') or not fnSmall.endswith('.maegz'):
                        fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                        fnOut = self._getExtraPath("%s.png" % fnRoot)
                        args = Plugin.getPluginHome('utils/rdkitUtils.py') + " draw %s %s" % (fnSmall, fnOut)

                        try:
                            Plugin.runRDKit(self, "python3", args)
                            smallMolecule._PDBLigandImage = pwobj.String(fnOut)
                        except:
                            smallMolecule._PDBLigandImage = pwobj.String("Not available")
                key2 = fnSmall.split("/")[-1]
                smallMolecule._Smiles = pwobj.String(smiles_dict[key2])
                outputSmallMolecules.append(smallMolecule)
            outputSmallMolecules.setProteinFile(protein)
            self._defineOutputs(**{'consensusSmallMolecules{}'.format(suffix): outputSmallMolecules})

    # --------------- INFO functions -------------------------

    def _citations(self):
        return [ "@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods = """Filtering choices:\n
        - Uniprot ID of the target protein.\n
        - Experimental method: structure determination method, among the possible options are X-RAY DIFFRACTION, SOLUTION NMR, ELECTRON
        MICROSCOPY, etc or None, in case the experimental technique is not specified.\n
        - Minimum molecular weight of the ligand associated with the structure.\n
        - Quality of the structure: the lower the angstrom resolution of the structure, the higher the quality.\n
        - Number of structure chains: it is recommended to set the value of this field to 1, as it simplifies the processing of the structure and the search.\n
        - Date of deposition: this field is enabled to ensure the reproducibility of the screening; the user sets a date after which any structure that is uploaded on the database will not be considered.\n
        - Number of final extracted ligands."""

        return methods

    # --------------------------- UTILS functions -----------------------------------

    def writeParamsFile(self, paramsFile):
        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('uniprotID: {}\n'.format(self.inputUniprotID.get()))
            f.write('experimentalMethod: {}\n'.format(self.getExperimentalMethodChoice()))
            f.write('molecularWeigth: {}\n'.format(self.molecularWeight.get()))
            f.write('maxResolution: {}\n'.format(self.maxResolution.get()))
            f.write('polymerCount: {}\n'.format(self.polymerCount.get()))
            f.write('topStructures: {}\n'.format(self.numberOfStructures.get()))
            f.write('date: {}\n'.format(self.date.get()))
            f.write('eps: {}\n'.format(self.eps.get()))
            f.write('min_samples: {}\n'.format(self.min_samples.get()))
        return paramsFile

    def getExperimentalMethodChoice(self):
        function = self.getEnumText('experimentalMethodChoice')
        return function

    def parseResults(self, outputFile):
        paths = []
        with open(outputFile) as tsv_file:
            list_cluster = []
            for row in tsv_file:
                row = row.replace('\n', '')
                if "Cluster" in row:
                    list_cluster = []
                elif "#" in row:
                    if list_cluster:
                        paths.append(list_cluster)
                else:
                    list_cluster.append(row)

        return paths


    def obtain_protein_file(self, pdb_file):
        with open(pdb_file) as tsv_file:
            for row in tsv_file:
                pdb_file1 = row
                #str(row).split("/")
                #pdb = pdb_file1[-1]
                out = "only_protein.pdb"
                out_protein = clean_PDB(pdb_file1, out, waters=True, HETATM=True)

        return out_protein


    def parse_csv(self, csv_file):
        dict_smile = {}
        with open(csv_file) as csv:
            for row in csv:
                if row[0] == "#":
                    pass
                else:
                    word = row.split(",")
                    key = str(word[0]) + "_lig.pdb"
                    dict_smile[key] = word[-1]

        return dict_smile



