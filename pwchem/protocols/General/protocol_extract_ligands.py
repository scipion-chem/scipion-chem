# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, re, csv, glob, json
import numpy as np
from urllib.request import urlopen
from Bio.PDB import PDBIO, Select
from sklearn.cluster import DBSCAN

import pwem.convert as emconv

from pyworkflow.utils.path import copyFile
from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem.constants import elements_mass

from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes
import pyworkflow.object as pwobj
from pwchem.utils import clean_PDB

class ResSelect(Select):
    def __init__(self, tResidue):
        super().__init__()
        self.targetResidue = tResidue

    def accept_residue(self, residue):
        if residue == self.targetResidue:
            return 1
        else:
            return 0

class ProtocolLigandsExtraction(EMProtocol):
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
        self._insertFunctionStep('fetchStep')
        self._insertFunctionStep('createOutputStep')

    def fetchStep(self):
        # Fetch UNIPROT info
        uniprot_id = self.inputUniprotID.get().strip()
        url = "http://www.uniprot.org/uniprot/%s.json"
        with urlopen(url % uniprot_id.strip()) as response:
            jsonDic = json.loads(response.read())

        # Extract PDB Ids
        pdbIds = []
        for cross in jsonDic['uniProtKBCrossReferences']:
            if cross['database'] == 'PDB':
                pdbIds.append(cross['id'])

        # Fetch PDB info
        ligNames, alignedFns = {}, {}
        for pdbId in pdbIds:
            ligNames[pdbId] = []
            url2 = "https://data.rcsb.org/rest/v1/core/entry/%s"
            with urlopen(url2 % pdbId) as response:
                jDic = json.loads(response.read())

            # Extract Ligands Ids
            ligIds = jDic['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
            for ligId in ligIds:
                url3 = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s" % (pdbId, ligId)
                with urlopen(url3) as response:
                    jLigDic = json.loads(response.read())
                    ligNames[pdbId].append(jLigDic['pdbx_entity_nonpoly']['comp_id'])

        # Align structures
        oriASH = emconv.AtomicStructHandler()
        atomStructPath = oriASH.readFromPDBDatabase(pdbIds[0], type='mmCif', dir=self._getExtraPath())
        oriASH.read(atomStructPath)

        alignedFns[pdbIds[0]] = atomStructPath
        tmpASH = emconv.AtomicStructHandler()
        for pdbId2 in pdbIds[1:]:
            pdbFile = tmpASH.readFromPDBDatabase(pdbId2, type='mmCif', dir=self._getTmpPath())
            outFn = os.path.abspath(self._getExtraPath(os.path.basename(pdbFile)))
            matrix, rms = oriASH.getTransformMatrix(pdbFile, outFn=outFn)
            alignedFns[pdbId2] = outFn

        # Save ligands
        ligandFiles = {}
        for pdbId in ligNames:
            s = MMCIFParser().get_structure(pdbId, alignedFns[pdbId])
            io = PDBIO()
            io.set_structure(s)
            for ligId in ligNames[pdbId]:
                ligi = 0
                for residue in s.get_residues():
                    # Several HETATM residues with same name might be found. Stored in different structROIs
                    if residue.get_resname() == ligId:
                        allLigId = 'g{}_{}_{}'.format(ligi, pdbId, ligId)
                        ligandFiles[allLigId] = self._getPath(allLigId + '.pdb')
                        io.save(ligandFiles[allLigId], ResSelect(residue))
                        ligi += 1


    def createOutputStep(self):
        cmass = []
        ligandFiles = glob.glob(self._getPath('*.pdb'))
        for ligFile in ligandFiles:
            ligId = os.path.basename(ligFile).split('.')[0]
            s = PDBParser().get_structure(ligId, ligFile)
            coords, ws = [], []
            print(ligFile)
            for i, residue in enumerate(s.get_residues()):
                for atom in residue:
                    coords.append(atom.get_coord())
                    ws.append(elements_mass[atom.element])
            cmass.append(np.average(coords, axis=0, weights=ws))

        model = DBSCAN(min_samples=1, eps=5)
        pred = model.fit_predict(np.array(cmass))

        clusters = [[] for i in range(max(pred)+1)]
        for i, ligFile in enumerate(ligandFiles):
            clusters[pred[i]] += [ligFile]

        clean_PDB(self.getRefProteinFile(), outFn=self.getOriginalReceptorFile(), waters=True, HETATM=True)

        for i, cluster in enumerate(clusters):
            outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),
                                                                suffix='outputSmallMolecules_{}'.format(i))
            for fnSmall in cluster:
                smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
                smallMolecule.setPoseFile(fnSmall)
                outputSmallMolecules.append(smallMolecule)
            outputSmallMolecules.setProteinFile(self.getOriginalReceptorFile())
            outputSmallMolecules.setDocked()
            self._defineOutputs(**{'outputSmallMolecules_{}'.format(i): outputSmallMolecules})

    # --------------- INFO functions -------------------------

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

    def getRefProteinFile(self):
        return glob.glob(self._getExtraPath('*.cif'))[0]

    def getOriginalReceptorFile(self):
        return self._getPath(self.inputUniprotID.get() + '.pdb')



