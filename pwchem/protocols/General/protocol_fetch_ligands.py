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

class ProtocolLigandsFetching(EMProtocol):
    """
    Perform the extraction of ligands that bind to a target protein in PDB database

    """
    # Nombre del protocolo (aparece en grande arriba)
    _label = 'Ligand fetching'
    _experimentalMethodChoices = ['Any', 'X-RAY DIFFRACTION', 'SOLUTION NMR', 'ELECTRON MICROSCOPY',
                                  'NEUTRON DIFRACTION', 'ELECTRON CRYSTALLOGRAPHY', 'SOLID-STATE NMR',
                                  'SOLUTION SCATTERING', 'THEORETICAL MODEL']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        group = form.addGroup('Input')
        group.addParam('idType', params.EnumParam,
                       label="ID database: ", choices=['UniProt', 'PDB'], default=1,
                       help='Select database the input ID is coming from. If Uniprot, every related PDB will be fetch')
        group.addParam('inputID', params.StringParam,
                       allowsNull=False, label="Input ID: ",
                       help='Select the input ID of the target molecule')

        group = form.addGroup('Clustering Descriptors')
        group.addParam('eps', params.FloatParam, default=0, label='Maximum distance between instances: ',
                       help="The Maximum distance between two samples for one to be considered as in the neighborhood "
                            "of the other.")

        group.addParam('min_samples', params.IntParam, default=1, label='Minimum samples of a cluster: ',
                       help="The number of samples in a neighborhood for a point to be considered as a core point.")

        form.addSection(label='Filtering')
        group = form.addGroup('Filters structure')
        group.addParam('experimentalMethodChoice', params.EnumParam, default=0, label='Experimental method: ',
                       choices=self._experimentalMethodChoices,
                       help="Experimental method used to determine the protein structure")

        group.addParam('maxResolution', params.FloatParam, default=-1, label='Maximum resolution : ',
                       help="The lower the resolution value, the higher is the quality of the structure."
                            "Use -1 (or any negative) for 'Any resolution'")

        group.addParam('dateBefore', params.EnumParam, default=0,
                       label='Date filtering : ', choices=['None', 'Before', 'After'],
                       help="Only consider structures that were deposited before / after given year (included)")

        group.addParam('date', params.IntParam, default=2020,
                       label='Year : ', condition='dateBefore!=0',
                       help="Year for filtering")

        group = form.addGroup('Filters ligand')
        group.addParam('molecularWeight', params.FloatParam, default=0, label='Minimum molecular weight (Daltons): ',
                       help="Ligand minimum molecular weight to be considered in Daltons")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fetchStep')
        self._insertFunctionStep('createOutputStep')

    def fetchStep(self):
        if self.idType.get() == 0:
            # Fetch UNIPROT info
            uniprot_id = self.inputID.get().strip()
            url = "http://www.uniprot.org/uniprot/%s.json"
            with urlopen(url % uniprot_id.strip()) as response:
                jsonDic = json.loads(response.read())

            # Extract PDB Ids
            pdbIds = []
            for cross in jsonDic['uniProtKBCrossReferences']:
                if cross['database'] == 'PDB':
                    pdbIds.append(cross['id'])

        elif self.idType.get() == 1:
            pdbIds = [self.inputID.get().strip()]

        if len(pdbIds) > 0:
            # Fetch PDB info
            ligNames, alignedFns = {}, {}
            for pdbId in pdbIds:
                ligNames[pdbId] = []
                url2 = "https://data.rcsb.org/rest/v1/core/entry/%s"
                with urlopen(url2 % pdbId) as response:
                    jDic = json.loads(response.read())

                if self.checkStructureFilters(jDic):
                    # Extract Ligands Ids
                    ligIds = jDic['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
                    for ligId in ligIds:
                        url3 = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s" % (pdbId, ligId)
                        with urlopen(url3) as response:
                            jLigDic = json.loads(response.read())

                        if self.checkLigandFilters(jLigDic):
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
                    ligi = 1
                    for residue in s.get_residues():
                        # Several HETATM residues with same name might be found. Stored in different structROIs
                        if residue.get_resname() == ligId:
                            allLigId = '{}_{}_frag{}'.format(pdbId, ligId, ligi)
                            ligandFiles[allLigId] = self._getPath(allLigId + '.pdb')
                            io.save(ligandFiles[allLigId], ResSelect(residue))
                            ligi += 1
        else:
            print('No related PDB information found')


    def createOutputStep(self):
        cmass = []
        ligandFiles = glob.glob(self._getPath('*.pdb'))
        for ligFile in ligandFiles:
            ligId = os.path.basename(ligFile).split('.')[0]
            s = PDBParser().get_structure(ligId, ligFile)
            coords, ws = [], []
            for i, residue in enumerate(s.get_residues()):
                for atom in residue:
                    coords.append(atom.get_coord())
                    ws.append(elements_mass[atom.element])
            cmass.append(np.average(coords, axis=0, weights=ws))

        if len(cmass) > 0:
            model = DBSCAN(min_samples=1, eps=5)
            pred = model.fit_predict(np.array(cmass))

            clusters = [[] for i in range(max(pred)+1)]
            for i, ligFile in enumerate(ligandFiles):
                clusters[pred[i]] += [ligFile]

            clean_PDB(self.getRefProteinFile(), outFn=self.getOriginalReceptorFile(), waters=True, HETATM=True)

            outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),
                                                                suffix='outputSmallMolecules')
            for i, cluster in enumerate(clusters):
                for fnSmall in cluster:
                    smallMolecule = SmallMolecule(smallMolFilename=fnSmall)
                    smallMolecule.setMolName(os.path.basename(fnSmall).split('_frag')[0])
                    smallMolecule.setPoseId(1)
                    smallMolecule.setPoseFile(fnSmall)
                    smallMolecule.setGridId(i+1)
                    smallMolecule.setDockId(self.getObjId())
                    outputSmallMolecules.append(smallMolecule)
                outputSmallMolecules.setProteinFile(self.getOriginalReceptorFile())
                outputSmallMolecules.setDocked()
            self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    # --------------------------- UTILS functions -----------------------------------

    def getExperimentalMethodChoice(self):
        function = self.getEnumText('experimentalMethodChoice')
        return function

    def getRefProteinFile(self):
        return glob.glob(self._getExtraPath('*.cif'))[0]

    def getOriginalReceptorFile(self):
        return self._getPath(self.inputID.get() + '.pdb')

    def checkStructureFilters(self, jDic):
        method = self.getEnumText('experimentalMethodChoice')
        if method != 'Any':
            check1 = False
            for item in jDic['exptl']:
                if item['method'] == method:
                    check1 = True
        else:
            check1 = True

        maxRes = self.maxResolution.get()
        if maxRes > 0:
            check2 = False
            for res in jDic['rcsb_entry_info']['resolution_combined']:
                if res < maxRes:
                    check2 = True
        else:
            check2 = True

        if self.dateBefore.get() != 0:
            check3 = False
            date = jDic['rcsb_accession_info']['deposit_date']
            year = int(date.split('-')[0])
            if year <= self.date.get() and self.dateBefore.get() == 1:
                check3 = True
            elif year >= self.date.get() and self.dateBefore.get() == 2:
                check3 = True
        else:
            check3 = True

        return check3 and check2 and check1

    def checkLigandFilters(self, jDic):
        weight = self.molecularWeight.get()
        check1 = False
        if jDic['rcsb_nonpolymer_entity']['formula_weight'] >= weight:
            check1 = True
        return check1

