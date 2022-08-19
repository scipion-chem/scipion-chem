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

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import fillEmptyAttributes
import pyworkflow.object as pwobj
from pwchem.utils import clean_PDB

PDB, CHEMBL = 0, 1

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
    _assaytypeChoices = ['Any', '(B) Binding', '(F) Functional', '(A) ADMET', '(T) Toxicity', '(P) Physicochemical',
                         '(U) Unclassified']
    _bioactivityChoices = ['Any', 'Ki', 'Kd', 'IC50', 'EC50']



    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        group = form.addGroup('Input')
        group.addParam('uniprotID', params.StringParam, label="Uniprot ID: ", default='',
                       help='Uniprot ID to look for related ligands in the selected database')
        
        group.addParam('inputDatabase', params.EnumParam, label="Fetch from database: ", 
                       choices=['PDB', 'ChEMBL', 'ZINC', 'PubChem', 'BindingDB'], default=0,
                       help='Select database for fetching. \nIn the case of PDB, the small molecules will be fetched '
                            'together with a representant of the target structure (docked). BioPython will be used '
                            'to try superimposing the structures but it might lead to problems in the position of the '
                            'ligands')
        
        group.addParam('targetID', params.StringParam, label="Target database ID: ", default='',
                       help='Target ID in the selected database. '
                            'Ligands associated with the target ID will be fetched if uniprot ID is empty')

        group.addParam('doAlignStructures', params.BooleanParam, label="Download target structure: ",
                       condition='inputDatabase==0', default=False,
                       help='Try aligning the target structures and outputing the ligands docked to one of those '
                            'aligned structures')
        
        group.addParam('ligandID', params.StringParam, label="Ligand database ID: ",
                       condition='inputDatabase!=0', default='',
                       help='Ligand ID in the selected database to be fetched')
        

        group = form.addGroup('Clustering Descriptors')
        group.addParam('eps', params.FloatParam, default=0, label='Maximum distance between instances: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="The maximum distance between two samples for one to be considered as in the neighborhood "
                            "of the other. Molecules clustered in structural ROIs as in docking")

        group.addParam('min_samples', params.IntParam, default=1, label='Minimum samples of a cluster: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="The number of samples in a neighborhood for a point to be considered as a core point."
                            "Molecules clustered in structural ROIs as in docking")

        form.addSection(label='Filtering')
        group = form.addGroup('Filters PDB', condition='inputDatabase==0')
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

        group = form.addGroup('Filters ChEMBL', condition='inputDatabase==1')
        group.addParam('assayTypeChoice', params.EnumParam, default=0, label='Assay type: ',
                       choices=self._assaytypeChoices,
                       help="Classification of the assay searched a. "
                            "Binding (B): it is based on the binding of ligand molecules to receptors, antibodies, "
                            "or other macromolecules. This type of assays obtain data for measuring binding of compound"
                            " to a molecular target, e.g. Ki, IC50, Kd.\n b. Functional (F): a set of systematic in "
                            "vivo experiments designed to measure the effect or biological role of a compound in a "
                            "cellular pathway or biological process, e.g., cell death in a cell line.\n c. ADMET (A): "
                            "assays related to the pharmacokinetic characteristics of a compound e.g., t1/2, oral "
                            "bioavailability.\n d. Toxicity (T): Data measuring toxicity of a compound, e.g., "
                            "cytotoxicity.\n e. Physicochemical (P): this type of assay is performed in the absence "
                            "of biological material and measures physicochemical properties of the compounds e.g., "
                            "chemical stability, solubility\n f. Unclassified (U): those assays which cannot be "
                            "included in the rest of type are classified into one of the above categories "
                            "e.g., ratio of binding vs efficacy.")

        group.addParam('bioactivityChoice', params.EnumParam, default=0, label='Bioactivity type: ',
                       choices=self._bioactivityChoices,
                       help="Bioactivity measure to be considered (if reported with other, it will be skipped)"
                            "If Any, the filter will not apply\n"
                            "a. IC50: (Half maximal inhibitory concentration) indicates the "
                            "required concentration of a compound to inhibit a biological process by 50%. "
                            "The lower the IC50 value, the higher the compound bioactivity.\n "
                            "b. EC50: (Half maximal effective concentration) is the drug concentration required to "
                            "produce half (50%) of the maximum effect of that compound. The lower the EC50 value, "
                            "the higher the compound bioactivity. EC50 is equivalent to IC50, but this last is used "
                            "in the case of working with inhibitors.\n c. Kd: dissociation constant, this term is "
                            "generic and describes the binding affinity between a molecule and an enzyme or a receptor."
                            "\n d. Ki: inhibition constant, the term is equivalent to the dissociation constant (Kd), "
                            "but it is used when the molecule for which bioactivity is measured is an inhibitor.")

        group.addParam('minBioactivity', params.FloatParam, default=0, label='Minimum bioactivity (nM) : ',
                       help="Minimum value for the selected bioactivity (for each compound, several activities can "
                            "be reported. If any of them passes this threshold, the compound is output.\n"
                            "The higher the bioactivity the higher the affinity of the ligand with the target.\n")

        group = form.addGroup('Filters ligand')
        group.addParam('molecularWeight', params.FloatParam, default=0,
                       label='Minimum molecular weight (Atomic Mass Units): ',
                       help="Ligand minimum molecular weight to be considered in Atomic Mass Units")
        group.addParam('minAtoms', params.IntParam, default=0,
                       label='Minimum number of heavy atoms: ',
                       help="Minimum number of non-hydrogen atoms to consider a ligand (e.g: for discarding salts)")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fetchStep')
        self._insertFunctionStep('createOutputStep')

    def fetchStep(self):
        targetIds = []
        uniprot_id, target_id, ligand_id = self.uniprotID.get(), self.targetID.get(), self.ligandID.get()
        if uniprot_id and uniprot_id.strip():
            self.addRelationToFile('Uniprot', uniprot_id)

            # Fetch UNIPROT info
            uniprot_id = self.uniprotID.get().strip()
            jsonDic = self.getJDic('Uniprot', None, uniprot_id)

            # function mapping uniprot ID to targets (PDBs, ChEMBL targets,...)
            targetIds = self.mapUniprot2TargetIds(jsonDic)
        elif target_id and target_id.strip():
            targetIds = [self.targetID.get().strip()]

        # from targetIDs to ligands function
        print(targetIds)
        if len(targetIds) > 0:
            ligandNames = self.fetchLigandsFromTargets(targetIds)

            if self.inputDatabase.get() == PDB:
                alignedFns = self.alignStructures(targetIds, align=self.doAlignStructures.get())
                ligandFiles = self.savePDBLigands(ligandNames, alignedFns)
            elif self.inputDatabase.get() == CHEMBL:
                ligandFiles = self.saveSMILigands(ligandNames)

        elif ligand_id and ligand_id.strip():
            ligandNames = self.fetchLigandFromID(ligand_id.strip())
            if self.inputDatabase.get() == CHEMBL:
                ligandFiles = self.saveSMILigands(ligandNames)

        else:
            print('No related target information found in ', self.getEnumText('inputDatabase'))


    def createOutputStep(self):
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),
                                                            suffix='outputSmallMolecules')
        if self.inputDatabase.get() == PDB:
            if self.doAlignStructures.get():
                clusters = self.clusterPDBLigands()
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

            else:
                ligandFiles = glob.glob(self._getPath('*.pdb'))
                for fnSmall in ligandFiles:
                    smallMolecule = SmallMolecule(smallMolFilename=fnSmall, molName='guess')
                    outputSmallMolecules.append(smallMolecule)

        elif self.inputDatabase.get() == CHEMBL:
            outDir = os.path.abspath(self._getPath())
            args = ' --multiFiles -iD "{}" --pattern {} -of mol2 --outputDir "{}"'. \
                format(os.path.abspath(self._getTmpPath()), '*.smi', outDir)
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

            for fnSmall in glob.glob(self._getPath('*.mol2')):
                smallMolecule = SmallMolecule(smallMolFilename=fnSmall, molName='guess')
                outputSmallMolecules.append(smallMolecule)

        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    # --------------------------- UTILS functions -----------------------------------
    def getIDsRelationsFile(self):
        return self._getExtraPath('idsRelations.txt')

    def addRelationToFile(self, name, id):
        mode = 'w'
        if os.path.exists(self.getIDsRelationsFile()):
            mode = 'a'

        with open(self.getIDsRelationsFile(), mode) as f:
            f.write('{} ID: {}\n'.format(name, id))


    def mapUniprot2TargetIds(self, jsonDic):
        if self.inputDatabase.get() in [PDB, CHEMBL]:
            # Extract target Ids from PDB, ChEMBL
            targetIds = []
            for cross in jsonDic['uniProtKBCrossReferences']:
                if cross['database'] == self.getEnumText('inputDatabase'):
                    targetIds.append(cross['id'])

        else:
            pass
            # Getting target ChEMBL IDs (take only first target or more?)
            #url = 'https://www.ebi.ac.uk/chembl/api/data/target.json?limit=30&target_components__accession=%s'
            #with urlopen(url % uniprot_id) as response:
             #   jDic = json.loads(response.read())

        return targetIds

    def fetchLigandsFromTargets(self, targetIds):
        if self.inputDatabase.get() == PDB:
            # Fetch PDB info
            ligNames = {}
            for pdbId in targetIds:
                self.addRelationToFile('\tPDB complex', pdbId)

                jDic = self.getJDic('PDB', 'entry', pdbId)

                if self.checkStructureFilters(jDic):
                    # Extract Ligands Ids
                    if 'non_polymer_entity_ids' in jDic['rcsb_entry_container_identifiers']:
                        ligIds = jDic['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
                        for ligId in ligIds:
                            jLigDic = self.getJDic('PDB', 'nonpolymer_entity', pdbId, id2=ligId)

                            if self.checkLigandFilters(jLigDic):
                                compId = jLigDic['pdbx_entity_nonpoly']['comp_id']
                                self.addRelationToFile('\t\tPDB compound', compId)
                                if not pdbId in ligNames:
                                    ligNames[pdbId] = {compId: ''}
                                else:
                                    ligNames[pdbId][compId] = ''


        elif self.inputDatabase.get() == CHEMBL:
            ligNames = {}
            for targetId in targetIds:
                self.addRelationToFile('\tChEMBL target', targetId)

                jDic = self.getJDic('ChEMBL', 'activity', targetId, limit=1)
                nLigs = jDic['page_meta']['total_count']
                jDic = self.getJDic('ChEMBL', 'activity', targetId, limit=nLigs)

                for jLigDic in jDic['activities']:
                    if self.checkStructureFilters(jLigDic):
                        chembl_id = jLigDic['molecule_chembl_id']
                        jMolDic = self.getJDic('ChEMBL', 'molecule', chembl_id)

                        if self.checkLigandFilters(jMolDic):
                            self.addRelationToFile('\t\tChEMBL compound', chembl_id)

                            if not targetId in ligNames:
                                ligNames[targetId] = {chembl_id: jLigDic['canonical_smiles']}
                            else:
                                ligNames[targetId][chembl_id] = jLigDic['canonical_smiles']

        return ligNames

    def fetchLigandFromID(self, ligandId):
        ligNames = {}
        if self.inputDatabase.get() == CHEMBL:
            jDic = self.getJDic('ChEMBL', 'molecule', ligandId)
            ligNames['noTarget'] = {ligandId: jDic['molecules'][0]['molecule_structures']['canonical_smiles']}
        return ligNames

    def clusterPDBLigands(self):
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

            clusters = [[] for i in range(max(pred) + 1)]
            for i, ligFile in enumerate(ligandFiles):
                clusters[pred[i]] += [ligFile]

            clean_PDB(self.getRefProteinFile(), outFn=self.getOriginalReceptorFile(), waters=True, HETATM=True)
        return clusters

    def alignStructures(self, pdbIds, align):
        # Align structures
        alignedFns = {}
        oriASH = emconv.AtomicStructHandler()
        atomStructPath = oriASH.readFromPDBDatabase(pdbIds[0], type='mmCif', dir=self._getExtraPath())
        oriASH.read(atomStructPath)

        alignedFns[pdbIds[0]] = atomStructPath
        tmpASH = emconv.AtomicStructHandler()
        for pdbId2 in pdbIds[1:]:
            pdbFile = tmpASH.readFromPDBDatabase(pdbId2, type='mmCif', dir=self._getTmpPath())
            outFn = os.path.abspath(self._getExtraPath(os.path.basename(pdbFile)))
            if align:
                try:
                    matrix, rms = oriASH.getTransformMatrix(pdbFile, outFn=outFn)
                    alignedFns[pdbId2] = outFn
                except:
                    print('{} could not be aligned'.format(pdbId2))
                    alignedFns[pdbId2] = pdbFile
            else:
                alignedFns[pdbId2] = pdbFile

        return alignedFns

    def savePDBLigands(self, ligNames, alignedFns):
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
                        if len(residue.get_atoms()) > self.minAtoms.get():
                            allLigId = '{}_{}_frag{}'.format(pdbId, ligId, ligi)
                            ligandFiles[allLigId] = self._getPath(allLigId + '.pdb')
                            io.save(ligandFiles[allLigId], ResSelect(residue))
                            ligi += 1
        return ligandFiles

    def saveSMILigands(self, ligNames, outDir=None):
        # Save ligands
        if not outDir:
            outDir = os.path.abspath(self._getTmpPath())
        ligandFiles = {}
        for targetId in ligNames:
            for ligId in ligNames[targetId]:
                with open(os.path.join(outDir, ligId) + '.smi', 'w') as f:
                    f.write(ligNames[targetId][ligId])
        return ligandFiles

    def getExperimentalMethodChoice(self):
        function = self.getEnumText('experimentalMethodChoice')
        return function

    def getRefProteinFile(self):
        return glob.glob(self._getExtraPath('*.cif'))[0]

    def getOriginalReceptorFile(self):
        return self._getPath(self.targetID.get() + '.pdb')

    def checkStructureFilters(self, jDic):
        checks = []
        if self.inputDatabase.get() == PDB:
            method = self.getEnumText('experimentalMethodChoice')
            if method != 'Any':
                checks.append(False)
                for item in jDic['exptl']:
                    if item['method'] == method:
                        checks[-1] = True
            else:
                checks.append(True)

            maxRes = self.maxResolution.get()
            if maxRes > 0:
                checks.append(False)
                for res in jDic['rcsb_entry_info']['resolution_combined']:
                    if res < maxRes:
                        checks[-1] = True
            else:
                checks.append(True)

            if self.dateBefore.get() != 0:
                checks.append(False)
                date = jDic['rcsb_accession_info']['deposit_date']
                year = int(date.split('-')[0])
                if year <= self.date.get() and self.dateBefore.get() == 1:
                    checks[-1] = True
                elif year >= self.date.get() and self.dateBefore.get() == 2:
                    checks[-1] = True
            else:
                checks.append(True)

        elif self.inputDatabase.get() == CHEMBL:
            assay = self.getEnumText('assayTypeChoice')
            if assay != 'Any':
                checks.append(False)
                if jDic['assay_type'] == assay.split(')')[0].split('(')[1]:
                    checks[-1] = True
            else:
                checks.append(True)

            bioMeasure = self.getEnumText('bioactivityChoice')
            if bioMeasure != 'Any':
                checks.append(False)
                if jDic['standard_type'] == bioMeasure and jDic['standard_units'] == 'nM':
                    if float(jDic['standard_value']) >= self.minBioactivity.get():
                        checks[-1] = True
            else:
                checks.append(True)

        return all(checks)

    def checkLigandFilters(self, jDic):
        checks = []
        weight = self.molecularWeight.get()
        minAtoms = self.minAtoms.get()

        if self.inputDatabase.get() == PDB:
            checks.append(False)
            if jDic['rcsb_nonpolymer_entity']['formula_weight'] >= weight:
                checks[-1] = True
            # Number of atoms not found in jDic. Number of atoms filter present when ligand is parsed

        elif self.inputDatabase.get() == CHEMBL:
            checks.append(False)
            if float(jDic['molecules'][0]['molecule_properties']['full_mwt']) >= weight:
                checks[-1] = True

            checks.append(False)
            if float(jDic['molecules'][0]['molecule_properties']['heavy_atoms']) >= minAtoms:
                checks[-1] = True

        return all(checks)

    def getJDic(self, database, data, id, id2=None, limit=-1):
        if database == 'Uniprot':
            url = "http://www.uniprot.org/uniprot/%s.json" % id
        elif database == 'PDB':
            if data == 'entry':
                url = "https://data.rcsb.org/rest/v1/core/entry/%s" % id
            elif data == 'nonpolymer_entity':
                url = "https://data.rcsb.org/rest/v1/core/nonpolymer_entity/%s/%s" % (id, id2)

        elif database == 'ChEMBL':
            if data == 'activity':
                url = 'https://www.ebi.ac.uk/chembl/api/data/activity.json?'
                if limit != -1:
                    url += 'limit=%s&'
                url += 'target_chembl_id=%s' % id

            elif data == 'molecule':
                url = 'https://www.ebi.ac.uk/chembl/api/data/molecule.json?'
                if limit != -1:
                    url += 'limit=%s&'
                url += 'molecule_chembl_id=%s' % id

        with urlopen(url) as response:
            jDic = json.loads(response.read())
        return jDic

