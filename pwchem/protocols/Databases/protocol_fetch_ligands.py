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
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import *
from pwchem.constants import elements_mass

from pwchem import Plugin as pwchemPlugin
from pwchem.utils import fillEmptyAttributes
import pyworkflow.object as pwobj
from pwchem.utils import clean_PDB

PDB, CHEMBL, BINDINGDB = 0, 1, 2
RDKIT, OPENBABEL = 0, 1

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

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')

        group = form.addGroup('Input')
        group.addParam('inputType', params.EnumParam, label="Get ligand from: ",
                       choices=['Uniprot entry', 'Target ID', 'Ligand ID'], default=0,
                       help='Select the type of input.\nUniprot entry: it will retrieve all targets related to the '
                            'entry in the selected database and then, all ligand related to those targets.\n'
                            'Target ID: retrieve the related ligand of the target from the selected database.\n'
                            'Ligand ID: retrieve the specified ligand from the selected database')

        group.addParam('inputDatabaseUniprot', params.EnumParam, label="Uniprot cross database: ",
                       choices=['PDB', 'ChEMBL', 'BindingDB'], default=0, condition='inputType == 0',
                       help='Select database for fetching of the target. \nIn the case of PDB, the small molecules '
                            'will be fetched together with a representant of the target structure (docked).\n'
                            'BioPython will be used to try superimposing the structures but it might lead to problems '
                            'in the position of the ligands')

        group.addParam('inputDatabaseTarget', params.EnumParam, label="Target database: ",
                       choices=['PDB', 'ChEMBL'], default=0, condition='inputType == 1',
                       help='Select database for fetching of the target. \nIn the case of PDB, the small molecules '
                            'will be fetched together with a representant of the target structure (docked).\n'
                            'BioPython will be used to try superimposing the structures but it might lead to problems '
                            'in the position of the ligands')

        group.addParam('inputDatabaseLigand', params.EnumParam, label="Ligand database: ",
                       choices=['ChEMBL', 'BindingDB', 'ZINC', 'PubChem'], default=0, condition='inputType in [2]',
                       help='Select database for fetching of the ligand')

        group.addParam('uniprotID', params.StringParam, label="Uniprot entry ID: ", default='',
                       condition='inputType == 0',
                       help='Uniprot ID to look for related targets in the selected database and then fetch ligands '
                            'on those targets')
        
        group.addParam('targetID', params.StringParam, label="Target database ID: ", default='',
                       condition='inputType == 1',
                       help='Target ID in the selected database. '
                            'Ligands associated with the target ID will be fetched if uniprot ID is empty')

        group.addParam('ligandID', params.StringParam, label="Ligand database ID: ", default='',
                       condition='inputType == 2',
                       help='Ligand ID in the selected database to be fetched')

        group = form.addGroup('Managing')
        group.addParam('doAlignStructures', params.BooleanParam, label="Download target structure: ",
                       condition='inputType == 0 and inputDatabaseUniprot==0', default=False,
                       help='Try aligning the target structures and outputing the ligands docked to one of those '
                            'aligned structures')

        group.addParam('structDatabase', params.EnumParam, label="Get ligand structure from: ",
                       choices=['PubChem', 'DrugBank', 'Smiles'], default=0,
                       condition='(inputType == 0 and inputDatabaseUniprot == 2) or '
                                 '(inputType == 2 and inputDatabaseLigand == 1)',
                       help='BindingDB does not provide directly the structures, but their IDs can be mapped to other '
                            'databases.\n'
                            'Select database for getting the ligand structures, mapping from BindingDB')

        group.addParam('useManager', params.EnumParam, default=0, label='Manage structure using: ',
                       choices=['RDKit', 'OpenBabel'],
                       help='Whether to manage the structure (parse and optimize if needed) using RDKit or OpenBabel')

        group = form.addGroup('Clustering Descriptors',
                              condition='inputType in [0, 1] and inputDatabaseTarget==0 and doAlignStructures')
        group.addParam('eps', params.FloatParam, default=0, label='Maximum distance between instances: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="The maximum distance between two samples for one to be considered as in the neighborhood "
                            "of the other. Molecules clustered in structural ROIs as in docking")

        group.addParam('min_samples', params.IntParam, default=1, label='Minimum samples of a cluster: ',
                       expertLevel=params.LEVEL_ADVANCED,
                       help="The number of samples in a neighborhood for a point to be considered as a core point."
                            "Molecules clustered in structural ROIs as in docking")

        form.addSection(label='Filtering')
        group = form.addGroup('Filters PDB', condition='inputType in [0, 1] and inputDatabaseTarget==0')
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

        group = form.addGroup('Filters ChEMBL', condition='inputType in [0, 1] and inputDatabaseTarget==1')
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

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('fetchStep')
        self._insertFunctionStep('createOutputStep')

    def fetchStep(self):
        mapDBD = False
        targetIds, ligandsDic = [], {}
        if self.inputType.get() == 0:
            # Download from Uniprot ID the target IDs or save ligands (DBD)
            iBase = self.inputDatabaseUniprot.get()
            uniprot_id = self.uniprotID.get().strip()
            self.addRelationToFile('Uniprot', uniprot_id)

            if iBase in [0, 1]:
                # Download related target IDs from PDB or ChEMBL
                jsonDic = self.getJDic('Uniprot', None, uniprot_id)
                targetIds = self.mapUniprot2TargetIds(jsonDic)
            else:
                # Download related ligands from BindingDB
                ligDic = self.mapUniprot2SmilesDic(uniprot_id)
                if self.structDatabase.get() == 2:
                    # Formatting DBD ligands from SMI
                    self.formatSMIs(ligDic, database='\t\tBindingDB')
                else:
                    # Mapping DBD ligands to other DB and download
                    mapDic = self.getDBDMapDic()
                    self.saveMappedDBDLigands(ligDic, mapDic)

        elif self.inputType.get() == 1:
            # Listing Target Ids
            targetIds = [self.targetID.get().strip()]

        elif self.inputType.get() == 2:
            # Listing and downloading ligand IDs
            ligandsDic = {self.ligandID.get().strip(): ''}
            iBase = self.inputDatabaseLigand.get()

            if iBase == 0:
                ligandNames = self.fetchFromLigandChemBL(ligandsDic)
                self.formatSMIs(self.getSMILigands(ligandNames), database='\t\tChEMBL')

            elif iBase == 1 and self.structDatabase.get() != 2:
                mapDic = self.getDBDMapDic()
                self.saveMappedDBDLigands(ligandsDic, mapDic)

            elif iBase == 1 and self.structDatabase.get() == 2:
                for ligId in ligandsDic:
                    smi = self.getDBDSmiles(ligId)
                    ligandsDic[ligId] = smi
                self.formatSMIs(ligandsDic, database='\t\tBindingDB')

            elif iBase == 2:
                ligandFiles = self.saveZINCLigands(ligandsDic)

            elif iBase == 3:
                ligandFiles = self.savePubChemLigands(ligandsDic)


        if len(targetIds) > 0:
            iBase = self.inputDatabaseTarget.get() if self.inputType.get() == 1 else iBase
            ligandNames = self.fetchLigandsFromTargets(targetIds, iBase)

            if iBase == 0:
                alignedFns = self.alignStructures(targetIds, align=self.doAlignStructures.get())
                self.savePDBLigands(ligandNames, alignedFns)
            elif iBase == 1:
                self.formatSMIs(self.getSMILigands(ligandNames), database=None)


    def createOutputStep(self):
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(),
                                                            suffix='outputSmallMolecules')
        if self.inputType.get() == 0 and self.inputDatabaseUniprot.get() == 0 and self.doAlignStructures.get():
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
            for fnSmall in glob.glob(self._getExtraPath('*')):
                if os.path.isfile(fnSmall) and not fnSmall.endswith('.txt'):
                    smallMolecule = SmallMolecule(smallMolFilename=fnSmall, molName='guess')
                    outputSmallMolecules.append(smallMolecule)

        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    # --------------------------- Summary functions -----------------------------------

    def _summary(self):
        summ = ''
        sumFile = self.getIDsRelationsFile()
        if os.path.exists(sumFile):
            with open(sumFile) as f:
                summ = f.read()

        err = ''
        errFile = self.getErrorsPath()
        if os.path.exists(errFile):
            with open(errFile) as f:
                err = f.read()
        return [summ, err]

    # --------------------------- UTILS functions -----------------------------------
    def getIDsRelationsFile(self):
        return self._getPath('idsRelations.txt')

    def addRelationToFile(self, name, id, name2=None, id2=None):
        mode = 'w'
        if os.path.exists(self.getIDsRelationsFile()):
            mode = 'a'

        with open(self.getIDsRelationsFile(), mode) as f:
            f.write('{} ID: {}'.format(name, id))
            if name2 and id2:
                f.write(' -> {} ID: {}'.format(name2, id2))
            f.write('\n')


    def mapUniprot2TargetIds(self, jsonDic):
        targetIds = []
        # Extract target Ids from PDB, ChEMBL
        for cross in jsonDic['uniProtKBCrossReferences']:
            if cross['database'] == self.getEnumText('inputDatabaseUniprot'):
                targetIds.append(cross['id'])

        return list(set(targetIds))

    def mapUniprot2SmilesDic(self, uniprot_id):
        ligDic = {}
        url = 'https://bindingdb.org/axis2/services/BDBService/getLigandsByUniprot?uniprot={}'.format(uniprot_id)
        with urlopen(url) as response:
            fullXML = response.read().decode('utf-8')
            ligIds = re.findall(r'<bdb:monomerid>\d+</bdb:monomerid>', fullXML)
            ligIds = [ligId.split('>')[1].split('<')[0] for ligId in ligIds]

            smiles = re.findall(r'<bdb:smiles>.+?</bdb:smiles>', fullXML)
            smiles = [smi.split('>')[1].split('<')[0].split()[0] for smi in smiles]

            for ligId, smi in zip(ligIds, smiles):
                ligDic[ligId] = smi

        return ligDic


    def fetchFromLigandChemBL(self, ligandsDic):
        ligNames = {}
        for ligandId in ligandsDic:
            jDic = self.getJDic('ChEMBL', 'molecule', ligandId)
            ligNames['noTarget'] = {ligandId: jDic['molecules'][0]['molecule_structures']['canonical_smiles']}
        return ligNames

    def fetchLigandsFromTargets(self, targetIds, iBase):
        if iBase == 0:
            # Fetch PDB info
            ligNames = {}
            for pdbId in targetIds:
                self.addRelationToFile('\tPDB complex', pdbId)
                jDic = self.getJDic('PDB', 'entry', pdbId)

                if self.checkStructureFilters(jDic, iBase):
                    # Extract Ligands Ids
                    if 'non_polymer_entity_ids' in jDic['rcsb_entry_container_identifiers']:
                        ligIds = jDic['rcsb_entry_container_identifiers']['non_polymer_entity_ids']
                        for ligId in ligIds:
                            jLigDic = self.getJDic('PDB', 'nonpolymer_entity', pdbId, id2=ligId)

                            if self.checkLigandFilters(jLigDic, iBase):
                                compId = jLigDic['pdbx_entity_nonpoly']['comp_id']
                                self.addRelationToFile('\t\tPDB compound', compId)
                                if not pdbId in ligNames:
                                    ligNames[pdbId] = {compId: ''}
                                else:
                                    ligNames[pdbId][compId] = ''

        elif iBase == 1:
            ligNames = {}
            for targetId in targetIds:
                self.addRelationToFile('\tChEMBL target', targetId)

                jDic = self.getJDic('ChEMBL', 'activity', targetId, limit=1)
                nLigs = jDic['page_meta']['total_count']
                jDic = self.getJDic('ChEMBL', 'activity', targetId, limit=nLigs)

                for jLigDic in jDic['activities']:
                    if self.checkStructureFilters(jLigDic, iBase):
                        chembl_id = jLigDic['molecule_chembl_id']
                        jMolDic = self.getJDic('ChEMBL', 'molecule', chembl_id)
                        if self.checkLigandFilters(jMolDic, iBase):
                            self.addRelationToFile('\t\tChEMBL compound', chembl_id)

                            if not targetId in ligNames:
                                ligNames[targetId] = {chembl_id: jLigDic['canonical_smiles']}
                            else:
                                ligNames[targetId][chembl_id] = jLigDic['canonical_smiles']

        return ligNames

    def clusterPDBLigands(self):
        cmass = []
        ligandFiles = glob.glob(self._getExtraPath('*.pdb'))
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
        atomStructPath = oriASH.readFromPDBDatabase(pdbIds[0], type='mmCif', dir=self._getPath())
        oriASH.read(atomStructPath)

        alignedFns[pdbIds[0]] = atomStructPath
        tmpASH = emconv.AtomicStructHandler()
        for pdbId2 in pdbIds[1:]:
            pdbFile = tmpASH.readFromPDBDatabase(pdbId2, type='mmCif', dir=self._getTmpPath())
            outFn = os.path.abspath(self._getTmpPath(os.path.basename(pdbFile)))
            if align:
                try:
                    matrix, rms = oriASH.getTransformMatrix(pdbFile, outFn=outFn)
                    alignedFns[pdbId2] = outFn
                except:
                    self.addToSummary('{} could not be aligned to {}'.format(pdbId2m, dbIds[0]))
                    alignedFns[pdbId2] = pdbFile
            else:
                alignedFns[pdbId2] = pdbFile

        return alignedFns

    def saveMappedDBDLigands(self, ligandsDic, mapDic):
        ligIds = []
        structDB = self.getEnumText('structDatabase')
        for ligId in ligandsDic:
            if ligId in mapDic:
                ligIds.append(mapDic[ligId])
                self.addRelationToFile('\t\tBindingDB', ligId, structDB, mapDic[ligId])
            else:
                self.addRelationToFile('\t\tBindingDB', ligId)
                self.addToSummary('BindingDB id {} could not be mapped to {} DB'.format(ligId, structDB))
        self.saveBDBLigands(ligIds)

    def savePDBLigands(self, ligNames, alignedFns):
        # Save ligands
        ligandFiles = {}
        ligIds = []
        for pdbId in ligNames:
            s = MMCIFParser().get_structure(pdbId, alignedFns[pdbId])
            io = PDBIO()
            io.set_structure(s)
            for ligId in ligNames[pdbId]:
                self.addRelationToFile('\t\tPDB ligand', ligId)
                ligi = 1
                for residue in s.get_residues():
                    # Several HETATM residues with same name might be found. Stored in different structROIs
                    if residue.get_resname() == ligId:
                        if len(list(residue.get_atoms())) > self.minAtoms.get():
                            if self.doAlignStructures:
                                allLigId = '{}_{}_frag{}'.format(pdbId, ligId, ligi)
                                ligandFiles[allLigId] = self._getExtraPath(allLigId + '.pdb')
                                io.save(ligandFiles[allLigId], ResSelect(residue))
                                ligi += 1
                            elif not ligId in ligIds:
                                ligandFiles[ligId] = self._getExtraPath(ligId + '.pdb')
                                io.save(ligandFiles[ligId], ResSelect(residue))
                                ligIds.append(ligId)

        return ligandFiles

    def getSMILigands(self, ligNames):
        ligandSMIs = {}
        for targetId in ligNames:
            for ligId in ligNames[targetId]:
                ligandSMIs[ligId] = ligNames[targetId][ligId]
        return ligandSMIs

    def formatSMIs(self, smiDic, database, outDir=None):
        '''smiDic:  {smiId: smi, ...}'''
        if database:
            for ligId in smiDic:
                self.addRelationToFile(database, ligId)
        
        if not outDir:
            outDir = os.path.abspath(self._getExtraPath())

        fnOut = os.path.abspath(self._getTmpPath('outputSMIs.smi'))
        with open(fnOut, 'w') as f:
            for smiId in smiDic:
                f.write('{} {}\n'.format(smiDic[smiId], smiId))

        args = ' -i "{}" -of sdf --outputDir "{}"'.format(fnOut, outDir)
        args += ' --make3D -nt {}'.format(self.numberOfThreads.get())

        # Formatting with RDKit (neccessary if they are maestro)
        if self.useManager.get() == RDKIT or fnOut.endswith(".mae") or fnOut.endswith(".maegz"):
            pwchemPlugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=outDir)

        # Formatting with OpenBabel
        elif self.useManager.get() == OPENBABEL:
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

    def saveBDBLigands(self, ligIds):
        if self.structDatabase.get() == 0:
            self.savePubChemLigands(ligIds, outDir=self._getExtraPath())
        elif self.structDatabase.get() == 1:
            self.saveDrugBankLigands(ligIds)


    def saveZINCLigands(self, ligand_ids, outDir=None):
        if not outDir:
            outDir = os.path.abspath(self._getExtraPath())

        for ligand_id in ligand_ids:
            self.addRelationToFile('\tZINC', ligand_id)
            if ligand_id.lower().startswith('zinc'):
                ligand_id = ligand_id[4:]

            zid = 'ZINC' + ligand_id.zfill(12)
            url = 'https://zinc.docking.org/substances/{}.sdf'.format(zid)
            with urlopen(url) as response:
                sdfStr = response.read().decode('utf-8')

            with open(os.path.join(outDir, zid) + '.sdf', 'w') as f:
                f.write(sdfStr)

    def savePubChemLigands(self, ligand_ids, outDir=None):
        if not outDir:
            outDir = os.path.abspath(self._getExtraPath())

        for ligand_id in ligand_ids:
            self.addRelationToFile('\tPubChem', ligand_id)
            outFile = os.path.abspath(os.path.join(outDir, '{}.sdf'.format(ligand_id)))
            try:
                self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(ligand_id), outFile),
                            cwd=self._getExtraPath())
            except:
                try:
                    self.runJob('wget', '"{}" -O {}'.format(self.getPubChemURL(ligand_id, dim=2), outFile),
                                cwd=self._getExtraPath())
                except:
                    self.addToSummary('Pubchem Compound with ID: {} could not be downloaded'.format(ligand_id))

    def saveDrugBankLigands(self, ligand_ids, outDir=None):
        if not outDir:
            outDir = os.path.abspath(self._getExtraPath())
        for ligand_id in ligand_ids:
            outFile = os.path.abspath(os.path.join(outDir, '{}.sdf'.format(ligand_id)))
            try:
                self.runJob('wget', '"{}" -O {}'.format(self.getDrugBankURL(ligand_id), outFile),
                            cwd=self._getExtraPath())
            except:
                self.addToSummary('DrugBank Compound with ID: {} could not be downloaded'.format(ligand_id))


    def getPubChemURL(self, pID, dim=3):
        return "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/CID/{}/record/SDF/?record_type={}d&response_type=save&response_basename=Conformer{}D_CID_{}".\
            format(pID, dim, dim, pID)

    def getDrugBankURL(self, pID):
        return 'https://go.drugbank.com/structures/small_molecule_drugs/{}.sdf'.format(pID)

    def getExperimentalMethodChoice(self):
        function = self.getEnumText('experimentalMethodChoice')
        return function

    def getRefProteinFile(self):
        cifFiles = glob.glob(self._getPath('*.cif'))
        if len(cifFiles) > 0:
            return cifFiles[0]
        else:
            return glob.glob(self._getPath('*.pdb'))[0]

    def getOriginalReceptorFile(self):
        return self._getPath(self.targetID.get() + '.pdb')

    def checkStructureFilters(self, jDic, iBase):
        checks = []
        if iBase == 0:
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

        elif iBase == 1:
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

    def countAtoms(self, molFormula, onlyHeavy=True):
        '''Count the number of atoms in a formula of the form C161H236N42O53'''
        count = 0
        letters, numbers = re.findall(r'[A-Z]+', molFormula), re.findall(r'\d+', molFormula)
        for let, num in zip(letters, numbers):
            if let == 'H' and onlyHeavy:
                continue
            else:
                count += int(num)
        return count

    def checkLigandFilters(self, jDic, iBase):
        checks = []
        weight = self.molecularWeight.get()
        minAtoms = self.minAtoms.get()

        if iBase == 0:
            checks.append(False)
            if jDic['rcsb_nonpolymer_entity']['formula_weight'] >= weight:
                checks[-1] = True
            # Number of atoms not found in jDic. Number of atoms filter present when ligand is parsed

        elif iBase == 1:
            checks.append(False)
            if float(jDic['molecules'][0]['molecule_properties']['full_mwt']) >= weight:
                checks[-1] = True

            checks.append(False)
            if 'heavy_atoms' in jDic['molecules'][0]['molecule_properties'] and \
                    jDic['molecules'][0]['molecule_properties']['heavy_atoms']:
                numAtoms = int(jDic['molecules'][0]['molecule_properties']['heavy_atoms'])
            else:
                numAtoms = self.countAtoms(jDic['molecules'][0]['molecule_properties']['full_molformula'])

            if numAtoms >= minAtoms:
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
                    url += 'limit=%s&' % limit
                url += 'target_chembl_id=%s' % id

            elif data == 'molecule':
                url = 'https://www.ebi.ac.uk/chembl/api/data/molecule.json?'
                if limit != -1:
                    url += 'limit=%s&' % limit
                url += 'molecule_chembl_id=%s' % id

        with urlopen(url) as response:
            jDic = json.loads(response.read())
        return jDic

    def getDBDSmiles(self, dbid):
        mapDic = {}
        url = 'https://www.bindingdb.org/rwd/bind/chemsearch/marvin/MolStructure.jsp?monomerid={}'.format(dbid)
        with urlopen(url) as response:
            fullHTML = response.read().decode('utf-8')

        smiLine = re.findall(r'<b>SMILES</b> <span class="darkgray">.+?</span>', fullHTML)[0]
        return smiLine.split('"darkgray">')[1].split('<')[0]

    def getDBDMapDic(self):
        mapDic = {}
        if self.structDatabase.get() == 0:
            mapFile = self._getTmpPath('DBD2PubChem.txt')
            url = "https://www.bindingdb.org/bind/BindingDB_CID.txt"
        elif self.structDatabase.get() == 1:
            mapFile = self._getTmpPath('DBD2DrugBank.txt')
            url = "https://www.bindingdb.org/bind/BindingDB_DrugBankID.txt"

        if not os.path.exists(mapFile):
            with urlopen(url) as response:
                with open(mapFile, 'w') as f:
                    f.write(response.read().decode('utf-8'))

        with open(mapFile) as f:
            for line in f:
                mapDic[line.split()[0]] = line.split()[1]

        return mapDic

    def fromUniprotDBD(self):
        return self.inputType.get() == 0 and self.inputDatabaseUniprot.get() == 2

    def getErrorsPath(self):
        return self._getPath('errors.txt')

    def addToSummary(self, sumStr):
        sumFile = self.getErrorsPath()
        w = 'w' if not os.path.exists(sumFile) else 'a'
        sumStr = sumStr if sumStr.endswith('\n') else sumStr + '\n'

        with open(sumFile, w) as f:
            f.write(sumStr)