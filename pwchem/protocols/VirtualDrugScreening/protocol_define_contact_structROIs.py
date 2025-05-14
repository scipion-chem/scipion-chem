# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
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


"""
This protocol is used to define structural regions from the contacts between a setofsmallmolecules and their receptor

"""
import os, pickle
import numpy as np
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import get_surface

from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import runInParallel, obabelMolConversion, performBatchThreading, clusterSurfaceCoords,\
    numberSort, createPocketFile, runOpenBabel, parseAtomStruct, isHet, getBaseName, removeNumberFromStr, \
    getMAEMoleculeFiles
from pwchem import Plugin
from pwchem.constants import MGL_DIC
from pwchem.protocols import ProtDefineStructROIs

RESIDUE, ATOM = 0, 1

class ProtDefineContactStructROIs(ProtDefineStructROIs):
    """
    Defines a set of structural ROIs ligand-receptor contacts
    """
    _label = 'Define contact structural ROIs'
    typeLabels = {'All': 'set', 'ROI': 'pocket', 'Molecule': 'molecule',
                  'Single Ligand': 'single'}

    def __init__(self, **kwargs):
        ProtDefineStructROIs.__init__(self, **kwargs)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSmallMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      label="Input Small Molecules: ",
                      help='Select the SetOfSmallMolecules you want to extract the contact ROIs from')

        group = form.addGroup('Molecules selection')
        group.addParam('selectionType', params.EnumParam,
                       choices=['All', 'ROI', 'Molecule', 'Single Ligand'], default=0,
                       label='Select a subset of ligands from: ',
                       help='Select a subset of ligands to be considered for the contacts.'
                            'All takes all ligands in the set')
        group.addParam('ligandSelection', params.StringParam,
                       condition="selectionType!=0", default='', label='Use ligands in group: ',
                       help='Only use the following subset of ligands to calculate the residue contacts')

        group = form.addGroup('Contacts definition')
        group.addParam('threshold', params.FloatParam, default=4.0, label='Distance threshold (A): ',
                      help="Distance threshold where atoms of ligand-residue will be considered in contact")
        group.addParam('contactLevelRec', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST, expertLevel=params.LEVEL_ADVANCED,
                      default=RESIDUE, label='Store receptor contacts at level: ', choices=['Residue', 'Atom'],
                      help="Whether to store the residue or the atom in contact for the receptor")
        group.addParam('contactLevelLig', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST, expertLevel=params.LEVEL_ADVANCED,
                      default=ATOM, label='Store ligand contacts at level: ', choices=['Residue', 'Atom'],
                      help="Whether to store the residue or the atom in contact for the ligand")

        group = form.addGroup('ROI definition')
        group = self._defineClusterParams(group)

        form.addParallelSection(threads=4, mpi=1)


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('convertInputStep')
        self._insertFunctionStep('getContactsStep')
        self._insertFunctionStep('definePocketsStep')
        self._insertFunctionStep('defineOutputStep')

    def convertInputStep(self):
        inMols = self.getSelectedInputMols()
        outDir = self.getInputMolsDir()
        os.mkdir(outDir)

        maeMols, otherMols = getMAEMoleculeFiles(inMols)
        if len(maeMols) > 0:
            try:
                from pwchemSchrodinger.utils.utils import convertMAEMolSet
                convertMAEMolSet(maeMols, outDir, self.numberOfThreads.get(), updateSet=False)
            except ImportError:
                print(
                    'Conversion of MAE input files could not be performed because schrodinger plugin is not installed')

        runInParallel(obabelMolConversion, '.pdb', outDir, True,
                      paramList=[item.clone() for item in otherMols], jobs=self.numberOfThreads.get())

        recFile = self.inputSmallMols.get().getProteinFile()
        self.convertReceptor2PDB(recFile)

    def performContactAnalysis(self, molFns, molLists, it, recFile, outDir):
        for ligFile in molFns:
            recAtoms, ligAtoms = self.parseRecLigAtoms(asFile=recFile, ligFile=ligFile)

            contactDic = self.getLigandContacts(recAtoms, ligAtoms)
            self.writeContactFile(contactDic, ligFile, outDir)

    def getContactsStep(self):
        recFile = self.getReceptorPDB()
        molFns = self.getInputMolFiles()

        outDir = self.getContactsDir()
        os.mkdir(outDir)

        performBatchThreading(self.performContactAnalysis, molFns, self.numberOfThreads.get(), cloneItem=False,
                              recFile=recFile, outDir=outDir)

    def definePocketsStep(self):
        conStrs, cDir = [], self.getContactsDir()
        for contFile in os.listdir(cDir):
          nConStrs, areRes = self.parseContacts(os.path.join(cDir, contFile))
          conStrs += nConStrs

        conStrs = numberSort(list(set(conStrs)))

        coordsClusters = []
        recFile = self.getReceptorPDB()
        pocketCoords = self.getContactCoords(conStrs, recFile, areRes)
        if self.surfaceCoords:
            surfStruct = self.getInputSurfaceStructure()
            pocketCoords = self.mapSurfaceCoords(pocketCoords, surfStruct)

        coordsClusters += clusterSurfaceCoords(pocketCoords, self.maxIntraDistance.get())

        with open(self.getClusterCoordsPickle(), 'wb') as f:
          pickle.dump(coordsClusters, f)

    def defineOutputStep(self):
        pdbFile = self.getReceptorPDB()
        outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))

        with open(self.getClusterCoordsPickle(), 'rb') as f:
          coordsClusters = pickle.load(f)

        for i, clust in enumerate(coordsClusters):
            pocketFile = createPocketFile(clust, i, outDir=self._getExtraPath())
            pocket = StructROI(pocketFile, pdbFile)
            pocket.calculateContacts()
            outPockets.append(pocket)

        if len(outPockets) > 0:
            outPockets.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outPockets)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
      validations = []
      pSet = self.inputSmallMols.get()
      if not pSet.isDocked():
        validations.append('Sets of input molecules must be docked first.\n'
                           'Set: {} has not been docked'.format(pSet))

      return validations

    # --------------------------- UTILS functions -----------------------------------
    # ----------- General ---------------
    def getInputMolsDir(self):
        return os.path.abspath(self._getExtraPath('inputMolecules'))

    def getReceptorPDB(self):
        outDir = self.getInputMolsDir()
        return os.path.abspath(os.path.join(outDir, 'receptor.pdb'))

    def getContactsDir(self):
        return os.path.abspath(self._getExtraPath('contacts'))

    def getClusterCoordsPickle(self):
      return self._getExtraPath('coordsCluster.pkl')

    # ----------- Convert step functions ---------------
    def getSelectedInputMols(self):
      '''Return the molecules present in the group selected in the form'''
      molSet = self.inputSmallMols.get()
      vType = self.typeLabels[self.getEnumText('selectionType')]

      groupDic = molSet.getGroupIndexes()[vType]
      if vType != 'set':
        indexes = groupDic[self.ligandSelection.get()]
      else:
        indexes = list(groupDic.values())[0]
      return molSet.getMolsFromIds(indexes)

    def convertReceptor2PDB(self, proteinFile):
        inExt = os.path.splitext(os.path.basename(proteinFile))[1]
        oFile = self.getReceptorPDB()
        if not os.path.exists(oFile):
          args = ' -i{} {} -opdb -O {}'.format(inExt[1:], os.path.abspath(proteinFile), oFile)
          runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())
        return oFile

    # ----------- Get contacts functions ---------------
    def getInputMolFiles(self):
        molFiles = []
        molDir = self.getInputMolsDir()
        for file in os.listdir(molDir):
            if 'receptor.pdb' not in file:
                molFiles.append(os.path.join(molDir, file))
        return molFiles

    def parseRecLigAtoms(self, asFile, ligFile=None, ligName=None):
        '''Return two list of the form with the atoms for the receptor and ligand or heteroatoms
        '''
        if ligFile:
            recAtoms, ligAtoms = self.parseSingleFile(asFile), self.parseSingleFile(ligFile)
        else:
            recAtoms, ligAtoms = self.parseCombinedFile(asFile, ligName=ligName)
        return recAtoms, ligAtoms

    def parseSingleFile(self, asFile):
        '''Parse all the atoms stored in a AtomStructFile'''
        atoms = []
        parser = parseAtomStruct(asFile)

        for model in parser:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms.append(atom)
        return atoms

    def parseCombinedFile(self, asFile, ligName=None):
        '''Parse a combined AtomStruct file trying to separate the receptor from the ligand atoms'''
        recAtoms, ligAtoms = [], []
        parser = parseAtomStruct(asFile)

        for model in parser:
            for residue in model.get_residues():
                if isHet(residue) and (residue.resname == ligName or not ligName):
                    for atom in residue:
                        ligAtoms.append(atom)
                else:
                    for atom in residue:
                        recAtoms.append(atom)
        return recAtoms, ligAtoms

    def getAtomsCoords(self, atoms):
        coords = {}
        for atom in atoms:
            coords[atom] = list(atom.get_coord())
        return coords

    def getLigandContacts(self, recAtoms, ligAtoms):
        '''Return a dictionary of the form:
        {ligAtomId: [recAtomIds]} with the atoms in contact
        '''
        contactAtoms = {}
        ligCoords, recCoords = self.getAtomsCoords(ligAtoms), self.getAtomsCoords(recAtoms)

        dists = distance.cdist(list(recCoords.values()), list(ligCoords.values()))
        idxs = np.argwhere(dists <= self.threshold.get())
        for idx in idxs:
            recAtom, ligAtom = list(recCoords.keys())[idx[0]], list(ligCoords.keys())[idx[1]]
            if ligAtom in contactAtoms:
                contactAtoms[ligAtom].append(recAtom)
            else:
                contactAtoms[ligAtom] = [recAtom]
        return contactAtoms

    def writeContactFile(self, contactDic, ligFile, outDir):
        ligLvlStr = 'RESIDUE' if self.contactLevelLig.get() == RESIDUE else 'ATOM'
        recLvlStr = 'RESIDUE' if self.contactLevelRec.get() == RESIDUE else 'ATOM'

        outFile = self.getContactFileName(ligFile, outDir)
        with open(outFile, 'w') as fOut:
            fOut.write('## RESIDUE representation: ResidueName.ResidueNumber(Number of residue in sequence)\n')
            fOut.write('## ATOM representation: AtomName.AtomSerialNumber(Number of atom in whole structure)\n')
            fOut.write(f'LIGFILE:: {os.path.abspath(ligFile)}\n')

            # Get ligand representations of atoms/residues and corresponding Atoms Bio objects
            ligDic = self.getLigEquivalenceDic(contactDic)
            ligStrs = numberSort(list(ligDic.keys()))

            for ligRep in ligStrs:
                # For each ligand rep, get receptor contact representantions
                recStrs = self.getRecContactReps(contactDic, ligDic, ligRep)

                fOut.write(f'LIGAND {ligLvlStr} :: {ligRep}\n')
                fOut.write(f'\tRECEPTOR {recLvlStr} :: {",".join(recStrs)}\n')
        return outFile

    def getContactFileName(self, ligFile, outDir):
        return os.path.join(outDir, getBaseName(ligFile) + '_contacts.txt')

    def getLigEquivalenceDic(self, cDic):
        '''Return a dictionary with {ligRep: [ligAtoms]}
        - ligRep: Str representing atom/residue of the ligand
        - ligAtom: BioPython Atom object contained by the ligRep
        '''
        ligDic = {}
        for ligAtom in cDic:
            ligRep = self.getNameFromAtom(ligAtom, self.contactLevelLig.get())
            if ligRep in ligDic:
                ligDic[ligRep].append(ligAtom)
            else:
                ligDic[ligRep] = [ligAtom]
        return ligDic

    def getNameFromAtom(self, atom, contactLevel=RESIDUE):
        if contactLevel == ATOM:
            name = f'{removeNumberFromStr(atom.name)}{atom.serial_number}'
        else:
            name = f'{atom.get_parent().resname}{atom.get_parent().id[1]}'
        return name

    def getRecContactReps(self, cDic, ligDic, ligRep):
        '''Given a ligRep, return all the recRep of contact receptor representant.
        - cDic: dic containing contact info {ligAtom: recAtom}
        - ligDic: dic containing ligand representations {ligRep: [ligAtoms]}
        - ligRep: ligand representation to check contacts
        '''
        recReps = []
        for ligAtom in ligDic[ligRep]:
            for recAtom in cDic[ligAtom]:
                recRep = self.getNameFromAtom(recAtom, self.contactLevelRec.get())
                recReps.append(recRep)
        recReps = numberSort(set(recReps))
        return recReps


    # ----------- Build ROIs functions ---------------
    def parseContacts(self, contactsFile):
      '''Return the contact string identifiers stored in the contacts file'''
      contacts, areResidues = [], True
      with open(contactsFile) as f:
        for line in f:
          if line.startswith('\tRECEPTOR'):
            cs = line.split('::')[1].strip().split(',')
            contacts += cs
            areResidues = line.split('::')[0].split()[-1].strip() == 'RESIDUE'

      return set(contacts), areResidues

    def getContactCoords(self, conStrs, recFile, areResidues=True):
      '''Return the coordinates of the residues/atoms in the contact strings'''
      structure = parseAtomStruct(recFile)
      coords = []
      for model in structure:
        if areResidues:
          for res in model.get_residues():
            resName = f'{res.resname}{res.id[1]}'
            if resName in conStrs:
              coords += [atom.get_coord() for atom in res]
        else:
          for atom in model.get_atoms():
            atomName = f'{removeNumberFromStr(atom.name)}{atom.serial_number}'
            if atomName in conStrs:
              coords += [atom.get_coord()]
      return coords

    def getInputStructure(self):
        inFile = self.getReceptorPDB()
        structure = parseAtomStruct(inFile)
        return structure

    def getInputSurfaceStructure(self):
      structure = self.getInputStructure()
      structSurface = get_surface(structure[0], MSMS=Plugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))
      return structSurface

    def mapSurfaceCoords(self, oCoords, surfStruct):
        sCoords = []
        for coord in oCoords:
            closerSCoords = self.closerSurfaceCoords(coord, surfStruct)
            for cCoord in closerSCoords:
                cCoord = list(cCoord)
                if cCoord not in sCoords:
                  sCoords.append(cCoord)

        return sCoords

    def closerSurfaceCoords(self, coord, surfStruct):
        distances = distance.cdist([coord], surfStruct)
        closestIndexes = distances < self.maxDepth.get()
        return list(surfStruct[closestIndexes[0]])




