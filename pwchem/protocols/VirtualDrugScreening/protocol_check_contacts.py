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
This protocol is to check the residues contacts from a set of small molecules docked into a receptor.
The user can specify a threshold, and all residues with atoms below that threshold are stored.
"""
import os, glob, time
from scipy import spatial

from pyworkflow.protocol import params
from pyworkflow.utils import Message, createLink
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem.constants import RDKIT_DIC
from pwchem import Plugin

RESIDUE, ATOM = 0, 1

class ProtocolContactsDocking(EMProtocol):
    """
    """
    _label = 'Docking contacts'

    def __init__(self, **args):
        super().__init__(**args)
        self.stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMolecules', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Docked Small Molecules: ",
                       help='Select the docked molecules whose contacts you want to get')

        form.addParam('threshold', params.FloatParam, default=5.0, label='Distance threshold (A): ',
                       help="Distance threshold where atoms of ligand-residue will be considered in contact")
        form.addParam('contactLevelRec', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST, expertLevel=params.LEVEL_ADVANCED,
                      default=0, label='Store receptor contacts at level: ', choices=['Residue', 'Atom'],
                      help="Whether to store the residue or the atom in contact for the receptor")
        form.addParam('contactLevelLig', params.EnumParam,
                      display=params.EnumParam.DISPLAY_HLIST, expertLevel=params.LEVEL_ADVANCED,
                      default=0, label='Store ligand contacts at level: ', choices=['Residue', 'Atom'],
                      help="Whether to store the residue or the atom in contact for the ligand")

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        cStep = self._insertFunctionStep('getContactsStep', prerequisites=[])
        self._insertFunctionStep('createOutputStep', prerequisites=[cStep])


    def performContactAnalysis(self, molFns, molLists, it, recFile, outDir):
        for ligFile in molFns:
            recAtoms, ligAtoms = self.parseRecLigAtoms(ASFile=recFile, ligFile=ligFile)

            contactDic = self.getLigandContacts(recAtoms, ligAtoms)
            # contactList = self.getAllContactPairs(contactDic)
            self.writeContactFile(contactDic, ligFile, outDir)


    def getContactsStep(self):
        recFile = self.inputMolecules.get().getProteinFile()
        molFns = [mol.getPoseFile() for mol in self.inputMolecules.get()]

        outDir = self.getOutputDir()
        os.mkdir(outDir)
        performBatchThreading(self.performContactAnalysis, molFns, self.numberOfThreads.get(), cloneItem=False,
                              recFile=recFile, outDir=outDir)


    def parseOutFiles(self):
        '''Return a dictionary of form: {ligandFile: {ligAtom: [recAtoms]}}
        with the contacts of atoms or residues stored in each of the output files.'''
        fileDic = {}
        outDir = self.getOutputDir()
        for file in os.listdir(outDir):
            with open(os.path.join(outDir, file)) as f:
                for line in f:
                    if line.startswith('LIGFILE'):
                        ligFile = line.split('::')[1].strip()
                        fileDic[ligFile] = {}
                    elif line.startswith('LIGAND'):
                        ligName = line.split('::')[1].strip()
                    elif line.startswith('\tRECEPTOR'):
                        recNames = line.split('::')[1].strip().split(',')
                        fileDic[ligFile][ligName] = recNames
        return fileDic


    def createOutputStep(self):
        inpSet = self.inputMolecules.get()
        outSet = inpSet.createCopy(self._getPath(), copyInfo=True)

        contactDic = self.parseOutFiles()
        for mol in inpSet:
            molFile = os.path.abspath(mol.getPoseFile())
            if molFile in contactDic:
                cDic = contactDic[molFile]
                contactFile = self.getContactFileName(molFile, self.getOutputDir())
                setattr(mol, "contactsFile", String(os.path.relpath(contactFile)))

                cStrs = []
                for ligAtom, recAtoms in cDic.items():
                    cStrs.append(f'{ligAtom}_{",".join(recAtoms)}')
                setattr(mol, "recContacts", String('-'.join(cStrs)))
            else:
                print(f'No contacts found for {molFile}')
                setattr(mol, "contactsFile", String(None))
                setattr(mol, "recContacts", String(None))

            outSet.append(mol)
        self._defineOutputs(outputSmallMolecules=outSet)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        pSet = self.inputMolecules.get()
        if not pSet.isDocked():
            validations.append('Sets of input molecules must be docked first.\n'
                               'Set: {} has not been docked'.format(pSet))

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    def getOutputDir(self):
        return os.path.abspath(self._getExtraPath('contacts'))

    # --------------------------- UTILS functions -----------------------------------

    def parseSingleFile(self, ASFile):
        '''Parse all the atoms stored in a AtomStructFile'''
        atoms = []
        parser = parseAtomStruct(ASFile)

        for model in parser:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atoms.append(atom)
        return atoms

    def parseCombinedFile(self, ASFile, ligName=None):
        '''Parse a combined AtomStruct file trying to separate the receptor from the ligand atoms'''
        recAtoms, ligAtoms = [], []
        parser = parseAtomStruct(ASFile)

        for model in parser:
            for chain in model:
                for residue in chain:
                    if is_het(residue) and (residue.resname == ligName or not ligName):
                        for atom in residue:
                            ligAtoms.append(atom)
                    else:
                        for atom in residue:
                            recAtoms.append(atom)
        return recAtoms, ligAtoms

    def parseRecLigAtoms(self, ASFile, ligFile=None, ligName=None):
        '''Return two list of the form with the atoms for the receptor and ligand or heteroatoms
        '''
        if ligFile:
            recAtoms, ligAtoms = self.parseSingleFile(ASFile), self.parseSingleFile(ligFile)
        else:
            recAtoms, ligAtoms = self.parseCombinedFile(ASFile, ligName=ligName)
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

        dists = spatial.distance.cdist(list(recCoords.values()), list(ligCoords.values()))
        idxs = np.argwhere(dists <= self.threshold.get())
        for idx in idxs:
            recAtom, ligAtom = list(recCoords.keys())[idx[0]], list(ligCoords.keys())[idx[1]]
            if ligAtom in contactAtoms:
                contactAtoms[ligAtom].append(recAtom)
            else:
                contactAtoms[ligAtom] = [recAtom]
        return contactAtoms

    def getNameFromAtom(self, atom, contactLevel=RESIDUE):
        if contactLevel == ATOM:
            name = f'{removeNumberFromStr(atom.name)}{atom.serial_number}'
        else:
            name = f'{atom.get_parent().resname}{atom.get_parent().id[1]}'
        return name

    def buildContact(self, ligAtom, recAtom):
        ligStr = self.getNameFromAtom(ligAtom, self.contactLevelLig.get())
        recStr = self.getNameFromAtom(recAtom, self.contactLevelRec.get())
        return f'{recStr}_{ligStr}'

    def getAllContactPairs(self, cDic):
        '''Build the String object that encodes the contacts'''
        contacts = []
        for ligAtom in cDic:
            for recAtom in cDic[ligAtom]:
                contact = self.buildContact(ligAtom, recAtom)
                if not contact in contacts:
                    contacts.append(contact)
        return contacts

    def getContactFileName(self, ligFile, outDir):
        return os.path.join(outDir, getBaseName(ligFile) + '_contacts.txt')

    def writeContactFile(self, contactDic, ligFile, outDir):
        ligLvlStr = 'RESIDUE' if self.contactLevelLig.get() == RESIDUE else 'ATOM'
        recLvlStr = 'RESIDUE' if self.contactLevelRec.get() == RESIDUE else 'ATOM'

        outFile = self.getContactFileName(ligFile, outDir)
        with open(outFile, 'w') as fOut:
            fOut.write(f'## RESIDUE representation: ResidueName.ResidueNumber(Number of residue in sequence)\n')
            fOut.write(f'## ATOM representation: AtomName.AtomSerialNumber(Number of atom in whole structure)\n')
            fOut.write(f'LIGFILE:: {os.path.abspath(ligFile)}\n')

            ligDic = {self.getNameFromAtom(ligAtom, self.contactLevelLig.get()): ligAtom for ligAtom in contactDic}
            ligStrs = natural_sort(list(ligDic.keys()))
            for ligStr in ligStrs:
                recStrs = [self.getNameFromAtom(recAtom, self.contactLevelRec.get())
                           for recAtom in contactDic[ligDic[ligStr]]]
                recStrs = natural_sort(set(recStrs))

                fOut.write(f'LIGAND {ligLvlStr} :: {ligStr}\n')
                fOut.write(f'\tRECEPTOR {recLvlStr} :: {",".join(recStrs)}\n')






