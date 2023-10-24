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

"""
import os
from Bio.PDB import PDBParser, MMCIFParser

from pyworkflow.utils import Message
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *

AS, MOLS = 0, 1

class ProtocolRMSDDocking(EMProtocol):
    """
    Calculates the RMSD between a docked molecule and another molecule in the same receptor, either from an AtomStruct
    or a docked SetOfSmallMolecules
    """
    _label = 'RMSD docking'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMolecules', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', label="Input Docked Small Molecules: ",
                       help='Select the docked molecules to measure the RMSD against the reference AtomStruct')

        form.addParam('refOrigin', params.EnumParam, default=0,
                      label='Reference molecule origin: ', choices=['AtomStruct', 'SmallMolecule'], 
                      help='Where to extract the reference molecule from.')
        
        form.addParam('refAtomStruct', params.PointerParam, condition='refOrigin==0',
                      pointerClass='AtomStruct', label="Input AtomStruct with ligand: ",
                      help='Select the AtomStruct with an attached ligand to measure the RMSD against the '
                           'input small molecules. ')
        form.addParam('refLigName', params.StringParam, condition='refOrigin==0', label="Reference molecule: ",
                      help='Select the reference molecule to measure the RMSD against')

        form.addParam('refSmallMolecules', params.PointerParam, condition='refOrigin==1',
                      pointerClass='SetOfSmallMolecules', label="Input reference molecules: ",
                      help='Select the SetOfMolecules where the reference is included')
        form.addParam('refMolName', params.StringParam, condition='refOrigin==1', label="Reference molecule: ",
                      help='Select the reference molecule to measure the RMSD against')

        form.addParam('onlyHeavy', params.BooleanParam, default=True,
                      label='Use only heavy atoms for RMSD: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Use only non-hydrogen atoms for RMSD calculation')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('rmsdStep')
        self._insertFunctionStep('createOutputStep')

    def rmsdStep(self):
        rmsdDic = {}
        if self.refOrigin.get() == AS:
            refPosDic = self.getLigandPosDic(self.refAtomStruct.get(), self.refLigName.get())
        else:
            for mol in self.refSmallMolecules.get():
                if self.refMolName.get() == mol.__str__():
                    refPosDic = self.getLigandPosDic(mol)
                    break
        
        for mol in self.inputSmallMolecules.get():
            posDic = self.getLigandPosDic(mol)
            if hasattr(mol, '_mappingFile'):
                posDic = self.mapLabels(mol, posDic)

            k1, k2 = list(posDic.keys()), list(refPosDic.keys())
            k1.sort(), k2.sort()
            if self.checkSameKeys(posDic, refPosDic):
                rmsd = calculateRMSDKeys(posDic, refPosDic)
                rmsdDic[mol.getPoseFile()] = rmsd
                print('Mol: {}. RMSD: {}'.format(os.path.basename(mol.getPoseFile()), rmsd))
            else:
                rmsd = 10000
        with open(self._getExtraPath('rmsd.tsv'), 'w') as f:
            for k in rmsdDic:
                f.write('{}\t{}\n'.format(k, rmsdDic[k]))

    def createOutputStep(self):
        rmsdDic = {}
        with open(self._getExtraPath('rmsd.tsv')) as f:
            for line in f:
                rmsdDic[line.split()[0]] = line.split()[1]

        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        for mol in self.inputSmallMolecules.get():
            #Specific attribute name for each score?
            k = mol.getPoseFile()
            if k in rmsdDic:
                setattr(mol, "_rmsdToRef", Float(rmsdDic[k]))
            else:
                setattr(mol, "_rmsdToRef", Float(10000))
            newMols.append(mol)
        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        if not self.inputSmallMolecules.get().isDocked():
            validations += ['Input Small Molecules is not docked yet\n']

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    ####################################### UTILS functions #############################

    def mapLabels(self, mol, posDic):
        '''Map the atom labels which were normally reorganized during the ligand preparation to those in the original
        molecule'''
        mapDic, pDic = mol.getMapDic(), {}
        for prevLabel in posDic:
            pDic[mapDic[prevLabel.upper()]] = posDic[prevLabel.upper()]
        return pDic

    def getLigandPosDic(self, item, molName=None):
        onlyHeavy = self.onlyHeavy.get()
        if issubclass(type(item), SmallMolecule):
            posDic = item.getAtomsPosDic(onlyHeavy)

        elif issubclass(type(item), AtomStruct):
            posDic = {}
            ASFile = item.getFileName()
            if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
                pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
                parser = PDBParser().get_structure(pdb_code, ASFile)
            elif ASFile.endswith('.cif'):
                pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
                parser = MMCIFParser().get_structure(pdb_code, ASFile)
            else:
                print('Unknown AtomStruct file format')

            for model in parser:
                for chain in model:
                    for residue in chain:
                        if isHet(residue) and residue.resname == molName:
                            for atom in residue:
                                atomId, coords = atom.get_id(), atom.get_coord()
                                if not atomId.startswith('H') or not onlyHeavy:
                                    posDic[atomId] = list(map(float, coords))

        return posDic

    def checkSameKeys(self, d1, d2):
        return set(d1.keys()) == set(d2.keys())
