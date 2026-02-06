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

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem.protocols import ProtocolConsensusDocking

AS, MOLS = 0, 1


class ResidueNameSelect(Select):
    def __init__(self, resname):
        self.resname = resname

    def accept_residue(self, residue):
        return residue.get_resname() == self.resname

class ProtocolRMSDDocking(ProtocolConsensusDocking):
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

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        cStep = self._insertFunctionStep(self.convertInputStep, prerequisites=[])
        self._insertFunctionStep(self.createOutputStep, prerequisites=[cStep])

    def createOutputStep(self):
        if self.refOrigin.get() == AS:
            refFile = self.writeRefMol(self.refAtomStruct.get().getFileName(), self.refLigName.get())
        else:
            for mol in self.refSmallMolecules.get():
                if self.refMolName.get() == mol.__str__():
                    refFile = mol.getPoseFile() if mol.getPoseFile() is not None else mol.getFileName()
                    refFile = os.path.abspath(refFile)
                    break

        rmsdDic = runRdkitRMSD(self, self.getAllInputMols(), referenceFile=refFile)


        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        for mol in self.inputSmallMolecules.get():
            molName = mol.getMolName()
            if molName in rmsdDic and len(rmsdDic[molName]) > 0:
                newRMSD = rmsdDic[molName].pop(0)
                setattr(mol, "_rmsdToRef", Float(newRMSD))
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

    def writeRefMol(self, ASFile, molName):
        if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
            pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
            structure = PDBParser().get_structure(pdb_code, ASFile)
        elif ASFile.endswith('.cif'):
            pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
            structure = MMCIFParser().get_structure(pdb_code, ASFile)
        else:
            print('Unknown AtomStruct file format')

        oFile = os.path.abspath(self._getExtraPath(f"{molName}.pdb"))
        io = PDBIO()
        io.set_structure(structure)
        io.save(oFile, ResidueNameSelect(molName))

        return oFile

    def checkSameKeys(self, d1, d2):
        return set(d1.keys()) == set(d2.keys())

    def getAllInputMols(self):
        mols = []
        convMolNames = self.getConvMolNames()
        for mol in self.inputSmallMolecules.get():
            newMol = mol.clone()
            molFileName = getBaseName(mol.getPoseFile())
            if molFileName in convMolNames:
                newMol.setPoseFile(os.path.relpath(convMolNames[molFileName]))
            mols.append(newMol)
        return mols

    def getOriginalReceptorFile(self):
        return self.inputSmallMolecules.get().getProteinFile()
