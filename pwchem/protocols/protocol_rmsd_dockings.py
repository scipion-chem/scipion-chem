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
This protocol is used to score docking positions obtained by several software using descriptors and score functions
available in the Open Drug Discovery Toolkit (ODDT, https://github.com/oddt/oddt)

"""
import os, re

from pyworkflow.utils import Message
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *

class ProtocolRMSDDocking(EMProtocol):
    """
    Executes the scoring of a set of molecules which have been previously docked.
    """
    _label = 'RMSD docking'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMols', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Docked Small Molecules: ",
                       help='Select the docked molecules to measure the RMSD against the reference AtomStruct')
        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct',
                      label="Input AtomStruct with ligand: ",
                      help='Select the AtomStruct with an attached ligand to measure the RMSD against the '
                           'input small molecules. ')

        form.addParam('onlyHeavy', params.BooleanParam, default=True,
                      label='Use only heavy atoms for RMSD: ', expertLevel=params.LEVEL_ADVANCED,
                      help='Use only non-hydrogen atoms for RMSD calculation')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('rmsdStep')
        self._insertFunctionStep('createOutputStep')

    def rmsdStep(self):
        rmsdDic = {}
        for mol in self.inputSmallMols.get():
            AS = self.inputAtomStruct.get()
            posDic1, posDic2 = self.getLigandPosDic(mol), self.getLigandPosDic(AS)
            k1, k2 = list(posDic1.keys()), list(posDic2.keys())
            k1.sort(), k2.sort()
            if self.checkSameKeys(posDic1, posDic2):
                print('Coords1: ', posDic1.values())
                print('Coords2: ', posDic2.values())
                rmsd = calculateRMSD(posDic1.values(), posDic2.values())
                print('RMSD: ', rmsd)
                rmsdDic[os.path.abspath(mol.getPoseFile())] = rmsd
                print('Mol: {}. RMSD: {}'.format(os.path.basename(mol.getPoseFile()), rmsd))
            else:
                rmsd = 10000

        with open(self._getExtraPath('rmsd.txt'), 'w') as f:
            for k in rmsdDic:
                f.write('{}\t{}\n'.format(k, rmsdDic[k]))

    def createOutputStep(self):
        rmsdDic = {}
        with open(self._getExtraPath('rmsd.txt')) as f:
            for line in f:
                rmsdDic[line.split()[0]] = line.split()[1]

        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMols.get(), self._getPath(), copyInfo=True)
        for mol in self.inputSmallMols.get():
            #Specific attribute name for each score?
            k = os.path.abspath(mol.getPoseFile())
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
        if not self.inputSmallMols.get().isDocked():
            validations += ['Input Small Molecules is not docked yet\n']

        return validations

    def _warnings(self):
        warnings = []
        return warnings

    ####################################### UTILS functions #############################

    def getLigandPosDic(self, item):
        onlyHeavy = self.onlyHeavy.get()
        if issubclass(type(item), SmallMolecule):
            posDic = item.getAtomsPosDic(onlyHeavy)
        elif issubclass(type(item), AtomStruct):
            posDic = {}
            if '.pdb' in item.getFileName():
                with open(item.getFileName()) as fIn:
                    for line in fIn:
                        if line.startswith('HETATM'):
                            elements = splitPDBLine(line)
                            atomId, atomType, coords = elements[2], elements[-1], elements[6:9]
                            if atomType != 'H' or not onlyHeavy:
                                posDic[atomId] = list(map(float, coords))
            # if cif
        return posDic

    def checkSameKeys(self, d1, d2):
        return set(d1.keys()) == set(d2.keys())
