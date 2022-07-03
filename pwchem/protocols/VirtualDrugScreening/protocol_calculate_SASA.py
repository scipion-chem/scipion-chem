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
This protocol is used to calculate the SASA of an AtomStruct

"""
import os, json

from pyworkflow.utils import Message
from pyworkflow.protocol import params
from pwem.convert.atom_struct import toPdb, toCIF, AtomicStructHandler, addScipionAttribute
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem import Plugin as pwchemPlugin
from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from pwchem.constants import MGL_DIC

class ProtCalculateSASA(EMProtocol):
    """
    Calculate SASA of an AtomStruct
    """
    _label = 'Calculate SASA'
    _ATTRNAME = 'SASA'
    _OUTNAME = 'outputAtomStruct'
    _possibleOutputs = {_OUTNAME: AtomStruct}

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      allowsNull=False, label="Input AtomStruct: ",
                      help='Select the AtomStruct object where the structural ROIs will be defined')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculateSASAStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateSASAStep(self):
        inputAS = self.inputAtomStruct.get().getFileName()
        outSASA = self.getSASAFile()
        args = '{} {}'.format(os.path.abspath(inputAS), outSASA)
        pwchemPlugin.runScript(self, 'calculate_SASA.py', args, env='plip', cwd=self._getExtraPath())

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        outStructFileName = self._getPath('outputStructureSASA.cif')
        # Write conservation in a section of the output cif file
        ASH = AtomicStructHandler()
        sasaDic = self.getSASADic()
        inpAS = toCIF(self.inputAtomStruct.get().getFileName(), self._getTmpPath('inputStruct.cif'))
        cifDic = ASH.readLowLevel(inpAS)
        cifDic = addScipionAttribute(cifDic, sasaDic, self._ATTRNAME)
        ASH._writeLowLevel(outStructFileName, cifDic)

        AS = AtomStruct(filename=outStructFileName)
        self._defineOutputs(outputAtomStruct=AS)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getSASAFile(self):
        return os.path.abspath(self._getPath('sasa.txt'))

    def getSASADic(self):
        sasaDic = {}
        with open(self.getSASAFile()) as fIn:
            for line in fIn:
                sasaDic[line.split()[0]] = float(line.split()[1])
        return sasaDic
