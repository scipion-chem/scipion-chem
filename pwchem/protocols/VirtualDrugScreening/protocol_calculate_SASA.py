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
from pwem.convert.atom_struct import toCIF, AtomicStructHandler, addScipionAttribute
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct

from pwchem.objects import SequenceChem
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

        group = form.addGroup('Sequence accessibility regions')
        group.addParam('extractSequence', params.BooleanParam, default=False,
                      label='Extract sequence with accesibility values: ',
                      help='Whether to output sequence with the accesibility values stored')
        group.addParam('chain_name', params.StringParam, condition='extractSequence',
                       label='Structure chain:', help="Chain of the structure to extract ")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('calculateSASAStep')
        self._insertFunctionStep('defineOutputStep')

    def calculateSASAStep(self, inProt=True):
        outSASA = self.getSASAFile(inProt)
        if not os.path.exists(outSASA):
            inputAS = self.inputAtomStruct.get().getFileName()

            if inProt:
                inFile = os.path.abspath(self._getExtraPath('inpdb.pdb'))
            else:
                inFile = os.path.abspath(self.getProject().getTmpPath('inpdb.pdb'))
            # Removes HETATM (sometimes yield error in Biopython)
            cleanPDB(inputAS, inFile, waters=True, hetatm=True)

            calculate_SASA(inFile, outSASA)

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

        if self.extractSequence.get():
            inChain = self.getChain()
            if inChain:
                inFile = inpStruct.getFileName()
                ASH.read(inFile)
                outStr = str(ASH.getSequenceFromChain(modelID=0, chainID=inChain))

                outSeq = SequenceChem(name='{}_{}'.format(getBaseName(inFile), inChain),
                                      sequence=outStr, id='{}_{}'.format(getBaseName(inFile), inChain),
                                      attributesFile=self._getExtraPath('sequenceAttributes.txt'))
                outSeq.addAttributes({'SASA': self.getSASAChainValues(inChain)})

                self._defineOutputs(outputSequence=outSeq)


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
    def getSASAFile(self, inProt=True):
        if inProt:
            inFile = self._getPath('sasa.txt')
        else:
            inFile = self.getProject().getTmpPath('sasa.txt')
        return os.path.abspath(inFile)

    def getSASAChainValues(self, chain):
        vals, sasaDic = [], self.getSASADic()
        for res, val in sasaDic.items():
            if res.split(':')[0] == chain:
                vals.append(val)
        return vals

    def getSASADic(self):
        sasaDic = {}
        with open(self.getSASAFile()) as fIn:
            for line in fIn:
                sasaDic[line.split()[0]] = float(line.split()[1])
        return sasaDic

    def getChain(self):
      chainJson = self.chain_name.get()
      if chainJson:
          chain = json.loads(chainJson)['chain']
      else:
          chain = ''
      return chain
