# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This protocol is used to obtain peptides from sequences.

"""
import glob

import requests
from pwem.objects import AtomStruct, Sequence, SetOfSequences, SetOfAtomStructs

import pwchem.constants as constants
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwchem import Plugin

URL="https://api.esmatlas.com/foldSequence/v1/pdb/"

class ProtPeptideFromSequence(EMProtocol):
    """
    Generates a pdb peptide file from a sequence with ESM Metagenomics Atlas (https://github.com/facebookresearch/esm).
    """
    _label = 'Get peptide from sequence'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('seqObject', params.BooleanParam, default=False,
                       label='Input sequence objects: ',
                       help='Choose whether to input sequence objects or string.')

        group.addParam('name', params.StringParam, default=' ',
                       label='Input sequences names: ', condition='not seqObject',allowsNull=True,
                       help='Input sequence name. Separate names with ",". Eg: NAME1, NAME2')
        group.addParam('inputSeq', params.StringParam, default='',
                       label='Input sequences: ', condition='not seqObject', allowsNull=True,
                       help='Input sequence to get peptide pdb. Separate sequences with ",". Eg: IKILAVR, paffaktsav')

        group.addParam('set', params.BooleanParam, default=False,
                       label='Input a set of sequences: ', condition='seqObject',
                       help='Choose whether to input a set of sequences or individual objects.')
        group.addParam('inputSeqObject', params.PointerParam, pointerClass=Sequence,
                       label='Input sequence: ', condition='seqObject and not set', allowsNull=True,
                       help='Input sequence to get peptide pdb.')
        group.addParam('inputSeqSet', params.PointerParam, pointerClass=SetOfSequences,
                       label='Input sequences: ', condition='seqObject and set', allowsNull=True,
                       help='Input sequences to get peptide pdb.')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('getESMpdbAPIStep')
        self._insertFunctionStep('createOutputStep')

    def getESMpdbAPIStep(self):
        if (not self.seqObject.get()):
            sequences = self.inputSeq.get()
            sequenceList = [s.strip().upper() for s in sequences.split(",") if s.strip()]

            names = self.name.get()
            namesList = [s.strip().upper() for s in names.split(",") if s.strip()]

            if len(sequenceList) != len(namesList):
                raise ValueError("Number of sequences and names must match.")

            for seq, name in zip(sequenceList, namesList):
                self.doRequest(seq,name)

        elif (self.seqObject.get() and self.set.get()):
            for seq in self.inputSeqSet.get():
                sequence = seq.getSequence()
                name = seq.getSeqName()
                self.doRequest(sequence,name)

        elif (self.seqObject.get() and not self.set.get()):
            seq = self.inputSeqObject.get().getSequence()
            name = self.inputSeqObject.get().getSeqName()
            self.doRequest(seq, name)

    def createOutputStep(self):
        pdbFilesDir = self.getPath()
        pdbFiles = glob.glob(os.path.join(pdbFilesDir, "*.pdb"))

        if not pdbFiles:
            self.warning("No PDB files found in:", pdbFilesDir)
            return

        if len(pdbFiles)>1:
            setOfAtomStructs = SetOfAtomStructs().create(self._getPath())

            for pdbFile in pdbFiles:
                atomStruct = AtomStruct(filename=pdbFile)
                setOfAtomStructs.append(atomStruct)

            self._defineOutputs(outputAtomStructs=setOfAtomStructs)
        else:
            pdbFile = pdbFiles[0]
            atomStruct = AtomStruct(filename=pdbFile)
            self._defineOutputs(outputAtomStruct=atomStruct)



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

    # --------------------------- UTILS functions ------------------------
    def doRequest(self, seq, name):
        response = requests.post(URL, data=seq)
        pdbFile = self.getPath(f'{name}.pdb')

        if response.status_code == 200:
            with open(pdbFile, "w") as f:
                f.write(response.text)
        else:
            self.error(f"Error: {response.status_code} {response.text}")