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
This protocol is used to define a set of sequences from AtomStruct, Sequence or SetOfSequences

"""

import os, json
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import SetOfSequences, Sequence

from pwchem.utils.utilsFasta import parseFasta


class ProtDefineSetOfSequences(EMProtocol):
    """
    AI Generated:

    Protocol: Define Set of Sequences

    This protocol generates a SetOfSequences object from multiple possible inputs:
    - Individual Sequence objects
    - AtomStruct objects (extracting chains)
    - External PDB codes

    It is designed to unify different sequence sources into a single
    Scipion-compatible SetOfSequences container for downstream analysis.

    Overview
    --------
    The protocol collects user-defined inputs and constructs a curated set of
    biological sequences, optionally specifying chain selection and positional
    metadata.

    Input options
    -------------
    inputOrigin:
        Defines the source type of the sequence:
        - Sequence: direct Sequence object
        - AtomStruct: extract sequence from a structure
        - PDB code: retrieve sequence from a PDB identifier

    inputSequence:
        Sequence object to be included (if inputOrigin = Sequence)

    inputAtomStruct:
        Structural object from which a chain sequence is extracted

    inputPDB:
        PDB identifier used to retrieve structural sequence

    inpChain:
        Chain identifier used when extracting sequence from structures or PDB

    inpPositions:
        Optional positional information associated with the sequence

    inputList:
        Internal list of encoded sequence entries to be processed

    Workflow
    --------
    1. User adds sequences from different origins
    2. Each input is stored in a structured JSON-like entry
    3. During execution, each entry is parsed
    4. FASTA files are read to reconstruct sequences
    5. Sequence objects are created and added to a SetOfSequences

    Output
    ------
    outputSequences:
        SetOfSequences object containing:
        - Sequence objects created from all selected inputs
        - Standardized identifiers (name/id)
        - Raw sequence strings extracted from FASTA representations

    Key Features
    ------------
    - Supports multiple input sources (Sequence, structure, PDB)
    - Converts heterogeneous inputs into a unified sequence set
    - Uses FASTA parsing for sequence reconstruction
    - Fully compatible with Scipion sequence workflows
    - Enables downstream comparative or alignment analysis

    Notes
    -----
    - Input list entries are expected in JSON format encoded in text lines
    - Sequence names are used as keys to extract FASTA entries
    - Only valid parsed entries are included in the final set
    """
    _label = 'Define set of sequences'

    # -------------------------- DEFINE param functions ----------------------
    def _addInputForm(self, form):
        form.addParam('inputSequence', params.PointerParam,
                      pointerClass='Sequence', allowsNull=True,
                      label="Input sequence: ", condition='inputOrigin==0',
                      help='Select the sequence object to add to the set')

        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input structure: ", condition='inputOrigin==1',
                      help='Select the AtomStruct object whose sequence to add to the set')

        form.addParam('inputPDB', params.StringParam, condition='inputOrigin==2',
                      label='Input PDB code: ', help='Specify the PDB code of the sequence to add')

        form.addParam('inpChain', params.StringParam,
                      label='Input chain: ', condition='inputOrigin in [1, 2]',
                      help='Specify the protein chain to use as sequence.')

        form.addParam('inpPositions', params.StringParam,
                      label='Input positions: ',
                      help='Specify the positions of the sequence to add in the output.')

        form.addParam('addInput', params.LabelParam,
                      label='Add input: ',
                      help='Add sequence to the output set')

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputOrigin', params.EnumParam, default=0,
                       label='Input origin: ', choices=['Sequence', 'AtomStruct', 'PDB code'],
                       help='Input sequence origin to add to the set')
        self._addInputForm(group)
        group.addParam('inputList', params.TextParam, width=100,
                      default='', label='List of inputs: ',
                      help='The list of input to use for the final output set.')

        form.addSection(label='Input Pointers')
        form.addParam('inputPointerLabels', params.LabelParam, important=True,
                      label='Records of inputs. Do not modificate manually',
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')
        form.addParam('inputPointers', params.MultiPointerParam, pointerClass="Sequence, AtomStruct",
                      label='Input Pointers: ', allowsNull=True,
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
      sequenceSet = SetOfSequences().create(outputPath=self._getPath())
      for inputLine in self.inputList.get().split('\n'):
          if inputLine.strip():
              inpJson = json.loads(inputLine.split(')')[1].strip())

              seqName = inpJson['name']
              seqDic = parseFasta(os.path.abspath(inpJson['seqFile']))
              seqObj = Sequence(sequence=seqDic[seqName], name=seqName, id=seqName)
              sequenceSet.append(seqObj)

      self._defineOutputs(outputSequences=sequenceSet)

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


