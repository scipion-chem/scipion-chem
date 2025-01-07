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
This protocol maps an attribute stored on AtomStruct or SequenceChem residues to the residues of a similar sequence
and calculates and stores the mean value of the attribute for each of its ROIs as attributes.

"""
import os, json
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.convert import cifToPdb
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import NAME, SPEC, VALUE
from pwem.viewers.viewer_localres import getStructureRecipient, getResiduePositions, getResidueAverage
from pwem.convert.atom_struct import AtomicStructHandler

from pwchem.utils import Float, runOpenBabel
from pwchem.utils.utilsFasta import pairwiseAlign

STRUCTURE, SEQUENCE = 0, 1
fastaName = "pairWise.fasta"

class ProtMapAttributeToSeqROIs(EMProtocol):
    """
    Maps a set of SequenceROIs to their respective structure ROIs in an AtomStruct
    """
    _label = 'Map attribute to sequence ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input sequence')
        group.addParam('inputSequenceROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       allowsNull=False, label="Input sequence ROIs: ",
                       help='Select the Input Sequence ROIs that will be labeled')

        group = form.addGroup('Map from')
        group.addParam('inputFrom', params.EnumParam, default=STRUCTURE,
                       label='Input from: ', choices=['AtomStruct', 'Sequence'],
                       help='Type of attributes input you want to use')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      label="Input AtomStruct: ", condition='inputFrom==0',
                      help='Select the AtomStruct object where the pockets will be defined')

        group.addParam('chain_name', params.StringParam, condition='inputFrom==0',
                      label='Chain of interest: ',
                      help='Specify the chain of the residue of interest')

        group.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       label="Input sequence: ", condition='inputFrom==1',
                       help='Input sequence where attributes are stored')

        group.addParam('attrName', params.StringParam, label='Attribute to map: ',
                       help='Specify the attribute of the input to map to the ROIs')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('alignSequencesStep')
        self._insertFunctionStep('defineOutputStep')

        self._mapWarning = 'Mapping of ROIs not posible, check the alignment in ', self._getPath(fastaName)


    def alignSequencesStep(self):
        seq, seqName = self.getInputSequence()
        seqAttr, seqNameAttr = self.getInputAttrSequence()

        # Alignment
        outFile = os.path.abspath(self._getPath(fastaName))
        pairwiseAlign(seq, seqAttr, outFile, seqName1=seqName, seqName2=seqNameAttr)

    def defineOutputStep(self):
        attrVals = self.getAttrVals()
        mapDic = self.mapResidues()   #{idxROISeq: idxAttrSeq}

        inpSet = self.inputSequenceROIs.get()
        outSet = inpSet.createCopy(self._getPath(), copyInfo=True)
        for roi in inpSet:
            attrCurVals = []
            roiIdxs = roi.getROIIdxs()
            for idx in range(roiIdxs[0], roiIdxs[1] + 1):
                if idx in mapDic:
                    attrCurVals.append(attrVals[mapDic[idx]])

            if attrCurVals:
                attrAv = sum(attrCurVals)/len(attrCurVals)
            else:
                attrAv = 0
            setattr(roi, self.attrName.get(), Float(attrAv))
            outSet.append(roi)

        self._defineOutputs(outputROIs=outSet)





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
    def _getStructureAttributes(self):
        '''Returns a list with the names of the attributes of the AtomStruct object'''
        obj = self.inputAtomStruct.get()
        ASH = AtomicStructHandler()
        cifDic = ASH.readLowLevel(obj.getFileName())
        return list(set(cifDic['_scipion_attributes.name']))

    def getAttrVals(self):
        '''Returnn the values of the attribute'''
        attrName = self.attrName.get()
        if self.inputFrom.get() == STRUCTURE:
            cifFile = self.inputAtomStruct.get().getFileName()
            cifDic = AtomicStructHandler().readLowLevel(cifFile)

            names, values, specs = np.array(cifDic[NAME]), np.array(cifDic[VALUE]), np.array(cifDic[SPEC])
            recipient = getStructureRecipient(cifDic, attrName)
            attrValues, attrSpecs = values[names == attrName], specs[names == attrName]
            if recipient == 'atoms':
                attrValues, attrSpecs = getResidueAverage(attrValues, attrSpecs)

            _, attrChainValues = getResiduePositions(attrSpecs, attrValues, self.getChainName())

            attrValues = list(map(float, attrChainValues))

        else:
            attrValues = self.inputSequence.get().getAttributesDic()[attrName]
        return attrValues

    def getInputAttributes(self):
        '''Return the names of the attributes of the input sequence to map on the ROIs
        '''
        if self.inputFrom.get() == STRUCTURE:
            attrs = self._getStructureAttributes()
        else:
            attrs = list(self.inputSequence.get().getAttributesDic().keys())
        return attrs

    def getASFileName(self):
        inpStruct = self.inputAtomStruct.get()
        inpFile = inpStruct.getFileName()
        basename = os.path.basename(inpFile).split('.')[0]
        if inpFile.endswith('.cif'):
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            cifToPdb(inpFile, inpPDBFile)

        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            inpStruct.convert2PDB(outPDB=inpPDBFile)

        elif inpFile.endswith('.pdbqt'):
            inpPDBFile = os.path.abspath(self._getExtraPath(basename + '.pdb'))
            args = ' -ipdbqt {} -opdb -O {}'.format(os.path.abspath(inpFile), inpPDBFile)
            runOpenBabel(protocol=self, args=args, cwd=self._getTmpPath())

        else:
            inpPDBFile = inpFile

        return inpPDBFile

    def getASName(self):
      return os.path.splitext(os.path.basename(self.getASFileName()))[0]

    def getInputSequence(self):
        inputObj = getattr(self, 'inputSequenceROIs').get()
        seq = inputObj.getSequence()
        seqName = inputObj.getSequenceObj().getId()
        if not seqName:
            seqName = inputObj.getSeqName()
        return seq, seqName

    def getInputAttrSequence(self):
        if self.inputFrom.get() == STRUCTURE:
            from pwem.convert.atom_struct import AtomicStructHandler
            asFile = self.getASFileName()
            seqName = os.path.basename(asFile)
            handler = AtomicStructHandler(asFile)
            chainName = getattr(self, 'chain_name').get()

            # Parse chainName for model and chain selection
            struct = json.loads(chainName)  # From wizard dictionary
            chainID, modelId = struct["chain"].upper().strip(), int(struct["model"])

            seq = str(handler.getSequenceFromChain(modelID=modelId, chainID=chainID))
        else:
            seqObj = self.inputSequence.get()
            seq, seqName= seqObj.getSequence(), seqObj.getSeqName()
        return seq, seqName

    def mapResidues(self):
        '''Returns a dictionary which maps the idxs of the residues of the ROIs sequence and
        the sequence from the attributes input
        {idxROISeq: idxAttrSeq}
        '''
        seq, seqAttr = self.parseAlignment()

        mapDic = {}
        i, j = 0, 0
        for k in range(len(seq)):
            if seq[k] == '-':
                i = i + 1
            if seqAttr[k] == '-':
                j = j + 1
            if seq[k] != '-' and seqAttr[k] != '-':
                mapDic[k-i+1] = k-j

        return mapDic

    def getChainName(self):
        return json.loads(self.chain_name.get())['chain']

    def parseAlignment(self):
        seq, seqAS = '', ''
        first = True
        with open(self._getPath(fastaName)) as f:
            f.readline()
            for line in f:
                if not line.startswith('>'):
                    if first:
                        seq += line.strip()
                    else:
                        seqAS += line.strip()
                else:
                    first = False
        return seq, seqAS







