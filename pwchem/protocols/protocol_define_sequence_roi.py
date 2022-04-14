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
This protocol is used to import a set of pockets (of fpocket, p2rank, autoligand) from some files

"""
import os, json
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, min_dist, residue_depth
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *
from pwchem import Plugin

class ProtDefineSeqROI(EMProtocol):
    """
    Defines a set of pockets from a set of coordinates / residues / predocked ligands
    """
    _label = 'Define sequence ROI'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                      allowsNull=False, label="Input sequence: ",
                      help='Select the sequence object where the ROI will be defined')

        group.addParam('resPosition', params.StringParam,
                      label='Residues of interest: ',
                      help='Specify the residue to define a region of interest.\n'
                           'You can either select a single residue or a range '
                           '(it will take into account the first and last residues selected)')
        group.addParam('descrip', params.StringParam,
                       label='ROI description: ',
                       help='Specify some description for this region of interest')

        group.addParam('addResidue', params.LabelParam,
                      label='Add defined residue',
                      help='Here you can define a residue or range which will be added to the list of residues below.')
        group.addParam('inResidues', params.TextParam, width=70, default='',
                      label='Input residues: ',
                      help='Input residues to define the ROI.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        inpSeq = self.inputSequence.get()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))

        residuesStr = self.inResidues.get().strip().split('\n')
        for rStr in residuesStr:
            resDic = json.loads(rStr)
            roi, resIdxs = resDic['residues'], resDic['index']
            idxs = [int(resIdxs.split('-')[0]), int(resIdxs.split('-')[1])]

            roiSeq = Sequence(sequence=roi, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                              description=resDic['desc'])
            seqROI = SequenceROI(sequence=inpSeq, seqROI=roiSeq, roiIdx=idxs[0], roiIdx2=idxs[1])
            outROIs.append(seqROI)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)


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


