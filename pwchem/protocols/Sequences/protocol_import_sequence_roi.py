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
This protocol is used to import a set of sequence ROIs (regions of interest) from some files

"""
import csv

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import SequenceHandler

from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *

class ProtImportSeqROI(EMProtocol):
    """
    Defines a set of sequence ROIs from a file containing those ROIs
    """
    _label = 'Import sequence ROI'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputFile', params.PathParam,
                       label='Input ROIs file: ', allowsNull=False,
                       help='Input file where the ROIs are stored in a determined format')

        group.addParam('fileFormat', params.EnumParam, default=1,
                      choices=['None of listed', 'IEDB'],
                      label='Input file format: ',
                      help="Format of the input file")
        #todo: if not IEDB, try to define by separator and columns of interest
        # group.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
        #               allowsNull=True, label="Input sequence: ",
        #               help='Select the sequence object where the ROI will be defined')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        if self.fileFormat.get() == 0:
            pass
        elif self.fileFormat.get() == 1:
            roiDic = self.parseIEDB()

        for uniprotID in roiDic:
            outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs_{}.sqlite'.format(uniprotID)))
            for roiObj in roiDic[uniprotID]['rois']:
                outROIs.append(roiObj)

            if len(outROIs) > 0:
                outDic = {'outputROIs_{}'.format(uniprotID): outROIs}
                self._defineOutputs(**outDic)


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

    def parseIEDB(self):
        '''Parses the ROIs epitopes file from IEDB and returns a dictionary with UniprotID of the
        parent sequence as keys for dictionaries with the ROI objets:
        {uniprotID1: {"sequence": parentSequence, "rois": [roiSequence1, roiSequence2,...]},
         uniprotID2: {...}
        }
        '''
        roiDic = {}
        with open(self.inputFile.get()) as f:
            f.readline()
            f.readline()
            csvData = list(csv.reader(f, quotechar='"', delimiter=',', quoting=csv.QUOTE_ALL,
                                      skipinitialspace=True))
            for sline in csvData:
                if 'Discontinuous peptide' in sline[1]:
                    roiStrs, roiIdxsPairs = self.parseDiscontinuous(sline)
                else:
                    roiStrs, roiIdxsPairs = [sline[2]], [sline[5:7]]

                for roiStr, roiIdxs in zip(roiStrs, roiIdxsPairs):
                    roiDic = self.addROI2Dic(roiDic, sline, roiStr, roiIdxs)

        return roiDic

    def parseDiscontinuous(self, sline):
        roiStr, roiID, roiIdxs, roiDesc = sline[2], sline[0], sline[5:7], sline[10]
        rois, roiIdxs, prevIdx = [], [], ''
        if len(roiStr.split(':')) > 1:
            return [], []

        for res in roiStr.split(','):
            res = res.strip()
            if res:
                aa, idx = res[0], int(res[1:])
                if prevIdx == '':
                    #Start ROI parsing
                    rois.append(aa)
                    roiIdxs.append([idx, ''])
                    prevIdx = idx
                else:
                    step = idx - prevIdx
                    if step > 6:
                        #ROI ends on previous
                        roiIdxs[-1][1] = prevIdx
                        rois.append(aa)
                        roiIdxs.append([idx, ''])
                    else:
                        #Grow ROI
                        rois[-1] += '-'*(step-1) + aa
                    prevIdx = idx

        roiIdxs[-1][1] = idx
        return rois, roiIdxs

    def addROI2Dic(self, roiDic, sline, roiStr, roiIdxs):
        roiID, roiDesc, seqDesc, uniprotID = sline[0], sline[10], sline[11], sline[12]
        if uniprotID.strip():
            if not uniprotID in roiDic:
                print('Looking for {} in UniProt'.format(uniprotID))
                seqDic, error = SequenceHandler().downloadSeqFromDatabase(uniprotID, dataBase='UniProt')
                if seqDic is None:
                    print("Error: ", error)
                else:
                    sequence = seqDic['sequence']
                    seq = Sequence(name=uniprotID, id=uniprotID, sequence=sequence,
                                   description=sline[11])
                    roiDic[uniprotID] = {'sequence': seq, 'rois': []}

            seq = roiDic[uniprotID]['sequence']
            roiSeq = Sequence(name=roiID, id=roiID, sequence=roiStr,
                              description=roiDesc)
            if roiIdxs[0]:
                roi = SequenceROI(sequence=seq, seqROI=roiSeq,
                                  roiIdx=roiIdxs[0], roiIdx2=roiIdxs[1])

                roiDic[uniprotID]['rois'].append(roi)
        return roiDic

