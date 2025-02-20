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
import json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.utils.utilsFasta import pairwiseAlign, parseFasta
from pwchem.objects import SequenceROI, SetOfSequenceROIs, Sequence
from pwchem.utils import *

class ProtDefineSeqROI(EMProtocol):
    """
    Defines a list of sequence ROIs, each of them from:\n'
        1) A residue or range of residues.\n'
        2) A predefined variant\n'
        3) One or several mutations
    """
    _label = 'Define sequence ROIs'
    _inputOptions = ['Sequence', 'SequenceVariants']
    _originOptions = ['Residues', 'SubSequences', 'Variant', 'Mutations']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('chooseInput', params.EnumParam, choices=self._inputOptions,
                       label='Define ROIs from: ', default=0, display=params.EnumParam.DISPLAY_HLIST,
                       help='Define sequence ROIs from a sequence or from a sequence variants object')
        group.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                      allowsNull=True, label="Input sequence: ", condition='chooseInput==0',
                      help='Select the sequence object where the ROI will be defined')
        group.addParam('inputSequenceVariants', params.PointerParam, pointerClass='SequenceVariants',
                       label='Input Sequence Variants:', condition='chooseInput==1', allowsNull=True,
                       help="Sequence containing the information about the variants and mutations")

        group = form.addGroup('Add ROI')
        group.addParam('whichToAdd', params.EnumParam, choices=self._originOptions,
                       display=params.EnumParam.DISPLAY_HLIST,
                       label='Add ROI from: ', default=0,
                       help='Add ROI from which definition (residues, variant or mutation)')
        #From residues
        group.addParam('resPosition', params.StringParam, label='Residues of interest: ',
                       condition='whichToAdd==0',
                       help='Specify the residue to define a region of interest.\n'
                            'You can either select a single residue or a range '
                            '(it will take into account the first and last residues selected)')

        # From a set of sequences (that must constain subsequences)
        group.addParam('inputSubsequences', params.PointerParam, label='SubSequences: ', condition='whichToAdd==1',
                       pointerClass='SetOfSequences', allowsNull=True,
                       help='Specify the set of sequences containing the subsequences that will be defined as ROIs')

        #From Variant
        group.addParam('selectVariant', params.StringParam, condition='chooseInput==1 and whichToAdd==2',
                       label='Select a predefined variant:',
                       help="Variant to use for defining the ROIs. Each mutation will be a different ROI")
        #From mutations
        group.addParam('selectMutation', params.StringParam,
                       label='Select some mutations: ', condition='chooseInput==1 and whichToAdd==3',
                       help="Mutations to be defined as sequence ROIs.\n"
                            "You can do multiple selection. Each mutation will be a different ROI")

        # Common for ROIs independent of the origin
        group.addParam('descrip', params.StringParam, default='',
                       label='ROI description: ', condition='whichToAdd in [0]',
                       help='Specify some description for this region of interest')

        group.addParam('addROI', params.LabelParam,
                       label='Add defined ROIs: ',
                       help='Add defined residues, variant or mutations to become a ROI')
        group.addParam('inROIs', params.TextParam, width=70, default='',
                      label='Input residues: ',
                      help='Input residues to define the ROI.')

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
        self._insertFunctionStep('defineOutputStep')

    def defineOutputStep(self):
        inpSeq = self.getInputSequence()
        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))

        residuesStr = self.inROIs.get().strip().split('\n')
        for rStr in residuesStr:
            roiInfo = ':'.join(rStr.split(':')[1:]).strip()

            # Residues origin
            if '{}:'.format(self._originOptions[0]) in rStr:
                resDic = json.loads(roiInfo)
                roiList, resIdxs = [resDic['residues']], resDic['index']
                idxsList = [[int(resIdxs.split('-')[0]), int(resIdxs.split('-')[1])]]
                descList = [resDic['desc']] if 'desc' in resDic else ['']

            elif '{}:'.format(self._originOptions[1]) in rStr:
                setIdx = json.loads(':'.join(rStr.split(':')[1:]))['PointerIdx']
                inSeqSet = self.inputPointers[int(setIdx)].get()
                roiList, idxsList, descList = self.roisFromSequences(inSeqSet)

            elif not 'Original' == roiInfo:
                # Variant origin
                if '{}:'.format(self._originOptions[2]) in rStr:
                    var2mutDic = self.inputSequenceVariants.get().getMutationsInLineage()
                    muts = var2mutDic[roiInfo]

                # Mutants origin
                elif '{}:'.format(self._originOptions[3]) in rStr:
                    muts = roiInfo.split(',')

                roiList, idxsList, descList = [], [], []
                for mut in muts:
                    mut = mut.strip()
                    roiList.append(mut[-1])
                    idxsList.append([int(mut[1:-1]), int(mut[1:-1])])
                    descList.append(mut)

            elif 'Original' == roiInfo:
                roiList = []

            for i in range(len(roiList)):
                roi, idxs, desc = roiList[i], idxsList[i], descList[i]
                roiSeq = Sequence(sequence=roi, name='ROI_{}-{}'.format(*idxs), id='ROI_{}-{}'.format(*idxs),
                                  description=desc)
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
        if not self.inputSequence.get() and not self.inputSequenceVariants.get():
            errors += ['You must specify an input Sequence or SequenceVariants']
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getInputSequence(self):
        if self.chooseInput.get() == 0:
            return self.inputSequence.get()
        elif self.chooseInput.get() == 1:
            return self.inputSequenceVariants.get()._sequence

    def parseROIFasta(self, fastaFile):
        roiStrList, roiIdxList = [], []
        fDic = parseFasta(fastaFile)
        roiLine = fDic[list(fDic.keys())[1]]
        for roiStr in roiLine.split('-'):
            if roiStr.strip():
                roiIdx = roiLine.index(roiStr)
                roiIdxs = [roiIdx + 1, roiIdx + len(roiLine) + 1]
                roiStrList.append(roiStr), roiIdxList.append(roiIdxs)

        return roiStrList, roiIdxList

    def roisFromSequences(self, seqs):
        roiList, idxsList, descList = [], [], []
        inputSeq = self.getInputSequence()
        for i, subseq in enumerate(seqs):
            subseqName = subseq.getSeqName()
            if not subseqName:
                subseqName = 'sequence_{}'.format(i)

            out_file = os.path.abspath(self._getExtraPath("pairWise_{}.fasta".format(subseqName)))
            pairwiseAlign(inputSeq, subseq, out_file, seqName1='Original', seqName2=subseqName)

            roiStrList, roiIdxsList = self.parseROIFasta(out_file)
            for i in range(len(roiStrList)):
                roiList.append(roiStrList[i]), idxsList.append(roiIdxsList[i]), descList.append(subseqName)
        return roiList, idxsList, descList


    # ADD WIZARD UTILS
    def getOriginLabel(self):
      if self.whichToAdd.get() == 0:
        inputLabel, sumLabel = 'resPosition', 'Residues'
      elif self.whichToAdd.get() == 2:
        inputLabel, sumLabel = 'selectVariant', 'Variant'
      elif self.whichToAdd.get() == 3:
        inputLabel, sumLabel = 'selectMutation', 'Mutations'
      elif self.whichToAdd.get() == 1:
        inputLabel, sumLabel = 'inputSubsequences', 'SubSequences'
      return inputLabel, sumLabel

    def getPrevPointersIds(self):
      ids = []
      for p in self.inputPointers:
        ids.append(p.get().getObjId())
      return ids

    def buildSumLine(self, type=None):
      inputLabel, sumLabel = self.getOriginLabel()

      if sumLabel == 'SubSequences':
        prevIds = self.getPrevPointersIds()
        newSet = self.inputSubsequences.get()
        newId = newSet.getObjId()

        prevPointers = self.inputPointers
        if newId not in prevIds:
          newIndex = len(prevPointers)
          prevPointers.append(Pointer(newSet))
        else:
          newIndex = prevIds.index(newId)
        setattr(self, 'inputPointers', prevPointers)

      roiInfo = getattr(self, inputLabel).get()
      if sumLabel == 'SubSequences':
        roiInfo = '{"PointerIdx": "%s", "Name": "%s"}' % (newIndex, roiInfo.__str__())
      elif self.descrip.get().strip():
        roiInfo = roiInfo.replace('}', ', "desc": "%s"}' % (self.descrip.get().strip()))
      return f'{sumLabel}: {roiInfo}'




