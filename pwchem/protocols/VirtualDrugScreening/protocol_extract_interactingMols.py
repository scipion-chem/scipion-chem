# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import numpy as np

from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from ...objects import SetOfSmallMolecules


class ProtExtractInteractingMols(EMProtocol):
  """Extract a subset of the interacting molecules stored in the set of sequences"""
  _label = 'extract interacting molecules'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequencesChem",
                    label='Input protein sequences: ',
                    help="Set of protein sequences containing interacting information (from ConPLex from example)")

    mGroup = form.addGroup('Filters')
    mGroup.addParam('chooseSeq', params.StringParam, label='Extract interacting molecules for protein: ', default='All',
                    help='Extract the interacting molecules for the selected protein sequence. '
                         'Use the wizard to get a list of the protein sequences in the input')
    mGroup.addParam('chooseMol', params.StringParam, label='Extract only from this subset of molecules: ',
                    default='All', expertLevel=params.LEVEL_ADVANCED,
                    help='Extract the only interacting molecules for the selected protein sequence in this subset. '
                         'Use the wizard to get a list of the interacting molecules in the input')
    mGroup.addParam('scThres', params.FloatParam, label='Score threshold: ', default=0.3,
                    help='Score threshold to filter interacting molecules below it for the selected protein sequences')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.createOutputStep)

  def createOutputStep(self):
    inSeqs = self.inputSequences.get()
    _, _, molNames = self.getFilteredOutput()

    seqName = None
    inSeqNames = self.chooseSeq.get().split(',')
    if len(inSeqNames) == 1 and inSeqNames[0].strip() != 'All':
      seqName = inSeqNames[0].strip()

    outMols = SetOfSmallMolecules().create(outputPath=self._getPath())
    for mol in self.getInputMols():
      molName = mol.getMolName()
      if molName in molNames:
          if seqName:
            mol._interactScore = pwobj.Float(inSeqs.getInteractScoresDic()[seqName][molName])
          outMols.append(mol)

    self._defineOutputs(outputSmallMolecules=outMols)


  ############## UTILS ########################
  def getInputMols(self):
    return self.inputSequences.get().getInteractMols()

  def getFilteredOutput(self):
    inSeqs = self.inputSequences.get()
    intDic = inSeqs.getInteractScoresDic()

    seqNames, molNames = inSeqs.getSequenceNames(), inSeqs.getInteractMolNames()
    seqNames, molNames = self.filterNames(seqNames, molNames)

    intAr = self.formatInteractionsArray(intDic, seqNames, molNames)
    intAr, seqNames, molNames = self.filterScores(intAr, seqNames, molNames)
    return intAr, seqNames, molNames

  def filterNames(self, seqNames, molNames):
    inSeqNames = self.chooseSeq.get().split(',')
    if 'All' not in inSeqNames:
      seqNames = [seqName for seqName in seqNames if seqName in inSeqNames]

    inMolNames = self.chooseMol.get().split(',')
    if 'All' not in inMolNames:
      molNames = [molName for molName in molNames if molName in inMolNames]

    return seqNames, molNames

  def filterScores(self, intAr, seqNames, molNames):
    ips, ims = [], []
    scThres = self.scThres.get()

    for ip, seqName in enumerate(seqNames):
      if any(intAr[ip, :] >= scThres):
        ips.append(ip)

    for im, molName in enumerate(molNames):
      if any(intAr[:, im] >= scThres):
        ims.append(im)

    if not len(seqNames) == len(ips):
      seqNames = list(np.array(seqNames)[ips])
      intAr = intAr[ips, :]

    if not len(molNames) == len(ims):
      molNames = list(np.array(molNames)[ims])
      intAr = intAr[:, ims]

    return intAr, seqNames, molNames

  def formatInteractionsArray(self, intDic, seqNames, molNames):
    intAr = np.zeros((len(seqNames), len(molNames)))
    for i, seqName in enumerate(seqNames):
      for j, molName in enumerate(molNames):
        intAr[i, j] = intDic[seqName][molName]
    return intAr


