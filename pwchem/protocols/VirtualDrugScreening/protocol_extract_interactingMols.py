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

from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules, SmallMoleculesLibrary
from pwchem.utils import getFilteredOutput


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

  def defineMolsOutput(self, intDic, molNames, seqName):
    outMols = SetOfSmallMolecules().create(outputPath=self._getPath())
    for mol in self.getInputMols():
      molName = mol.getMolName()
      if molName in molNames:
        if seqName:
          mol._interactScore = pwobj.Float(intDic[seqName][molName])
        outMols.append(mol)

    self._defineOutputs(outputSmallMolecules=outMols)

  def defineLibraryOutput(self, intMols, intDic, molNames, seqName):
    inFile, oFile = intMols.getFileName(), self._getPath('outputLibrary.smi')
    with open(inFile) as fIn:
      with open(oFile, 'w') as fO:
        for line in fIn:
          smi, smiName = line.split()[0].strip(), line.split()[1].strip()
          if smiName in molNames:
            if seqName:
              score = intDic[seqName][smiName]
            fO.write(f'{smi}\t{smiName}\t{score}\n')

    intMols.setFileName(oFile)
    intMols.calculateLength()
    self._defineOutputs(outputLibrary=intMols)

  def createOutputStep(self):
    inSeqs = self.inputSequences.get()
    filtSeqNames, filtMolNames = self.chooseSeq.get().strip().split(','), self.chooseMol.get().strip().split(',')
    molNames = getFilteredOutput(inSeqs, filtSeqNames, filtMolNames, self.scThres.get())[-1]

    seqName = None
    if len(filtSeqNames) == 1 and filtSeqNames[0].strip() != 'All':
      seqName = filtSeqNames[0].strip()

    intMols = self.getInputMols()
    intDic = inSeqs.getInteractScoresDic()
    if isinstance(intMols, SetOfSmallMolecules):
      self.defineMolsOutput(intDic, molNames, seqName)

    elif isinstance(intMols, SmallMoleculesLibrary):
      self.defineLibraryOutput(intMols, intDic, molNames, seqName)


  ############## UTILS ########################
  def getInputMols(self):
    return self.inputSequences.get().getInteractMols()

