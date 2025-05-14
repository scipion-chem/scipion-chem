# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
# *             Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
This protocol is used to rank a set of docked molecules coming from different origins.
These molecules could have been docked using several different protocols.
The scores will be combined for the molecules, depending on their ranking.

"""

import os, json

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects import Float

from pwchem.objects import SetOfSmallMolecules


def normalizeToRange(iterable, normRange=[0, 1]):
  maxIt, minIt = max(iterable), min(iterable)
  return [((normRange[1] - normRange[0]) * (i - minIt)) / (maxIt - minIt) + normRange[0] for i in iterable]

class ProtocolRankDocking(EMProtocol):
    """
    Executes the rank scoring to combine different origin docked molecules scores.
    The rank score will be defined by the ranking from the different docking scores for each set.
    """
    _label = 'Rank docking score'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputMoleculesSets', params.MultiPointerParam, pointerClass='SetOfSmallMolecules',
                        label="Input Docked Small Molecules: ",
                        help='Select the docked molecules to be ranked')
        group.addParam('useScore', params.BooleanParam, label="Weight rank by item score: ", default=False,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Use the selected score to weight the contribution of each vote.'
                            'The scores are normalized prior to be used as weight')

        group = form.addGroup('Define scores')
        group.addParam('defineScores', params.BooleanParam, label="Define scores: ", default=False,
                       help='Define the scores used from each set and its weight.\n'
                            'If true, the inputs and scores in the summary list will be used. '
                            'Inputs not present in the list will be treated as default.\n'
                            'By default, the inputs in the multipointer param, _energy attribute, '
                            ' weight 1 and small values are better will be used.')
        group.addParam('defineInput', params.StringParam, label="Select input: ", default='', condition='defineScores',
                       help='Select an input to define its score')
        group.addParam('defineScore', params.StringParam, label="Select score: ", default='', condition='defineScores',
                       help='Select a score for the selected input to be used')
        group.addParam('defineWeight', params.FloatParam, label="Select weight: ",
                       default=1.0, condition='defineScores',
                       help='Select a weight for the input with respect to the others defined')
        group.addParam('defineDirection', params.BooleanParam, label="Small is good?: ",
                       default=True, condition='defineScores',
                       help='Define the direction of the score, whether small values are prefered')

        group.addParam('defineSummary', params.TextParam, width=70, label="Scores summary: ",
                       default='', condition='defineScores',
                       help='Summary of the inputs that will be used, their respective scores and their weights')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.createOutputStep)

    def rankVoting(self, mols, scoreAttribute, smallIsGood=True):
        scDic = {mol.clone(): getattr(mol, scoreAttribute).get() for mol in mols}
        scDic = dict(sorted(scDic.items(), key=lambda x:x[1], reverse=not smallIsGood))
        if self.useScore.get():
          normSc = normalizeToRange(scDic.values())
          for i, mol in enumerate(scDic):
            scDic[mol] = normSc[i]

        # The better the position in the ranking and the more appereances of a molecule in it will generate a higher
        # voting score.
        vote = {mol.clone().getMolName(): 0.0 for mol in mols}
        nDocks = len(scDic)
        for i, (mol, score) in enumerate(scDic.items()):
          curVote = 1.0 - i/nDocks
          if self.useScore.get():
            score = 1 - score if smallIsGood else score
            curVote *= score
          vote[mol.getMolName()] += curVote

        # Normalize by the total number of dockings
        for molName, score in vote.items():
            vote[molName] = score / nDocks
        return vote

    def getBestMols(self, molSet, attr='_energy', smallIsGood=True, bestDic={}):
      for mol in molSet:
        energy = getattr(mol, attr)
        molName = mol.getMolName()
        if molName not in bestDic or \
                (smallIsGood and energy < getattr(bestDic[molName], attr)) or\
                (not smallIsGood and energy > getattr(bestDic[molName], attr)):
          bestDic[molName] = mol.clone()
      return bestDic

    def createOutputStep(self):
        self.voteDic, bestMols = {}, {}
        inpSum = self.getInputSummary()
        for inpIdx, inpDic in inpSum.items():
            molSet = self.inputMoleculesSets[inpIdx].get()
            bestMols = self.getBestMols(molSet, inpDic['Score'], inpDic['Small'], bestDic=bestMols)

            curVote = self.rankVoting(molSet, inpDic['Score'], inpDic['Small'])
            for molName, score in curVote.items():
                if molName in self.voteDic:
                    self.voteDic[molName] += score * inpDic['Weight']
                else:
                    self.voteDic[molName] = score * inpDic['Weight']

        with open(self.getResultsFile(), 'w') as f:
            for molName, score in self.voteDic.items():
                f.write(f'{molName}\t{score}\n')

        outMols = SetOfSmallMolecules(filename=self._getPath('outputSmallMolecules.sqlite'))
        bestMols = self.fillEmptyAttributes(bestMols.values())

        usedIds = []
        for mol in bestMols:
          setattr(mol, 'rankScore', Float(self.voteDic[mol.getMolName()]))
          molId = mol.getObjId()
          if molId in usedIds:
            molId = self.getMinObjId(usedIds)
            mol.setObjId(molId)
          outMols.append(mol)
          usedIds.append(molId)

        outMols.proteinFile.set(self.getOriginalReceptorFile())
        outMols.setDocked(True)
        outMols.saveGroupIndexes()
        self._defineOutputs(outputSmallMolecules=outMols)

    def getOriginalReceptorFile(self):
        return self.inputMoleculesSets[0].get().getProteinFile()

    def fillEmptyAttributes(self, inSet):
        '''Fill all items with empty attributes'''
        attributes = self.getAllAttributes(self.inputMoleculesSets)
        for item in inSet:
            for attr in attributes:
                if not hasattr(item, attr):
                    item.__setattr__(attr, attributes[attr])
        return inSet

    def getAllAttributes(self, inSets):
        '''Return a dic with {attrName: ScipionObj=None}'''
        attributes = {}
        for inSet in inSets:
            item = inSet.get().getFirstItem()
            attrKeys = item.getObjDict().keys()
            for attrK in attrKeys:
                if not attrK in attributes:
                    value = item.__getattribute__(attrK)
                    attributes[attrK] = value.clone()
                    attributes[attrK].set(None)
        return attributes

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if os.path.exists(self.getResultsFile()):
            resDic = self.parseResults()
            scDic = dict(sorted(resDic.items(), key=lambda x: x[1], reverse=True))
            summary += ['Mol Name\tRanking score']
            summary += [f'{molName} {score}' for molName, score in scDic.items()]
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def parseResults(self):
      resDic = {}
      resFile = self.getResultsFile()
      with open(resFile) as f:
          for line in f:
              resDic[line.split()[0]] = float(line.split()[1])
      return resDic

    def getMinObjId(self, doneObjIds):
      '''Return the int not included in the input list'''
      maxListed = max(doneObjIds)
      notListed = set(range(1, maxListed)).difference(set(doneObjIds))
      return min(notListed) if len(notListed) > 0 else maxListed + 1

    def getInputSummary(self):
      inDic = {}
      if self.defineScores.get():
        for scoreLine in self.defineSummary.get().split('\n'):
          if scoreLine.strip():
            lDic = json.loads(scoreLine)
            inDic[lDic['InputIndex']] = {'Score': lDic['Score'], 'Weight': lDic['Weight'],
                                         'Small': bool(lDic['Small'])}

      for i, inPointer in enumerate(self.inputMoleculesSets):
        if i not in inDic:
          inDic[i] = {'Score': '_energy', 'Weight': 1.0, 'Small': True}
      return inDic


    def createElementLine(self):
        inputIndex = self.defineInput.get().split('//')[0]
        sumStr = f'{{"InputIndex": {inputIndex}, ' \
                 f'"Score": "{self.defineScore.get()}", "Weight": {self.defineWeight.get()}, ' \
                 f'"Small": "{self.defineDirection.get()}"}}\n'
        return sumStr

    def getResultsFile(self):
      return self._getPath('results.tsv')






