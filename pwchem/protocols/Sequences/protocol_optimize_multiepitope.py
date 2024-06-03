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
This protocol is used to optimize a multiepitope object based on some predictor scores using genetic algorithms
"""
import numpy, json, random, time, os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.utils import runInParallel
from pwchem.utils.utilsDEAP import *

from deap import base, creator, tools, algorithms

REM, ADD, MOD = 0, 1, 2

def ensureMinSize(mEpList, minSize):
  zerosIdxs = numpy.where(mEpList == -1)
  if minSize > len(mEpList) - len(zerosIdxs):
    nActive = len(mEpList) - len(zerosIdxs)
    nActivate = minSize - nActive
    activateIdxs = numpy.random.choice(zerosIdxs, nActivate, False)
    mEpList[activateIdxs] = 1
  return mEpList

def buildMultiEpitopeZeros(epScores, pZero, maxSize, minSize):
  '''Build a list of epitope indexes representing a multiepitope.
  - epScores: weight for each of the scores to sample from randomly
  - pZero: probability of zero epitope (-1) (this way, we get different number of epitopes)
  - maxSize: maximum number of epitopes (including zeroEpitope) to form the multiEpitope
  '''
  # Populate multiepitope with -1 (no epitopes)
  multiEpList = numpy.random.choice([-1, 1], maxSize, True, [pZero, 1 - pZero])
  multiEpList = ensureMinSize(multiEpList, minSize)

  # Randomly choose the necessary number of epitopes, weighting them by the scores
  epIdxs, epScores = list(epScores.keys()), list(epScores.values())
  epScores = numpy.array(epScores) / sum(epScores)
  nEps = numpy.count_nonzero(multiEpList == 1)
  epitopes = numpy.random.choice(epIdxs, nEps, False, epScores)

  # Populate multiepitope with chosen epitopes
  i = 0
  for j in range(len(multiEpList)):
    if multiEpList[j] == 1:
      multiEpList[j] = epitopes[i]
      i += 1
  return multiEpList

def mutDifferentIdx(ind, low, high, indpb):
  '''Mutates an individual by chaging on of the non-zero ints in their sequence
  by another in range(low, high+1) which is not already present in the sequence'''
  for i in range(len(ind)):
    if random.random() < indpb:
      possibleIdxs = set(range(low, high + 1)) - set(ind)
      if len(possibleIdxs) == 0:
        raise IndexError("Individual cannot be mutated because there are no index options left")
      ind[i] = random.sample(sorted(possibleIdxs) + [-1], 1)[0]
  return ind,

def ensureDifferentIds(nEpitopes):
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      for child in offspring:
        idxs = []
        for i in range(len(child)):
          idx = child[i]
          if idx != -1:
            if idx in idxs:
              possibleIdxs = set(range(nEpitopes)) - set(idxs)
              child[i] = random.sample(sorted(possibleIdxs), 1)[0]
              idxs.append(child[i])
            else:
              idxs.append(idx)
      return offspring

    return wrapper

  return decorator

def checkSameInd(ind1, ind2):
  '''Check that two multiEpitopes are not formed by the same epitopes'''
  cInd1, cInd2 = ind1.copy(), ind2.copy()
  cInd1.sort(), cInd2.sort()
  return cInd1 == cInd2

def reportPoolStatus(poolDic):
  '''Check the status of the AsynPool objects stored as values of the dictionary and reports when they finish
  '''
  ready = []
  while len(ready) < len(poolDic):
    time.sleep(5)
    for evalSoft, po in poolDic.items():
      if po.ready() and evalSoft not in ready:
        ready.append(evalSoft)
        print(f'{evalSoft} execution finished ({len(ready)} / {len(poolDic)})')

SIMP, MUPLUS, MUCOMMA = 0, 1, 2
P1, P2, UNI = 0, 1, 2
TOUR, ROUL, BEST = 0, 1, 2

IEDB, IIITD, DDG = 'IEDB-Coverage', 'IIITD', 'DDG'

class ProtOptimizeMultiEpitope(EMProtocol):
    """
    Optimize a MultiEpitope object based on some scores using Genetic Algorithms
    """
    _label = 'Optimize multiepitope'

    _agTypes = ['eaSimple', 'eaMuPlusLambda', 'eaMuCommaLambda']
    # _mutTypes = ['Random', 'Shuffle']
    _crossTypes = ['One point', 'Two points', 'Uniform']
    _selectionTypes = ['Tournament', 'Roulette', 'Best']

    _evalParameters = {IEDB: ['mhc', 'pop'],
                       IIITD: ['chooseIIITDEvaluator'],
                       DDG: ['chooseDDGEvaluator'],
                      }

    _gaAlgorithms = {SIMP: myEaSimple, MUPLUS: myEaMuPlusLambda, MUCOMMA: myEaMuCommaLambda}
    _cxAlgorithms = {P1: tools.cxOnePoint, P2: tools.cxTwoPoint, UNI: tools.cxUniform}
    _selAlgorithms = {TOUR: tools.selTournament, ROUL: tools.selRoulette, BEST: tools.selBest}

    _indMemory = {}

    # -------------------------- DEFINE param functions ----------------------
    def getEvaluationProtocols(self):
      origins = {}
      try:
        from iedb.protocols import ProtMHCIIPopulationCoverage
        origins['IEDB - Coverage'] = ProtMHCIIPopulationCoverage
      except:
        pass

      try:
        from iiitd.protocols import ProtIIITDEvaluations
        origins['IIITD'] = ProtIIITDEvaluations
      except:
        pass

      try:
        from ddg.protocols import ProtDDGEvaluations
        origins['DDG'] = ProtDDGEvaluations
      except:
        pass

      return origins

    def getEvaluationOrigins(self):
      oriDic = self.getEvaluationProtocols()
      return list(oriDic.keys())

    def _defineMHCEvalParams(self, evalGroup, allCond):
      evalGroup.addParam('mhc', params.EnumParam, label='MHC type: ', display=params.EnumParam.DISPLAY_HLIST,
                         default=2, choices=['I', 'II', 'combined'], condition=allCond,
                         help="Type of MHC fpor the prediction.")

      evalGroup.addParam('pop', params.StringParam, label='Population: ', default='Area', condition=allCond,
                         help="Select population by geographic area or ethnicity. If several, comma-separated. "
                              "For all the geographic areas, use 'Area', for ethnicities use 'Ethnicity'."
                              "For the complete list of populations, you can check http://tools.iedb.org/population/")
      return evalGroup

    def _defineEvalParams(self, evalGroup, originParam):
      evalProtDic = self.getEvaluationProtocols()
      for i, (evalLabel, evalProtocol) in enumerate(evalProtDic.items()):
        allCond = f'{originParam}=={i}'
        if evalLabel == 'IEDB - Coverage':
          self._defineMHCEvalParams(evalGroup, allCond)
        elif hasattr(evalProtocol, '_defineEvalParams'):
          evalProtocol()._defineEvalParams(evalGroup, allCond=allCond)
      return evalGroup


    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label="Input sequence ROIs: ",
                       help='Select the Set of sequence ROIs where to sample the epitopes')

        group = form.addGroup('Multiepitope definition')
        line = group.addLine('Number of epitopes: ', important=True, help='Range of possible number of epitopes')
        line.addParam('minEp', params.IntParam, label='Min: ', default=3)
        line.addParam('maxEp', params.IntParam, label='Max: ', default=10)
        group.addParam('pZero', params.FloatParam, label='Zero probability: ', default=0.3,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Probability of empty epitope when building a multiepitope. The higher the probability,'
                            'newly created multiepitopes will more likely have smaller number of epitopes')

        # line = group.addLine('Number of residues: ', help='Range of possible total number of residues')
        # line.addParam('minRes', params.IntParam, label='Min: ', default=10)
        # line.addParam('maxRes', params.IntParam, label='Max: ', default=1000)


        form.addSection(label='Epitope sampling')
        group = form.addGroup('Epitope sampling scores')
        group.addParam('inScore', params.StringParam, label='Select a score: ', default='', important=True,
                       help='Choose a score to define the epitope sampling probability')
        group.addParam('scoreWeight', params.FloatParam, label='Score weight: ', default=1, important=True,
                       help='Define a weight for the score. The bigger the absolute value of the weight, the more '
                            'important it will be. If negative, the filter will look for smaller values.')
        group.addParam('addScore', params.LabelParam, label='Add score: ',
                       help='Add defined score')

        group.addParam('scoreSummary', params.TextParam, width=100, label='Score summary:', default='',
                       help='Summary of the epitope sampling scores')

        group.addParam('normScores', params.BooleanParam, label='Normalize scores: ', default=True,
                       expertLevel=params.LEVEL_ADVANCED,
                       help='Normalize defined scores in range [0, 1]')

        group = form.addGroup('Manual sampling')
        group.addParam('inROI', params.StringParam, label='Select a ROI: ', default='',
                       help='Choose a ROI to be included manually as epitope')
        group.addParam('manualProb', params.FloatParam, label='Manual weight: ', default=1,
                       help='Define the probability of including the selected epitope in the multiepitope.'
                            'If 1, the epitope will always be included')
        group.addParam('addManual', params.LabelParam, label='Add score: ',
                       help='Add defined epitope manual sampling')
        group.addParam('manualSummary', params.TextParam, width=100, label='Manual summary:', default='',
                       help='Summary of the manual sampling epitopes. They will be added randomly independently of '
                            'the sampling epitopes based on scores')

        form.addSection(label='Linkers')
        group = form.addGroup('Multiepitope linkers')
        # group.addParam('linkerAction', params.EnumParam, label='Linker strategy: ', default=0,
        #                choices=['Random', 'Per input Set'],
        #                help='Choose a strategy for the linkers to be chosen')
        group.addParam('inLinker', params.StringParam, label='Define linker sequence: ', default='',
                       help='Define the sequence of the linker to be used')
        # group.addParam('linkerWeight', params.FloatParam, label='Linker weight: ', default=1,
        #                help='Weight for the defined linker')
        # group.addParam('addLinker', params.LabelParam, label='Add linker: ',
        #                help='Add defined linker')
        # group.addParam('linkerSummary', params.TextParam, width=100, label='Linker summary:', default='',
        #                help='Summary of the defined linkers.')

        form.addSection(label='Multiepitope evaluations')
        group = form.addGroup('Multiepitope evaluations')
        group.addParam('multiEval', params.EnumParam, label='Select evaluation origin: ', default=0,
                       choices=self.getEvaluationOrigins(), important=True,
                       help='Select the multiepitope evaluation to add')
        self._defineEvalParams(group, originParam='multiEval')

        group.addParam('evalWeight', params.FloatParam, label='Evaluation weight: ', default=1,
                       help='Define a weight for the evaluation. The bigger the absolute value of the weight, the more '
                            'important it will be.')
        group.addParam('addEval', params.LabelParam, label='Add conditional filter: ',
                       help='Add defined conditional filter')

        group = form.addGroup('Evaluations summary')
        group.addParam('evalSummary', params.TextParam, width=100, label='Evaluations summary:', default='',
                       help='Summary of the selected conditional filters')

        form.addSection(label='Genetic Algorithm')
        group = form.addGroup('General parameters')
        group.addParam('nGen', params.IntParam, label='Number of generations: ', default=10, important=True,
                       help='Number of generations for the GA')
        group.addParam('nPop', params.IntParam, label='Population size: ', default=20, important=True,
                       help='Number of individuals in the population (per generation)')

        group.addParam('gaType', params.EnumParam, label='GA type: ', default=SIMP, choices=self._agTypes,
                       help='Choose a strategy for the linkers to be chosen. For more information: '
                            'https://deap.readthedocs.io/en/master/api/algo.html#complete-algorithms')
        group.addParam('varProp', params.FloatParam, label='Proportion of varied individuals: ', default=0.25,
                       condition='gaType!=0',
                       help='Defines the number of varied individuals as the proportion with the number of individuals '
                            'in the population.')

        group = form.addGroup('Mutation')
        # group.addParam('mutType', params.EnumParam, label='Mutation type: ', default=0, choices=self._mutTypes,
        #                help='Choose a strategy for the mutation:'
        #                     '\n- Random will replace an epitope for a random sampled from the input.'
        #                     '\n- Shuffle will exchange the position of an epitope with another present in the '
        #                     'individual')
        group.addParam('mutProb', params.FloatParam, label='Mutate probability: ', default=0.4,
                       help='Probability for a individual to be mutated for the next generation.')
        group.addParam('mutIndProb', params.FloatParam, label='Epitope mutation probability: ', default=0.25,
                      help='Probability for each epitope in the individual to be mutated independently')

        group = form.addGroup('Crossover')
        group.addParam('crossType', params.EnumParam, label='Crossover type: ', default=0, choices=self._crossTypes,
                       help='Choose a strategy for the crossover. For more information: '
                            'https://deap.readthedocs.io/en/master/api/tools.html#crossover')
        group.addParam('crossProb', params.FloatParam, label='Crossover probability: ', default=0.4,
                       help='Probability for a pair of individuals to be crossed for the next generation.')
        group.addParam('crossIndProb', params.FloatParam, label='Epitope crossover probability: ', default=0.25,
                       condition='crossType==2',
                       help='Probability for each epitope in the individual to be crossed independently')

        group = form.addGroup('Selection')
        group.addParam('selType', params.EnumParam, label='Selection type: ', default=0, choices=self._selectionTypes,
                       help='Choose a strategy for the selection. For more information: '
                            'https://deap.readthedocs.io/en/master/api/tools.html#selection')
        group.addParam('tournSize', params.IntParam, label='Tournament size: ', default=3, condition='selType==0',
                       help='The number of individuals participating in each tournament')

        group = form.addGroup('Hall of Fame')
        group.addParam('hallSize', params.IntParam, label='Hall of fame size: ', default=5, important=True,
                       help='The number of best scored individuals that will be saved outside of the population')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
      self._insertFunctionStep(self.optimizeMultiepitopeStep)
      self._insertFunctionStep(self.defineOutputStep)

    def optimizeMultiepitopeStep(self):
      self.inputROISeqs = {roi.getObjId(): roi.getROISequence() for roi in self.inputROIs.get()}
      deapToolbox = self.defineDeapToolbox()

      nPop, nGen = self.nPop.get(), self.nGen.get()
      mutProb, crossProb = self.mutProb.get(), self.crossProb.get()

      pop = deapToolbox.population(n=nPop)
      hallOfFame = tools.HallOfFame(self.hallSize.get(), checkSameInd)

      kwargs = {'cxpb': crossProb, 'mutpb': mutProb, 'ngen': nGen, 'halloffame': hallOfFame, 'verbose': False}
      if self.gaType.get() != SIMP:
        lamb = round(nPop * self.varProp.get())
        kwargs.update({'mu': nPop, 'lambda_': int(lamb)})

      gaAlgorithm = self._gaAlgorithms[self.gaType.get()]
      pop, info = gaAlgorithm(pop, deapToolbox, **kwargs)
      print('Final pop: ', pop)
      print('Scores: ', [ind.fitness.values for ind in pop])
      print('Fame: ', hallOfFame.items)
      print('Scores: ', [ind.fitness.values for ind in hallOfFame.items])



    def defineOutputStep(self):
      # todo: make functional: evaluations should be callable to extract the scores from another protocol
      pass

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        if len(self.inputROIs.get()) < self.maxEp.get():
          errors.append('The maximum number of epitopes cannot be higher than the number of input epitopes')
        evalDics = self.getEvalDics()
        if IEDB in evalDics:
          mhcs = [self.getEnumText('mhc')]
          mhcs = ['I', 'II'] if mhcs == 'combined' else mhcs
          for mhc in mhcs:
            if not hasattr(self.inputROIs.get().getFirstItem(), f'_allelesMHC{mhc}'):
              errors.append(f'{IEDB} analysis cannot be performed because no information about the MHC{mhc} alleles '
                            f'is found')

        if self.gaType.get() == MUCOMMA and self.varProp.get() < 1:
          errors.append(f'The Proportion of varied individuals {self.varProp.get()} must be >= 1 '
                        f'for the {MUCOMMA} GA type')
        return errors

    def _warnings(self):
        ws = []
        return ws

    # --------------------------- DEAP functions -----------------------------------
    def defineDeapToolbox(self):
      epScores = self.getEpitopeScores()
      indOptions = len(self.inputROIs.get())

      toolbox = base.Toolbox()

      creator.create("FitnessMin", base.Fitness, weights=(1.0,))
      creator.create("Individual", list, fitness=creator.FitnessMin)

      toolbox.register("indices", buildMultiEpitopeZeros, epScores, self.pZero.get(),
                       self.maxEp.get(), self.minEp.get())
      toolbox.register("individual", tools.initIterate, creator.Individual, toolbox.indices)
      toolbox.register("evaluate", self.multiEpitopeEval)

      toolbox.register("population", tools.initRepeat, list, toolbox.individual)

      toolbox = self.registerMate(toolbox)
      toolbox.decorate("mate", ensureDifferentIds(indOptions))
      toolbox.register("mutate", mutDifferentIdx, low=1, high=indOptions, indpb=self.mutIndProb.get())
      toolbox.decorate("mutate", ensureDifferentIds(indOptions))
      toolbox = self.registerSelection(toolbox)
      return toolbox

    def registerMate(self, toolbox):
      if self.crossType.get() == P1:
        toolbox.register("mate", tools.cxOnePoint)
      elif self.crossType.get() == P2:
        toolbox.register("mate", tools.cxTwoPoint)
      elif self.crossType.get() == UNI:
        toolbox.register("mate", tools.cxUniform, indpb=self.crossIndProb.get())
      return toolbox

    def registerSelection(self, toolbox):
      if self.selType.get() == TOUR:
        toolbox.register("select", tools.selTournament, tournsize=self.tournSize.get())
      elif self.selType.get() == TOUR:
        toolbox.register("select", tools.selRoulette)
      elif self.selType.get() == BEST:
        toolbox.register("select", tools.selBest)
      return toolbox

    def getMultiEpitopeSeqs(self, individuals):
      '''Returns the set of multiepitope sequences codified into a set of DEAP individuals'''
      seqs = {}
      for i, ind in enumerate(individuals):
        curSeqs = []
        for roiIdx in ind:
          if roiIdx != -1:
            curSeqs += [self.inputROISeqs[roiIdx]]

        links = self.getLinkers(len(curSeqs))
        seqs[i] = self.mergeEpsLinks(curSeqs, links)
      return seqs


    # --------------------------- WIZARD functions -----------------------------------
    def buildScoreSumLine(self):
      sLine = f'{{"Score": "{self.inScore.get()}", "Weight": {self.scoreWeight.get()}}}'
      return sLine

    def buildManualSumLine(self):
      sLine = f'{{"Epitope": "{self.inROI.get()}", "Prob": {self.manualProb.get()}}}'
      return sLine

    def buildLinkerSumLine(self):
      sLine = f'{{"Linker": "{self.inLinker.get()}", "Weight": {self.linkerWeight.get()}}}'
      return sLine

    def buildEvalSumLine(self):
      evalParams = self.getEvalParams()
      sLine = f'{{"Evaluation": "{self.getEnumText("multiEval")}", "Weight": {self.evalWeight.get()}'
      for parName in evalParams:
        sLine += f', "{parName}": "{self.getParamValue(parName)}"'
      return sLine + '}'

    # --------------------------- EVALUATION functions -----------------------------------
    def multiEpitopeEval(self, individuals):
      '''Performs the evaluations of a set of DEAP individuals.
      To be used in batch with a modified GA algoritms that allows so.
      Returns a list of individual fitnesses: [(f1,), ...]
      '''
      resultsDic = self.performEvaluations(individuals)
      weigths = self.getEvalWeights()
      fitnesses = numpy.zeros(len(individuals))
      for (evalKey, softName), res in resultsDic.items():
        scores = numpy.array(list(map(float, res)))
        fitnesses += scores * weigths[evalKey]

      return [(f,) for f in fitnesses]

    def performEvaluations(self, individuals):
      '''Generalize caller to the evaluation functions.
      - individuals: list of DEAP individuals -> list of ROI indexes
      - evalDics: dictionary as {(evalKey, softwareName): {parameterName: parameterValue}}
      - jobs: int, number of jobs for parallelization
      Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
      '''
      jobs = self.numberOfThreads.get()
      evalDics = self.getEvalDics()
      sequences = self.getMultiEpitopeSeqs(individuals)

      resultsDic = {}
      disJobs = self.getDistributedJobs(evalDics, jobs)
      for source, sDics in evalDics.items():
        nt = disJobs[source]
        if source == IEDB:
          epiDic = self.performCoverageAnalysis(individuals, sDics, nt)

        elif source == IIITD:
          from iiitd import Plugin as iiitdPlugin
          epiDic = iiitdPlugin.performEvaluations(sequences, sDics, nt, iiitdPlugin.getBrowserData(), verbose=False)

        elif source == DDG:
          from ddg import Plugin as ddgPlugin
          epiDic = ddgPlugin.performEvaluations(sequences, sDics, nt, ddgPlugin.getBrowserData(), verbose=False)

        else:
          continue

        resultsDic.update(epiDic)
      return resultsDic

    def performCoverageAnalysis(self, individuals, sDics, nt):
      from iedb import Plugin as iedbPlugin
      from iedb.utils import translateArea, buildMHCCoverageArgs, parseCoverageResults

      resDic = {}
      # todo: make functional
      inROIs = self.inputROIs.get()
      for (evalKey, soft), argDic in sDics.items():
        epFile = self._getExtraPath(f'inputEpitopes_{evalKey}.tsv')
        oDir = os.path.abspath(self._getExtraPath())

        selPops = [p.strip() for p in argDic['pop'].split(',')]
        selPops = translateArea(selPops)

        coveArgs = buildMHCCoverageArgs(inROIs, epFile, selPops, argDic['mhc'], oDir, separated=False)
        runInParallel(iedbPlugin.runPopulationCoverage, None,
                      paramList=[coveArg for coveArg in coveArgs], jobs=nt)

        oFile = os.path.join(oDir, f'inputEpitopes_{evalKey}_All_results.tsv')
        coveDic = parseCoverageResults(oFile)['average']
        resDic[(evalKey, soft)] = {'Score': coveDic['coverage'], 'outFile': oFile}
      return resDic

    def getEvalDics(self):
      '''Returns a dictionary for each of the sources chosen in the evaluations (IEDB, IIITD, DDG)
      This dictionary are of the form: {source: {(evalKey, softwareName): {parameterName: parameterValue}}}
      '''
      mapFuncs = {}

      evalDics = {}
      for i, sline in enumerate(self.evalSummary.get().split('\n')):
        if sline.strip():
          sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
          source, weight = sDic.pop('Evaluation'), sDic.pop('Weight')
          if not source in evalDics:
            evalDics[source] = {}

          if source == IEDB:
            evalKey, softName = f'{source}_{i+1}', source

          elif source == IIITD:
            softName = sDic.pop("chooseIIITDEvaluator")
            evalKey = f'{source}-{softName}_{i + 1}'
            sDic.update({'software': softName})
            if not source in mapFuncs:
              from iiitd.utils import mapEvalParamNames as iiitdMap
              mapFuncs[source] = iiitdMap

          elif source == DDG:
            softName = sDic.pop("chooseDDGEvaluator")
            evalKey = f'{source}-{softName}_{i + 1}'
            sDic.update({'software': softName})
            if not source in mapFuncs:
              from ddg.utils import mapEvalParamNames as ddgMap
              mapFuncs[source] = ddgMap

          else:
            continue

          evalDics[source].update({evalKey: sDic})

      for source, evalDic in evalDics.items():
        evalDics[source] = mapFuncs[source](evalDic)

      return evalDics

    # --------------------------- PARAMETERS functions -----------------------------------
    def getParamValue(self, paramName):
      if isinstance(self.getParam(paramName), params.EnumParam):
        value = self.getEnumText(paramName)
      else:
        value = getattr(self, paramName).get()
      return value

    def getDistributedJobs(self, evalDics, totalJobs):
      jobs = {source: 0 for source in evalDics}
      for source, sDics in evalDics.items():
        if source == IEDB:
          # not parallelized
          jobs[source] = 1
        else:
          jobs[source] = len(sDics)

      i=0
      sources = list(jobs.keys())
      # Decrease the number of threads per source until the totaljobs is met or all 1
      while sum(jobs.values()) > totalJobs and not all(jobs.values()) == 1:
        if jobs[sources[i % len(sources)]] > 1:
          jobs[sources[i % len(sources)]] -= 1
      return jobs

    def getEvalParams(self):
      '''Return the parameters involved in a defined evaluation'''
      evalOrigin = self.getEnumText("multiEval")
      params = [p for p in self._evalParameters[evalOrigin]]
      if evalOrigin == IIITD:
        try:
          from iiitd.protocols import ProtIIITDEvaluations
          params += ProtIIITDEvaluations._softParams[self.getEnumText('chooseIIITDEvaluator')]
        except:
          pass

      elif evalOrigin == DDG:
        try:
          from ddg.protocols import ProtDDGEvaluations
          params += ProtDDGEvaluations._softParams[self.getEnumText('chooseDDGEvaluator')]
        except:
          pass
      return params

    def mergeEpsLinks(self, eps, links):
      '''Concatenates the epitopes and liners into a single string'''
      res = [0] * (len(eps) + len(links))
      res[::2], res[1::2] = eps, links
      return ''.join(res)

    def getLinkers(self, nEps):
      '''Returns an iterable with the linkers given the selected strategy and number of epitopes'''
      nLinkers = nEps - 1
      outLinks = [self.inLinker.get().strip()] * nLinkers

      # todo: linkers other than random
      # todo: linkers condified in individual so no different every evaluation, for now 1 linker
      # linkDic = self.parseLinkers()
      # linkers, linkWs = list(linkDic.keys()), list(linkDic.values())
      # outLinks = numpy.random.choice(linkers, nLinkers, True, linkWs)
      return outLinks

    def parseLinkers(self):
      '''Parse the linkers summary and
      returns a dictionary: {LinkerSeq: LinkerWeight}'''
      linkDic = {}
      for linkLine in self.linkerSummary.get().split('\n'):
        if linkLine.strip():
          sDic = json.loads(')'.join(linkLine.split(')')[1:]).strip())
          linkDic.update({sDic['Linker']: sDic['Weight']})
      return linkDic

    def getEpitopeScores(self):
      '''Return a dictionary of the form {ROIId: SamplingScore} for all the input ROIs working as epitopes
      with the defined sampling scores'''
      sampDic = {}
      wDic = self.getSamplingScores()
      for roi in self.inputROIs.get():
        roiId, roiScore = roi.getObjId(), 0
        for scoreKey, scoreW in wDic.items():
          roiScore += scoreW * getattr(roi, scoreKey).get()
        sampDic[roiId] = roiScore
      return sampDic

    def getSamplingScores(self):
      weights = {}
      for sline in self.scoreSummary.get().split('\n'):
        if sline.strip():
          sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
          weights.update({sDic['Score']: sDic['Weight']})
      return weights

    def getEvalWeights(self):
      evalWs = {}
      for i, sline in enumerate(self.evalSummary.get().split('\n')):
        if sline.strip():
          sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
          source, weight = sDic.pop('Evaluation'), sDic.pop('Weight')
          if source == IEDB:
            evalKey = f'{source}_{i+1}'
          elif source == IIITD:
            evalKey = f'{source}-{sDic.pop("chooseIIITDEvaluator")}_{i + 1}'
          elif source == DDG:
            evalKey = f'{source}-{sDic.pop("chooseDDGEvaluator")}_{i + 1}'
          evalWs[evalKey] = weight
      return evalWs