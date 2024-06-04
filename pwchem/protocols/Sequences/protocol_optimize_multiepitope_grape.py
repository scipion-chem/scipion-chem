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
import numpy, json, random, time, os, re, sys

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.utils import runInParallel
from pwchem.utils.utilsDEAP import *

from grape.grape import random_initialisation
from deap import base, creator, tools, algorithms



REM, ADD, MOD = 0, 1, 2

REPORT_ITEMS = ['gen', 'invalid', 'avg', 'std', 'min', 'max',
                'fitness_test',
                'best_ind_length', 'avg_length',
                'best_ind_nodes', 'avg_nodes',
                'best_ind_depth', 'avg_depth',
                'avg_used_codons', 'best_ind_used_codons',
                'structural_diversity', 'fitness_diversity',
                'selection_time', 'generation_time']


def buildGrammarFile(gFile, protDic, linkDic, compDic):
  '''Build the multiepitope grammar file with information about the number of proteins and epitopes'''
  with open(gFile, 'w') as f:
    # <multiEp> ::= <prot> <linkRoot> <prot> <linkRoot> <prot>
    prod1 = ' <linkRoot> '.join(["<prot>"] * len(protDic))
    f.write(f'<multiEp> ::= {prod1}\n\n')

    # <prot> ::= <p1Init> | <p2Init> | <p3Init>
    protNames = [f'<{protName}Init>' for protName in protDic]
    prod2 = ' | '.join(protNames)
    f.write(f'<prot> ::= {prod2}\n')

    # <linkRoot> ::= <p1Link> | <p2Link> | <p3Link>
    linkNames = [f'<{protName}Link>' for protName in protDic]
    prod3 = ' | '.join(linkNames)
    f.write(f'<linkRoot> ::= {prod3}\n\n')

    # <p1Init> ::= <p1EpComp> | <p1EpComp> <p1Link> <p1> | <p1> <p1Link> <p1EpComp>
    for protName, nComp in compDic.items():
      if nComp > 0:
        prod4 = f'<{protName}Comp> | <{protName}Comp> <{protName}Link> <{protName}> | <{protName}> <{protName}Link> <{protName}Comp>'
      else:
        prod4 = f'<{protName}>'
      f.write(f'<{protName}Init> ::= {prod4}\n')
    f.write('\n')

    # <p1> ::= <p1Ep> | <p1Ep> <p1Link> <p1>
    for protName in protDic:
      prod4 = f'<{protName}Ep> | <{protName}Ep> <{protName}Link> <{protName}>'
      f.write(f'<{protName}> ::= {prod4}\n')
    f.write('\n')

    # <p1EpComp> ::= p1C[0] <p1Link> p1C[1]
    for protName, nComp in compDic.items():
      if nComp > 0:
        prod5 = f' <{protName}Link> '.join([f'{protName}C[{i}]' for i in range(nComp)])
        f.write(f'<{protName}Comp> ::= {prod5}\n')
    f.write('\n')

    # <p1Ep> ::= p1[0]|p1[1]|p1[2]|p1[3]|p1[4]
    for protName, nProts in protDic.items():
      prod5 = '|'.join([f'{protName}[{i}]' for i in range(nProts)])
      f.write(f'<{protName}Ep> ::= {prod5}\n')
    f.write('\n')

    # <p1Link> ::= p1L[0]|p1L[1]|p1L[2]
    for protName, nLinks in linkDic.items():
      prod6 = '|'.join([f'{protName}L[{i}]' for i in range(nLinks)])
      f.write(f'<{protName}Link> ::= {prod6}\n')


def countEpitopes(ind, protNames):
  n = 0
  for pName in protNames:
    n += len(re.findall(fr'{pName}C?\[\d+\]', ind.phenotype))
  return n


def ensureEpitopeNumber(epMin, epMax, bnfGrammar):
  def decorator(func):
    def wrapper(*args, **kargs):
      offspring = func(*args, **kargs)
      protNames = getProteinNames(bnfGrammar)
      for child in offspring:
        nEps = countEpitopes(child, protNames)
        if nEps < epMin or nEps > epMax:
          child.invalid = True
      return offspring

    return wrapper

  return decorator


def checkSameInd(ind1, ind2):
  '''Check that two multiEpitopes are not formed by the same epitopes'''
  return ind1.phenotype == ind2.phenotype


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

IEDB, IIITD, DDG = 'IEDB', 'IIITD', 'DDG'


class ProtOptimizeMultiEpitopeGrape(EMProtocol):
  """
  Optimize a MultiEpitope object based on some scores using Genetic Algorithms
  """
  _label = 'Optimize multiepitope grape'

  _agTypes = ['Simple', 'MuPlusLambda', 'MuCommaLambda']
  # _mutTypes = ['Random', 'Shuffle']
  _crossTypes = ['Epitope cluster swap']
  _selectionTypes = ['Tournament', 'Roulette', 'Best']

  _evalParameters = {IEDB: ['mhc', 'pop'],
                     IIITD: ['chooseIIITDEvaluator'],
                     DDG: ['chooseDDGEvaluator'],
                     }

  _gaAlgorithms = {SIMP: ge_eaSimpleChem, MUPLUS: ge_eaMuPlusLambdaChem, 
                   MUCOMMA: ge_eaMuCommaLambdaChem}
  _cxAlgorithms = {P1: tools.cxOnePoint, P2: tools.cxTwoPoint, UNI: tools.cxUniform}
  _selAlgorithms = {TOUR: tools.selTournament, ROUL: tools.selRoulette, BEST: tools.selBest}

  _indMemory = {}

  # -------------------------- DEFINE param functions ----------------------
  def getEvaluationProtocols(self):
    origins = {}
    try:
      from iedb.protocols import ProtMHCIIPopulationCoverage
      origins[IEDB] = ProtMHCIIPopulationCoverage
    except:
      pass

    try:
      from iiitd.protocols import ProtIIITDEvaluations
      origins[IIITD] = ProtIIITDEvaluations
    except:
      pass

    # DDG websites closed to programmatic access
    # try:
    #   from ddg.protocols import ProtDDGEvaluations
    #   origins[DDG] = ProtDDGEvaluations
    # except:
    #   pass

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
      if evalLabel == IEDB:
        self._defineMHCEvalParams(evalGroup, allCond)
      elif hasattr(evalProtocol, '_defineEvalParams'):
        evalProtocol()._defineEvalParams(evalGroup, allCond=allCond)
    return evalGroup

  def _defineParams(self, form):
    """ """
    form.addSection(label=Message.LABEL_INPUT)
    group = form.addGroup('Input')
    group.addParam('inputROISets', params.MultiPointerParam, pointerClass='SetOfSequenceROIs',
                   label="Input epitope sets: ",
                   help='Select the Sets of sequence ROIs where to sample the epitopes')

    group = form.addGroup('Multiepitope definition')
    line = group.addLine('Number of epitopes: ',  help='Range of possible number of epitopes')
    line.addParam('minEp', params.IntParam, label='Min: ', default=3)
    line.addParam('maxEp', params.IntParam, label='Max: ', default=10)

    form.addParam('randomSeed', params.IntParam, label='Random seed: ', default=0,
                  expertLevel=params.LEVEL_ADVANCED,
                  help='Choose a random seed to ensure reproducibility')

    form.addSection(label='Epitope sampling')
    group = form.addGroup('Epitope sampling scores')
    group.addParam('inSet', params.StringParam, label='Select a input set: ', default='',
                   help='Choose a score to define the epitope sampling probability')
    group.addParam('inScore', params.StringParam, label='Select a score: ', default='',
                   help='Choose a score to define the epitope sampling probability')
    group.addParam('scoreWeight', params.FloatParam, label='Score weight: ', default=1,
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
    group.addParam('inputCompROISets', params.MultiPointerParam, pointerClass='SetOfSequenceROIs',
                   label="Input compulsory epitope sets: ", allowsNull=True,
                   help='Select the Sets of sequence ROIs with compulsory epitopes')
    # group.addParam('manualProb', params.FloatParam, label='Manual weight: ', default=1,
    #                help='Define the probability of including the selected epitope in the multiepitope.'
    #                     'If 1, the epitope will always be included')
    # group.addParam('addManual', params.LabelParam, label='Add score: ',
    #                help='Add defined epitope manual sampling')
    # group.addParam('manualSummary', params.TextParam, width=100, label='Manual summary:', default='',
    #                help='Summary of the manual sampling epitopes. They will be added randomly independently of '
    #                     'the sampling epitopes based on scores')

    form.addSection(label='Linkers')
    group = form.addGroup('Multiepitope linkers')
    group.addParam('linkerAction', params.EnumParam, label='Linker strategy: ', default=0,
                   choices=['Single', 'Random', 'Per input Set'],
                   help='Choose a strategy for the linkers to be chosen')
    group.addParam('inLinker', params.StringParam, label='Define linker sequence: ', default='',
                   condition='linkerAction==0',
                   help='Define the sequence of the linker to be used')

    group.addParam('inLinkerSet', params.PointerParam, pointerClass='SetOfSequences',
                   condition='linkerAction in [1, 2]', label='Define linker set: ', allowsNull=True,
                   help='Define the set of sequences of the linker to be chosen randomly or per protein')
    group.addParam('seleLinker', params.StringParam, label='Selected linker(s): ', default='',
                   condition='linkerAction==2',
                   help='Select the linker(s) to associate with the selected input protein')

    group.addParam('linkProtSet', params.StringParam, label='Select a input protein set: ', default='',
                   condition='linkerAction==2',
                   help='Choose a input protein set to associate with the linkers')
    group.addParam('addLinker', params.LabelParam, label='Add linker: ', condition='linkerAction==2',
                   help='Add defined linker(s) for protein')
    group.addParam('linkerSummary', params.TextParam, width=100, label='Linker summary:', default='',
                   condition='linkerAction==2',
                   help='Summary of the defined linkers.')

    form.addSection(label='Multiepitope evaluations')
    group = form.addGroup('Multiepitope evaluations')
    group.addParam('multiEval', params.EnumParam, label='Select evaluation origin: ', default=0,
                   choices=self.getEvaluationOrigins(),
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
    group.addParam('nGen', params.IntParam, label='Number of generations: ', default=10,
                   help='Number of generations for the GA')
    group.addParam('nPop', params.IntParam, label='Population size: ', default=20,
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
    line = group.addLine('Mutation probabilities: ', help='Probabilities for the different types of mutation')
    line.addParam('mutEProb', params.FloatParam, label='Epitope replacement: ', default=0.33)
    line.addParam('mutSProb', params.FloatParam, label='Epitope cluster reorder: ', default=0.33)
    line.addParam('mutLProb', params.FloatParam, label='Linker replacement: ', default=0.33)

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
    group.addParam('hallSize', params.IntParam, label='Hall of fame size: ', default=5,
                   help='The number of best scored individuals that will be saved outside of the population')



    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- STEPS functions ------------------------------
  def _insertAllSteps(self):
    self._insertFunctionStep(self.optimizeMultiepitopeStep)
    self._insertFunctionStep(self.defineOutputStep)

  def optimizeMultiepitopeStep(self):
    if self.randomSeed.get() != 0:
      random.seed(self.randomSeed.get())
      np.random.seed(self.randomSeed.get())

    epiDic = self.getEpitopesPerProt(self.inputROISets)
    linkDic = self.getLinkersPerProt()
    compDic = self.getCompEpitopesPerProt()
    nEpiDic, nLinkDic, nCompDic = self.makeCountDic(epiDic), self.makeCountDic(linkDic), self.makeCountDic(compDic)
    wDic = self.getEpitopeScores()
    self.seqsMapper = self.getROIsMapper(epiDic, linkDic, compDic)

    nPop, nGen = self.nPop.get(), self.nGen.get()
    mutProb, crossProb = self.mutProb.get(), self.crossProb.get()

    grammarFile = os.path.abspath(self._getPath('grammar.txt'))
    buildGrammarFile(grammarFile, nEpiDic, nLinkDic, nCompDic)
    nonRepNTs = ['<prot>'] + [f'<{pName}Ep>' for pName in nEpiDic]
    bnfGrammar = RepGrammar(grammarFile, nonRepNTs, self.convertWDic(wDic))

    deapToolbox = self.defineDeapToolbox(bnfGrammar)
    # todo: fix: when cleaning phenotypes, ensureall prots in: probably better to discard and redo the ind
    population = deapToolbox.populationCreator(pop_size=nPop, bnf_grammar=bnfGrammar,
                                               min_init_genome_length=95, max_init_genome_length=115, max_init_depth=35,
                                               codon_size=255, codon_consumption='weighted',
                                               genome_representation='list')
    print('Initial pop: ', [ind.phenotype for ind in population])

    hof = tools.HallOfFame(self.hallSize.get(), checkSameInd)
    kwargs = {'cxpb': crossProb, 'mutpb': mutProb, 'ngen': nGen,
              'mutEpb': self.mutEProb.get(), 'mutSpb': self.mutSProb.get(), 'mutLpb': self.mutLProb.get(),
              'halloffame': hof, 'verbose': True}
    if self.gaType.get() != SIMP:
      lamb = round(nPop * self.varProp.get())
      kwargs.update({'mu': nPop, 'lambda_': int(lamb)})

    gaAlgorithm = self._gaAlgorithms[self.gaType.get()]
    print('gaAlgorithm: ', gaAlgorithm)
    population, logbook = gaAlgorithm(population, deapToolbox, bnfGrammar, **kwargs)

    print('Final pop: ', [ind.phenotype for ind in population])
    print('Scores: ', [ind.fitness.values for ind in population])
    print('Fame: ', [ind.phenotype for ind in hof.items])
    print('Scores: ', [ind.fitness.values for ind in hof.items])

    best = hof.items[0]
    print("Best individual: \n", best.phenotype, best.fitness.values)

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
    return errors

  def _warnings(self):
    ws = []
    return ws

  # --------------------------- DEAP functions -----------------------------------
  def defineDeapToolbox(self, bnfGrammar):
    toolbox = base.Toolbox()

    # define a single objective, minimising fitness strategy:
    creator.create("FitnessMin", base.Fitness, weights=(1.0,))
    creator.create('Individual', RepIndividual, fitness=creator.FitnessMin)
    toolbox.register("populationCreator", random_initialisation, creator.Individual)

    toolbox.register("evaluate", self.multiEpitopeEval)

    # Tournament selection:
    toolbox.register("select", tools.selTournament, tournsize=self.tournSize.get())

    # Single-point crossover:
    toolbox.register("mate", crossover_multiepitope)
    toolbox.decorate("mate", ensureEpitopeNumber(self.minEp.get(), self.maxEp.get(), bnfGrammar))

    # Flip-int mutation:
    toolbox.register("mutate", mutationEpitope)
    return toolbox

  def idxsFromCode(self, code):
    pName = code.split('[')[0]
    roiIdx = int(code.split('[')[1].split(']')[0])
    return pName, roiIdx

  def getMultiEpitopeSeqs(self, individuals):
    '''Returns the set of multiepitope sequences codified into a set of DEAP individuals'''
    seqs = {}
    for i, ind in enumerate(individuals):
      curSeqs = []
      for roiCode in ind.phenotype.split():
        pName, roiIdx = self.idxsFromCode(roiCode)
        curSeqs += [self.seqsMapper[pName][roiIdx]]
      seqs[i] = ''.join(curSeqs)

    return seqs

  def getInputProteinNames(self):
    return list(self.getInputROIsMapper().keys())

  def getROIsPerProt(self, multiPointer):
    '''Return mapper with the input ROI sequences as:
        {protName: [rois]}'''
    mapDic = {}
    for roiPointer in multiPointer:
      roiSet = roiPointer.get()
      protName = roiSet.getSequenceObj().getId()
      mapDic[protName] = []
      for roi in roiSet:
        mapDic[protName].append(roi)
    return mapDic

  def getEpitopesPerProt(self, multiPointer):
    '''Return mapper with the input ROI sequences as:
    {protName: [roiSequences]}'''
    mapDic = {}
    for roiPointer in multiPointer:
      roiSet = roiPointer.get()
      protName = roiSet.getSequenceObj().getId()
      mapDic[protName] = []
      for roi in roiSet:
        mapDic[protName].append(roi.getROISequence())
    return mapDic

  def getLinkersPerProt(self):
    '''Return a dictionary with the available linkers per protein
    {pNameL: [linkers]}'''
    linkDic = {}
    protNames = self.getInputProteinNames()
    if self.linkerAction.get() == 0:
      linkDic = {pName: [self.inLinker.get()] for pName in protNames}
    elif self.linkerAction.get() == 1:
      for pName in protNames:
        linkDic[pName] = [roi.getROISequence() for roi in self.inLinkerSet.get()]
    elif self.linkerAction.get() == 2:
      linkDic = self.parseLinkerSummary()

    return linkDic

  def getCompEpitopesPerProt(self):
    '''Return a dictionary with the compulsory epitopes per protein
    {pNameC: [comp epis]}'''
    inDic = self.getEpitopesPerProt(self.inputCompROISets)

    cEpiDic = {}
    protNames = self.getInputProteinNames()
    for pName in protNames:
      if pName in inDic:
        cEpiDic[pName] = inDic[pName]
      else:
        cEpiDic[pName] = []

    return cEpiDic

  def getInputROIsMapper(self):
    '''Returns a dict as {pName: roiSet}'''
    return {roiPointer.get().getSequenceObj().getId(): roiPointer.get() for roiPointer in self.inputROISets}

  def getEpitopeScores(self):
    '''Return a dictionary of the form {roiSetIdx: {ROIId: SamplingScore}} for all the input ROIs working as epitopes
    with the defined sampling scores'''
    wDic = self.parseSamplingScores()
    sampDic = {protName: {} for protName in wDic}
    inputROIsMapper = self.getInputROIsMapper()
    for protName in wDic:
      for roi in inputROIsMapper[protName]:
        roiId, roiScore = roi.getObjId(), 0
        for scoreKey, scoreW in wDic[protName].items():
          roiScore += scoreW * getattr(roi, scoreKey).get()
        sampDic[protName][roiId] = roiScore
    return sampDic

  def convertWDic(self, wDic):
    d = {}
    for pName in wDic:
      pKey = f'<{pName}Ep>'
      d[pKey] = [score for score in wDic[pName].values()]
    return d

  def makeCountDic(self, d):
    return {key: len(d[key]) for key in d}

  def getROIsMapper(self, epsD, linkD, compD):
    mapper = epsD
    mapper.update({f'{pName}L': linkD[pName] for pName in linkD})
    mapper.update({f'{pName}C': compD[pName] for pName in compD})
    return mapper

  # --------------------------- WIZARD functions -----------------------------------
  def buildScoreSumLine(self):
    sLine = f'{{"Set": "{self.inSet.get().split("//")[0]}", "Score": "{self.inScore.get()}", ' \
            f'"Weight": {self.scoreWeight.get()}}}'
    return sLine

  def buildManualSumLine(self):
    sLine = f'{{"Epitope": "{self.inROI.get()}", "Prob": {self.manualProb.get()}}}'
    return sLine

  def getLinkerSeq(self, linkerStr):
    seq = ''
    for linkSeq in self.inLinkerSet.get():
      if linkSeq.__str__().strip() == linkerStr.strip():
        seq = linkSeq.getSequence()
    return seq

  def buildLinkerSumLine(self):
    linkSeqs = [self.getLinkerSeq(linkerStr) for linkerStr in self.seleLinker.get().split(', ')]
    linkers = ' // '.join(linkSeqs)
    sLine = f'{{"Set": "{self.linkProtSet.get().split("//")[0]}", "Linker(s)": "{linkers}" }}'
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
    print('ResultDic. ', resultsDic, len(resultsDic))
    weigths = self.getEvalWeights()
    print('ws: ', weigths, len(weigths))
    fitnesses = numpy.zeros(len(individuals))
    print('fitnesses: ', fitnesses, len(fitnesses))
    for (evalKey, softName), res in resultsDic.items():
      print('res: ', res, len(res))
      # todo: toxinpred returning only 4 scores
      scores = numpy.array(list(map(float, res)))
      print('Scores: ', scores)
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
    print('Individuals: ', individuals, len(individuals))
    sequences = self.getMultiEpitopeSeqs(individuals)
    print('Sequences: ', sequences, len(sequences))

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

  def getROIsFromInd(self, ind, prot2ROIsDic):
    '''Return the ROIs present in a ind'''
    curEps = []
    pNames = self.getInputProteinNames()
    for roiCode in ind.phenotype.split():
      pName, roiIdx = self.idxsFromCode(roiCode)
      if pName in pNames:
        curEps += [prot2ROIsDic[pName][roiIdx]]
    return curEps

  def removeEmptyAlleleROIs(self, rois, mhc):
    from iedb.utils import getMHCAlleles
    nRois = []
    for roi in rois:
      alleles = getMHCAlleles(roi, mhc)
      if alleles:
        nRois.append(roi)
    return nRois

  def performCoverageAnalysis(self, individuals, sDics, nt):
    from iedb import Plugin as iedbPlugin
    from iedb.utils import translateArea, buildMHCCoverageArgs, parseCoverageResults
    resDic = {}
    roisMapDic = self.getROIsPerProt(self.inputROISets)
    for evalKey, evalDic in sDics.items():
      oDir = os.path.abspath(self._getExtraPath())

      selPops = [p.strip() for p in evalDic['pop'].split(',')]
      selPops = translateArea(selPops)

      coveArgs = []
      for i, ind in enumerate(individuals):
        epFile = self._getExtraPath(f'inputEpitopes_{evalKey}_{i}.tsv')
        inROIs = self.getROIsFromInd(ind, roisMapDic)
        inROIs = self.removeEmptyAlleleROIs(inROIs, evalDic['mhc'])
        if inROIs:
          coveArgs += buildMHCCoverageArgs(inROIs, epFile, selPops, evalDic['mhc'], oDir, separated=False)
        else:
          oFile = os.path.join(oDir, f'inputEpitopes_{evalKey}_{i}_All_results.tsv')
          if os.path.exists(oFile):
            os.remove(oFile)

      if coveArgs:
        runInParallel(iedbPlugin.runPopulationCoverage, None,
                      paramList=[coveArg for coveArg in coveArgs], jobs=nt)

      coverages = []
      for i in range(len(individuals)):
        oFile = os.path.join(oDir, f'inputEpitopes_{evalKey}_{i}_All_results.tsv')
        if os.path.exists(oFile):
          coverage = parseCoverageResults(oFile)['average']['coverage']
        else:
          coverage = 0
        coverages.append(coverage)
      resDic[(evalKey, IEDB)] = coverages
    return resDic

  def getEvalDics(self):
    '''Returns a dictionary for each of the sources chosen in the evaluations (IEDB, IIITD, DDG)
    This dictionary are of the form: {source: {(evalKey, softwareName): {parameterName: parameterValue}}}
    '''
    mapFuncs, evalDics = {}, {}
    for i, sline in enumerate(self.evalSummary.get().split('\n')):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        source, weight = sDic.pop('Evaluation'), sDic.pop('Weight')
        if not source in evalDics:
          evalDics[source] = {}

        if source == IEDB:
          evalKey, softName = f'{source}_{i + 1}', source

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
      if source != IEDB:
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
      jobs[source] = len(sDics)

    i = 0
    sources = list(jobs.keys())
    # Decrease the number of threads per source until the totaljobs is met or all 1
    while sum(jobs.values()) > totalJobs and not all(jobs.values()) == 1:
      if jobs[sources[i % len(sources)]] > 1:
        jobs[sources[i % len(sources)]] -= 1
      i += 1

    i = 0
    while sum(jobs.values()) < totalJobs:
      jobs[sources[i % len(sources)]] += 1
      i += 1

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

  def parseSamplingScores(self):
    '''Return {protName: {scoreName: weight}}'''
    weights = {}
    protNames = self.getInputProteinNames()
    for sline in self.scoreSummary.get().split('\n'):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        pName = protNames[int(sDic['Set'])]
        if not pName in weights:
          weights[pName] = {}
        weights[pName].update({sDic['Score']: sDic['Weight']})
    return weights

  def parseLinkerSummary(self):
    '''Return a dictionary {setId: [linkers]}'''
    linkDic = {}
    protNames = self.getInputProteinNames()
    for sline in self.linkerSummary.get().split('\n'):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        pName = protNames[int(sDic['Set'])]
        if pName not in linkDic:
          linkDic[pName] = []
        linkers = sDic['Linker(s)'].split(' // ')
        linkDic[pName] += linkers
    return linkDic

  def getEvalWeights(self):
    evalWs = {}
    for i, sline in enumerate(self.evalSummary.get().split('\n')):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        source, weight = sDic.pop('Evaluation'), sDic.pop('Weight')
        if source == IEDB:
          evalKey = f'{source}_{i + 1}'
        elif source == IIITD:
          evalKey = f'{source}-{sDic.pop("chooseIIITDEvaluator")}_{i + 1}'
        elif source == DDG:
          evalKey = f'{source}-{sDic.pop("chooseDDGEvaluator")}_{i + 1}'
        evalWs[evalKey] = weight
    return evalWs