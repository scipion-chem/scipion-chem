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
import json, os, shutil, re, time, random
import numpy as np

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.objects.data import SetOfSequences, Sequence

from pwchem.utils import runInParallel
from pwchem.utils.utilsDEAP import getProteinNames, ge_eaSimpleChem, ge_eaMuPlusLambdaChem, ge_eaMuCommaLambdaChem, \
  RepGrammar, RepIndividual, randomInitialisation, crossoverMultiepitope, mutationEpitope

from deap import base, creator, tools


REM, ADD, MOD = 0, 1, 2

REPORT_ITEMS = ['gen', 'invalid', 'avg', 'std', 'min', 'max',
                'fitness_test',
                'best_ind_length', 'avg_length',
                'best_ind_nodes', 'avg_nodes',
                'best_ind_depth', 'avg_depth',
                'avg_used_codons', 'best_ind_used_codons',
                'structural_diversity', 'fitness_diversity',
                'selection_time', 'generation_time']

LINKERSET = 'linkerAction==2'

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

    # <p1Init> ::= <p1C> | <p1C> <p1Link> <p1> | <p1> <p1Link> <p1C>
    for protName, nComp in compDic.items():
      if nComp > 0:
        prod4 = f'<{protName}C> | <{protName}C> <{protName}Link> <{protName}> | <{protName}> <{protName}Link> <{protName}C>'
      else:
        prod4 = f'<{protName}>'
      f.write(f'<{protName}Init> ::= {prod4}\n')
    f.write('\n')

    # <p1> ::= <p1Ep> | <p1Ep> <p1Link> <p1>
    for protName in protDic:
      prod4 = f'<{protName}Ep> | <{protName}Ep> <{protName}Link> <{protName}>'
      f.write(f'<{protName}> ::= {prod4}\n')
    f.write('\n')

    # <p1C> ::= <p1EpComp> | <p1EpComp> <p1Link> <p1C>
    for protName in protDic:
      prod4 = f'<{protName}EpComp> | <{protName}EpComp> <{protName}Link> <{protName}C>'
      f.write(f'<{protName}C> ::= {prod4}\n')
    f.write('\n')

    # <p1EpComp> ::= p1C[0] | p1C[1]
    for protName, nComp in compDic.items():
      if nComp > 0:
        prod5 = '|'.join([f'{protName}C[{i}]' for i in range(nComp)])
        f.write(f'<{protName}EpComp> ::= {prod5}\n')
    f.write('\n')

    # <p1Ep> ::= p1[0]|p1[1]|p1[2]|p1[3]|p1[4]
    for protName, nProts in protDic.items():
      prod6 = '|'.join([f'{protName}[{i}]' for i in range(nProts)])
      f.write(f'<{protName}Ep> ::= {prod6}\n')
    f.write('\n')

    # <p1Link> ::= p1L[0]|p1L[1]|p1L[2]
    for protName, nLinks in linkDic.items():
      prod7 = '|'.join([f'{protName}L[{i}]' for i in range(nLinks)])
      f.write(f'<{protName}Link> ::= {prod7}\n')


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

IEDB, IIITD, DDG, VAXIGN = 'IEDB', 'IIITD', 'DDG', 'VAXIGN'


class ProtOptimizeMultiEpitope(EMProtocol):
  """
  Optimize a MultiEpitope object based on some scores using Genetic Algorithms
  """
  _label = 'Optimize multiepitope'

  _agTypes = ['Simple', 'MuPlusLambda', 'MuCommaLambda']
  _crossTypes = ['Epitope cluster swap']
  _selectionTypes = ['Tournament', 'Roulette', 'Best']

  _evalParameters = {IEDB: ['mhc', 'pop'],
                     IIITD: ['chooseIIITDEvaluator'],
                     DDG: ['chooseDDGEvaluator'],
                     VAXIGN: ['organ']
                     }

  _gaAlgorithms = {SIMP: ge_eaSimpleChem, MUPLUS: ge_eaMuPlusLambdaChem,
                   MUCOMMA: ge_eaMuCommaLambdaChem}
  _cxAlgorithms = {P1: tools.cxOnePoint, P2: tools.cxTwoPoint, UNI: tools.cxUniform}
  _selAlgorithms = {TOUR: tools.selTournament, ROUL: tools.selRoulette, BEST: tools.selBest}

  # -------------------------- DEFINE param functions ----------------------
  def getEvaluationProtocols(self):
    origins = {}
    try:
      from iedb.protocols import ProtMHCPopulationCoverage
      origins[IEDB] = ProtMHCPopulationCoverage
    except:
      pass

    try:
      from immuno.protocols import ProtIIITDEvaluations
      origins[IIITD] = ProtIIITDEvaluations
    except:
      pass

    try:
      from immuno.protocols import ProtVaxignMLEpitopeEvaluation
      origins[VAXIGN] = ProtVaxignMLEpitopeEvaluation
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
    group = form.addGroup('Fixed input')
    group.addParam('inputCompROISets', params.MultiPointerParam, pointerClass='SetOfSequenceROIs',
                   label="Input compulsory epitope sets: ", allowsNull=True,
                   help='Select the Sets of sequence ROIs with compulsory epitopes that'
                        ' will always be included in the multiepitope')
    group.addParam('cProb', params.FloatParam, label='Compulsory probability: ', default=0.9,
                   help='Probability of compulsory epitopes to appear in the created individuals.')

    group = form.addGroup('Sampling input')
    group.addParam('inputROISets', params.MultiPointerParam, pointerClass='SetOfSequenceROIs',
                   label="Input epitope sets: ",
                   help='Select the Sets of sequence ROIs where to sample the epitopes that will be variable '
                        'in the multiepitopes')

    group = form.addGroup('Multiepitope definition')
    line = group.addLine('Number of epitopes: ',  help='Range of possible number of epitopes')
    line.addParam('minEp', params.IntParam, label='Min: ', default=3)
    line.addParam('maxEp', params.IntParam, label='Max: ', default=10)

    form.addParam('randomSeed', params.IntParam, label='Random seed: ', default=0,
                  expertLevel=params.LEVEL_ADVANCED,
                  help='Choose a random seed to ensure reproducibility')

    group = form.addGroup('Output')
    group.addParam('hallSize', params.IntParam, label='Number of output multiepitopes: ', default=5,
                   help='The number of best scored individuals that will be saved outside of the population')

    form.addSection(label='Epitope sampling')
    group = form.addGroup('Default epitope sampling scores')
    group.addParam('inScoreDef', params.StringParam, label='Select a default score: ', default='',
                   help='Choose a score to define the epitope sampling probability')
    group.addParam('scoreWeightDef', params.FloatParam, label='Default score weight: ', default=1,
                   help='Define a default weight for the score. The bigger the absolute value of the weight, the more '
                        'important it will be. If negative, the filter will look for smaller values.')
    group.addParam('addScoreDef', params.LabelParam, label='Add default score: ',
                   help='Add defined default score')

    group.addParam('scoreSummaryDef', params.TextParam, width=100, label='Default score summary:', default='',
                   help='Summary of the epitope sampling scores')

    group = form.addGroup('Specific epitope sampling scores')
    group.addParam('inSet', params.StringParam, label='Select a input set: ', default='',
                   help='Select the input set for the specific sampling scores')
    group.addParam('inScore', params.StringParam, label='Select a score: ', default='',
                   help='Select a score from the input set')
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
                   condition=f'{LINKERSET}',
                   help='Select the linker(s) to associate with the selected input protein')

    group.addParam('linkProtSet', params.StringParam, label='Select a input protein set: ', default='',
                   condition=f'{LINKERSET}',
                   help='Choose a input protein set to associate with the linkers')
    group.addParam('addLinker', params.LabelParam, label='Add linker: ', condition=f'{LINKERSET}',
                   help='Add defined linker(s) for protein')
    group.addParam('linkerSummary', params.TextParam, width=100, label='Linker summary:', default='',
                   condition=f'{LINKERSET}',
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
    group.addParam('gaType', params.EnumParam, label='GA type: ', default=MUPLUS, choices=self._agTypes,
                   help='Choose a strategy for the linkers to be chosen. For more information: '
                        'https://deap.readthedocs.io/en/master/api/algo.html#complete-algorithms')

    group.addParam('nReRun', params.IntParam, label='Number of runs: ', default=2,
                   help='Run the GA search n times and store the best individuals coming from each of them')
    group.addParam('nGen', params.IntParam, label='Number of generations: ', default=10,
                   help='Number of generations for the GA')

    line = group.addLine('Population size: ',
                         help='Number of individuals for each population group.\n'
                         'Parents: number of individuals selected each generation\n'
                         'Offspring: number of individuals generated from parents prior to selection\n'
                         'Migrants: number of individuals randomly generated prior to selection')
    line.addParam('nPop', params.IntParam, label='Parents: ', default=20)
    line.addParam('nOffs', params.IntParam, label='Offspring: ', default=20, condition='gaType!=0')
    line.addParam('nMigr', params.IntParam, label='Migrants: ', default=10, condition='gaType!=0')

    group = form.addGroup('Mutation', expertLevel=params.LEVEL_ADVANCED,)
    group.addParam('mutProb', params.FloatParam, label='Mutate probability: ', default=0.4,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Probability for a individual to be mutated for the next generation.')
    line = group.addLine('Mutation probabilities: ', expertLevel=params.LEVEL_ADVANCED,
                         help='Probabilities for the different types of mutation')
    line.addParam('mutEProb', params.FloatParam, label='Epitope replacement: ', default=0.33)
    line.addParam('mutSProb', params.FloatParam, label='Epitope cluster reorder: ', default=0.33)
    line.addParam('mutLProb', params.FloatParam, label='Linker replacement: ', default=0.33)

    group = form.addGroup('Crossover', expertLevel=params.LEVEL_ADVANCED,)
    group.addParam('crossType', params.EnumParam, label='Crossover type: ', default=0, choices=self._crossTypes,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Choose a strategy for the crossover. For more information: '
                        'https://deap.readthedocs.io/en/master/api/tools.html#crossover')
    group.addParam('crossProb', params.FloatParam, label='Crossover probability: ', default=0.4,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Probability for a pair of individuals to be crossed for the next generation.')
    group.addParam('crossIndProb', params.FloatParam, label='Epitope crossover probability: ', default=0.25,
                   condition='crossType==2', expertLevel=params.LEVEL_ADVANCED,
                   help='Probability for each epitope in the individual to be crossed independently')

    group = form.addGroup('Selection', expertLevel=params.LEVEL_ADVANCED,)
    group.addParam('selType', params.EnumParam, label='Selection type: ', default=0, choices=self._selectionTypes,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Choose a strategy for the selection. For more information: '
                        'https://deap.readthedocs.io/en/master/api/tools.html#selection')
    group.addParam('tournSize', params.IntParam, label='Tournament size: ', default=3, condition='selType==0',
                   expertLevel=params.LEVEL_ADVANCED,
                   help='The number of individuals participating in each tournament')



    form.addParallelSection(threads=4, mpi=1)

  # --------------------------- STEPS functions ------------------------------
  def _insertAllSteps(self):
    self._insertFunctionStep(self.optimizeMultiepitopeStep)
    self._insertFunctionStep(self.defineOutputStep)

  def optimizeMultiepitopeStep(self):
    #  Setting random seed
    if self.randomSeed.get() != 0:
      random.seed(self.randomSeed.get())
      np.random.seed(self.randomSeed.get())

    # Preparing utils and mappers
    epiDic = self.getEpitopesPerProt(self.inputROISets)
    linkDic = self.getLinkersPerProt()
    compDic = self.getCompEpitopesPerProt()
    nEpiDic, nLinkDic, nCompDic = self.makeCountDic(epiDic), self.makeCountDic(linkDic), self.makeCountDic(compDic)
    wDic = self.getEpitopeScores()
    self.createROIsMapper(epiDic, linkDic, compDic)

    # Parameters
    nPop, nGen = self.nPop.get(), self.nGen.get()
    nOffs, nMigr = self.nOffs.get(), self.nMigr.get()
    mutProb, crossProb = self.mutProb.get(), self.crossProb.get()

    # Generate grammar
    grammarFile = os.path.abspath(self._getPath('grammar.txt'))
    buildGrammarFile(grammarFile, nEpiDic, nLinkDic, nCompDic)
    noRepNTs = ['<prot>'] + [f'<{pName}Ep>' for pName in nEpiDic] + [f'<{pName}EpComp>' for pName in nEpiDic]
    bnfGrammar = RepGrammar(grammarFile, noRepNTs, self.convertWDic(wDic))

    # Define arguments
    deapToolbox = self.defineDeapToolbox(bnfGrammar)
    self.hof = tools.HallOfFame(self.hallSize.get(), checkSameInd)
    kwargs = {'cxpb': crossProb, 'mutpb': mutProb, 'ngen': nGen, 'minEps': self.minEp.get(), 'maxEps': self.maxEp.get(), 'cProb': self.cProb.get(),
              'mutEpb': self.mutEProb.get(), 'mutSpb': self.mutSProb.get(), 'mutLpb': self.mutLProb.get(),
              'halloffame': self.hof, 'verbose': True}
    if self.gaType.get() != SIMP:
      kwargs.update({'mu': nPop, 'lambda_': nOffs, 'nMigr': nMigr})
    gaAlgorithm = self._gaAlgorithms[self.gaType.get()]
    
    for _ in range(self.nReRun.get()):
      population = deapToolbox.populationCreator(popSize=nPop, bnfGrammar=bnfGrammar,
                                                 minEps=self.minEp.get(), maxEps=self.maxEp.get(), cProb=self.cProb.get(),
                                                 minInitGenomeLength=95, maxInitGenomeLength=115, maxInitDepth=35,
                                                 codonSize=255, codonConsumption='weighted',
                                                 genomeRepresentation='list')

      gaAlgorithm(population, deapToolbox, bnfGrammar, **kwargs)

    self.saveHOF()


  def defineOutputStep(self):
    self.writeOutputCodes()
    hof = self.parseHOF()

    sequenceSet = SetOfSequences().create(outputPath=self._getPath())
    for i, phenotype in enumerate(hof):
      seq = self.getMultiEpitopeSeqs([phenotype], phen=True)[0]
      seqObj = Sequence(sequence=seq, name=f'MultiEpitope_{i}', id=f'MultiEpitope_{i}')
      seqObj._multiepitopeScore = params.Float(hof[phenotype])
      sequenceSet.append(seqObj)

    self._defineOutputs(outputMultiepitopes=sequenceSet)


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
    toolbox.register("populationCreator", randomInitialisation, creator.Individual)

    toolbox.register("evaluate", self.multiEpitopeEval)

    # Tournament selection:
    toolbox.register("select", tools.selTournament, tournsize=self.tournSize.get())

    # Single-point crossover:
    toolbox.register("mate", crossoverMultiepitope)
    toolbox.decorate("mate", ensureEpitopeNumber(self.minEp.get(), self.maxEp.get(), bnfGrammar))

    # Flip-int mutation:
    toolbox.register("mutate", mutationEpitope)
    return toolbox

  def idxsFromCode(self, code):
    pName = code.split('[')[0]
    roiIdx = int(code.split('[')[1].split(']')[0])
    return pName, roiIdx

  def getMultiEpitopeSeqs(self, individuals, phen=False):
    '''Returns the set of multiepitope sequences codified into a set of DEAP individuals'''
    seqs = {}
    seqMapper = self.loadROIsMapper()
    for i, ind in enumerate(individuals):
      curSeqs = []
      phenotype = ind.phenotype if not phen else ind
      for roiCode in phenotype.split():
        pName, roiIdx = self.idxsFromCode(roiCode)
        curSeqs += [seqMapper[pName][roiIdx]]
      seqs[i] = ''.join(curSeqs)

    return seqs

  def getInputProteinNames(self):
    return list(self.getROIsPerProt(self.inputROISets).keys())

  def getSet2ProtName(self, multiPointer=None):
    '''Return a dic that maps the input sets to their respective protName'''
    if not multiPointer:
      multiPointer = self.inputROISets
    mapDic = {}
    for i, roiPointer in enumerate(multiPointer):
      mapDic[str(i)] = roiPointer.get().getSequenceObj().getId()
    return mapDic

  def getROIsPerProt(self, multiPointer):
    '''Return mapper with the input ROI sequences as:
        {protName: [rois]}'''
    mapDic = {}
    for roiPointer in multiPointer:
      roiSet = roiPointer.get()
      protName = roiSet.getSequenceObj().getId()
      if protName not in mapDic:
        mapDic[protName] = []
      for roi in roiSet:
        mapDic[protName].append(roi.clone())
    return mapDic

  def getEpitopesPerProt(self, multiPointer):
    '''Return mapper with the input ROI sequences as:
    {protName: [roiSequences]}'''
    epDic = {}
    mapDic = self.getROIsPerProt(multiPointer)
    for protName in mapDic:
      if protName not in epDic:
        epDic[protName] = []
      for roi in mapDic[protName]:
        epDic[protName].append(roi.getROISequence())

    return epDic

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

  def getEpitopeScores(self):
    '''Return a dictionary of the form {roiSetIdx: {ROIId: SamplingScore}} for all the input ROIs working as epitopes
    with the defined sampling scores'''
    sDic = self.parseSamplingScores()
    defScores = self.parseDefSamplingScores()

    inputROIsMapper = self.getROIsPerProt(self.inputROISets)
    sampDic = {protName: {} for protName in inputROIsMapper}
    for protName in inputROIsMapper:
      spScores = []
      for i, roi in enumerate(inputROIsMapper[protName]):
        roiScore = 0
        # Adding specific sampling scores
        if protName in sDic:
          for scoreKey, scoreW in sDic[protName].items():
            if hasattr(roi, scoreKey):
              roiScore += scoreW * getattr(roi, scoreKey).get()
              spScores.append(scoreKey)

        # Adding default sampling scores (excluding specific score attributes used)
        for scoreKey, scoreW in defScores.items():
          if hasattr(roi, scoreKey) and scoreKey not in spScores:
            roiScore += scoreW * getattr(roi, scoreKey).get()

        sampDic[protName][i] = roiScore
    return sampDic

  def convertWDic(self, wDic):
    d = {}
    for pName in wDic:
      pKey = f'<{pName}Ep>'
      d[pKey] = [score for score in wDic[pName].values()]
    return d

  def makeCountDic(self, d):
    return {key: len(d[key]) for key in d}

  def getROIsMapperPath(self):
    return self._getExtraPath('codeMapper.json')

  def createROIsMapper(self, epsD, linkD, compD):
    mapper = epsD
    mapper.update({f'{pName}L': linkD[pName] for pName in linkD})
    mapper.update({f'{pName}C': compD[pName] for pName in compD})

    with open(self.getROIsMapperPath(), "w") as outfile:
      json.dump(mapper, outfile)
    return mapper

  def getHOFPath(self):
    return self._getPath('outputMultiepitopeCodes.txt')

  def saveHOF(self):
    with open(self.getHOFPath(), 'w') as f:
      for mep in self.hof:
        f.write(f'{mep.phenotype}\t: {mep.fitness.values[0]}\n')

  def parseHOF(self):
    hof = {}
    with open(self.getHOFPath()) as f:
      for line in f:
        sline = line.split(':')
        hof[sline[0].strip()] = float(sline[1].strip())
    return hof

  # --------------------------- WIZARD functions -----------------------------------
  def buildScoreSumLine(self):
    sLine = None
    if self.inSet.get() and self.inScoreDef.get().strip() and str(self.scoreWeightDef.get()).strip():
      sLine = f'{{"Set": "{self.inSet.get().split("//")[0]}", "Score": "{self.inScore.get()}", ' \
              f'"Weight": {self.scoreWeight.get()}}}'
    return sLine

  def buildScoreSumLineDef(self):
    sLine = None
    if self.inScoreDef.get().strip() and str(self.scoreWeightDef.get()).strip():
      sLine = f'{{"Score": "{self.inScoreDef.get()}", "Weight": {self.scoreWeightDef.get()}}}'
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

  def getAllInputScores(self):
    scoresAttrs = []
    for inPoint in self.inputROISets:
      inSet = inPoint.get()
      item = inSet.getFirstItem()
      for key, attr in item.getAttributesToStore():
        scoresAttrs.append(key)
    scoresAttrs = list(set(scoresAttrs))
    scoresAttrs.sort()
    return scoresAttrs

  # --------------------------- EVALUATION functions -----------------------------------
  def multiEpitopeEval(self, individuals):
    '''Performs the evaluations of a set of DEAP individuals.
    To be used in batch with a modified GA algoritms that allows so.
    Returns a list of individual fitnesses: [(f1,), ...]
    '''
    resultsDic = self.performEvaluations(individuals)
    weigths = self.getEvalWeights()
    fitnesses = np.zeros(len(individuals))
    for (evalKey, softName), res in resultsDic.items():
      scores = np.array(list(map(float, res)))
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
    outDir = os.path.abspath(self._getExtraPath())

    resultsDic = {}
    for source, sDics in evalDics.items():
      if source == IEDB:
        epiDic = self.performCoverageAnalysis(individuals, sDics, jobs)

      elif source == IIITD:
        from immuno import Plugin as iiitdPlugin
        epiDic = iiitdPlugin.performEvaluations(sequences, sDics, jobs, iiitdPlugin.getBrowserData(),
                                                outDir=outDir, verbose=False)

      elif source == DDG:
        from ddg import Plugin as ddgPlugin
        epiDic = ddgPlugin.performEvaluations(sequences, sDics, jobs, ddgPlugin.getBrowserData(), verbose=False)

      elif source == VAXIGN:
        epiDic = self.performVaxignML(sequences, sDics, jobs)

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

  def performCoverageAnalysis(self, individuals, sDics, nt, norm=True):
    '''Run the Coverage analysis on a set of sequences
    - inSeqs: dict like {seqId: sequence}
    - sDics: dictionary as {evalKey: {parameterName: parameterValue}}
    - nt: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
    '''
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
          coverage = parseCoverageResults(oFile, norm)['average']['coverage']
        else:
          coverage = 0
        coverages.append(coverage)
      resDic[(evalKey, IEDB)] = coverages
    return resDic

  def addNullScores(self, vaxDic, inSeqs):
    nVaxDic = {}
    for seqId in inSeqs:
      if str(seqId) not in vaxDic:
        nVaxDic[seqId] = 0.0
      else:
        nVaxDic[seqId] = float(vaxDic[str(seqId)])
    return nVaxDic

  def performVaxignML(self, inSeqs, sDics, nt, norm=True):
    '''Run Vaxign-ML on a set of sequences
    - inSeqs: dict like {seqId: sequence}
    - sDics: dictionary as {evalKey: {parameterName: parameterValue}}
    - nt: int, number of jobs for parallelization
    Returns a dictionary of the form: {(evalKey, softwareName): [scores]}
    '''
    from immuno import Plugin as immunoPlugin
    from immuno.protocols.protocol_vaxignML_prediction import filterSequences, writeFasta, parseResults

    resDic = {}
    keyName = 'inputSequencesVax'
    filtSeqs = filterSequences(inSeqs)
    if len(filtSeqs) > 0:
      oDir = os.path.abspath(self._getExtraPath('vaxResults'))
      if not os.path.exists(oDir):
        os.mkdir(oDir)
      inFasta = writeFasta(filtSeqs, self._getExtraPath(f'{keyName}.fa'))

      for evalKey, evalDic in sDics.items():
        curDir = os.path.join(oDir, evalKey)
        if os.path.exists(curDir):
          shutil.rmtree(curDir)

        kwargs = {'i': inFasta, 'o': curDir, 't': evalDic['organ'], 'p': nt}
        immunoPlugin.runVaxignML(self, kwargs, cwd=oDir)

        oFile = os.path.join(curDir, f'{keyName}.result.tsv')
        vaxDic = parseResults(oFile)
        vaxDic = self.addNullScores(vaxDic, inSeqs)
        resDic[(evalKey, VAXIGN)] = list(vaxDic.values())
    else:
      for evalKey, evalDic in sDics.items():
        resDic[(evalKey, VAXIGN)] = [0.0] * len(inSeqs)

    return resDic


  def getEvalDics(self):
    '''Returns a dictionary for each of the sources chosen in the evaluations (IEDB, IIITD, DDG)
    This dictionary are of the form: {source: {(evalKey, softwareName): {parameterName: parameterValue}}}
    '''
    mapFuncs, evalDics = {}, {}
    for i, sline in enumerate(self.evalSummary.get().split('\n')):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        source = sDic.pop('Evaluation')
        if source not in evalDics:
          evalDics[source] = {}

        if source in [IEDB, VAXIGN]:
          evalKey, softName = f'{source}_{i + 1}', source

        elif source == IIITD:
          softName = sDic.pop("chooseIIITDEvaluator")
          evalKey = f'{source}-{softName}_{i + 1}'
          sDic.update({'software': softName})
          if source not in mapFuncs:
            from immuno.utils import mapEvalParamNames as iiitdMap
            mapFuncs[source] = iiitdMap

        elif source == DDG:
          softName = sDic.pop("chooseDDGEvaluator")
          evalKey = f'{source}-{softName}_{i + 1}'
          sDic.update({'software': softName})
          if source not in mapFuncs:
            from ddg.utils import mapEvalParamNames as ddgMap
            mapFuncs[source] = ddgMap

        else:
          continue

        evalDics[source].update({evalKey: sDic})

    for source, evalDic in evalDics.items():
      if source not in [IEDB, VAXIGN]:
        evalDics[source] = mapFuncs[source](evalDic)

    return evalDics
  
  def loadROIsMapper(self):
    with open(self.getROIsMapperPath()) as fIn:
      seqMapper = json.loads(fIn.read())
    return seqMapper
  
  def writeOutputCodes(self):
    seqMapper = self.loadROIsMapper()
    with open(self._getPath('codes2Seq.txt'), 'w') as f:
      for pName in seqMapper:
        for roiIdx, seq in enumerate(seqMapper[pName]):
          f.write(f'{pName}[{roiIdx}]\t{seq}\n')
        f.write('\n')

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
    while sum(jobs.values()) > totalJobs and all(jobs.values()) != 1:
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
        from immuno.protocols import ProtIIITDEvaluations
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
    protNamesDic = self.getSet2ProtName()
    for sline in self.scoreSummary.get().split('\n'):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        pName = protNamesDic[sDic['Set']]
        if pName not in weights:
          weights[pName] = {}
        weights[pName].update({sDic['Score']: sDic['Weight']})
    return weights

  def parseDefSamplingScores(self):
    '''Return {scoreName: weight}'''
    defScores = {}
    for sline in self.scoreSummaryDef.get().split('\n'):
      if sline.strip():
        sDic = json.loads(')'.join(sline.split(')')[1:]).strip())
        defScores.update({sDic['Score']: sDic['Weight']})
    return defScores

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
        if source in [IEDB, VAXIGN]:
          evalKey = f'{source}_{i + 1}'
        elif source == IIITD:
          evalKey = f'{source}-{sDic.pop("chooseIIITDEvaluator")}_{i + 1}'
        elif source == DDG:
          evalKey = f'{source}-{sDic.pop("chooseDDGEvaluator")}_{i + 1}'
        evalWs[evalKey] = weight
    return evalWs