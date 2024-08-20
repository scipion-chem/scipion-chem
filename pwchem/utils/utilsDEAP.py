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
import re, random, copy, time, math, warnings
import numpy as np

from deap import base, creator, tools, algorithms

from grape.grape import Individual, Grammar, mapper_eager, mapper_lazy

# GRAPE UTILS

NOREP, WEIGHT = 'noRepeated', 'weighted'
PLUS, COMMA = 'plus', 'comma'

class RepIndividual(Individual):
    '''Individual including the option to not repeated nodes'''
    def __init__(self, genome, grammar, max_depth, codon_consumption, minEps=None, maxEps=None, cProb=0.9):
        self.genome = genome
        if codon_consumption == 'lazy':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_lazy(genome, grammar, max_depth)
        elif codon_consumption == 'eager':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_eager(genome, grammar, max_depth)
        elif codon_consumption == NOREP:
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = generalMapper(genome, grammar, max_depth, minEps, maxEps, mapType=NOREP, compProb=cProb)
        elif codon_consumption == WEIGHT:
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = generalMapper(genome, grammar, max_depth, minEps, maxEps, mapType=WEIGHT, compProb=cProb)
        else:
            raise ValueError("Unknown mapper")

class RepGrammar(Grammar):
    def __init__(self, file_address, noRepNTs, weightDic):
        super().__init__(file_address)
        self.noRepNTs = noRepNTs
        self.weightDic = weightDic

def random_initialisation(ind_class, pop_size, bnf_grammar, minEps, maxEps, cProb,
                          min_init_genome_length, max_init_genome_length,
                          max_init_depth, codon_size, codon_consumption,
                          genome_representation):
    """

    """
    population = []
    for i in range(pop_size):
      init_genome_length = random.randint(min_init_genome_length, max_init_genome_length)
      genome = [random.randint(0, codon_size) for i in range(init_genome_length)]
      ind = ind_class(genome, bnf_grammar, max_init_depth, codon_consumption, minEps, maxEps, cProb)
      while not ind.phenotype:
        genome = [random.randint(0, codon_size) for i in range(init_genome_length)]
        ind = ind_class(genome, bnf_grammar, max_init_depth, codon_consumption, minEps, maxEps, cProb)
      population.append(ind)

    if genome_representation == 'list':
      return population
    elif genome_representation == 'numpy':
      for ind in population:
        ind.genome = np.array(ind.genome)
      return population
    else:
      raise ValueError("Unkonwn genome representation")

def getNonRepeatedIndex(curIdx, repIdxs, nRules):
    '''Checks whether the index for that production rule has already been used and choses next available index
    - curIdx (int): index currently chosen
    - repIdxs (list): list with already used indexes
    - nRules (int): possible index range

    '''
    if len(repIdxs) >= nRules: #if not enough rules, invalid individual
        return None
    while curIdx in repIdxs:
        curIdx += 1
        if curIdx >= nRules:
            curIdx = 0
    return curIdx

def normalizeWeights(pWs):
  '''Return weights from 0 to 1 summing 1 in total'''
  if sum(pWs) == 0:
    pWs = [1] * len(pWs)
  else:
    if min(pWs) < 0:
      pWs = [w - min(pWs) + 0.1 for w in pWs]
  pWs = np.array(pWs) / sum(pWs)
  return pWs

def getWeightedIndex(genomeSeed, repIdxs, nRules, weights):
    if len(repIdxs) >= nRules: #if not enough rules, invalid individual
        return None
    pIdxs, pWs = [], []
    for i in range(nRules):
        if i not in repIdxs:
            pIdxs.append(i)
            if weights:
                pWs.append(weights[i])
            else:
                pWs.append(1)

    pWs = normalizeWeights(pWs)
    np.random.seed(genomeSeed)
    return np.random.choice(pIdxs, 1, p=pWs)[0]


def checkPhenotype(phen):
  '''Check for unmatched production rules from the phenotype that have not been transformed to terminals'''
  unPhenots = re.findall(r" ?\<\w+\> ?", phen)
  if len(unPhenots) > 0:
    return None
  return phen

def getNEpisPerProt(protsDic, nMin, nMax):
  '''Return the number of epitopes that each protein will contribute to the multiepitope so the epitope number limits
  (nMin, nMax) are fulfilled. The contribution for each protein is chosen randomly.
  - protsDic: {protName: possibleEpitopeNumber} contains the info for the number of epitopes each protein have
  :return: {protName: nEps} with the number of epitopes to contribute
  '''
  protNames = list(protsDic.keys())
  curNEps, nEps = 0, random.randint(nMin, nMax)
  epsDic = {protName: 0 for protName in protNames}
  while curNEps < nEps:
    curPName = random.choice(protNames)
    if epsDic[curPName] < protsDic[curPName]:
      epsDic[curPName] += 1
      curNEps += 1
  return epsDic

def getNCEpisPerProt(protsDic, compProb):
  '''Return the number of compulsory epitopes that each protein will contribute to the multiepitope being compProb the
  probability for each epitope to appear.
  The contribution for each protein is chosen randomly.
  - protsDic: {protName: possibleEpitopeNumber} contains the info for the number of epitopes each protein have
  :return: {protName: nEps} with the number of epitopes to contribute
  '''
  epsDic = {}
  protNames = list(protsDic.keys())
  for pName in protNames:
    epsDic[pName] = sum([1 for i in range(protsDic[pName]) if random.random() < compProb])

  return epsDic

def getProtsDic(grammar, comp=False):
  '''Returns a dictionary with the number of possible epitopes per protein described in a grammar
  :param grammar: grammar object
  :return: {pName: nEpitopes}
  '''
  protDic = {}
  noRepNTs = grammar.noRepNTs
  for protNT in noRepNTs:
    if protNT != '<prot>' and protNT.endswith('Comp>') == comp:
      prodName = protNT[1:-3] if not comp else protNT[1:-7] + 'C'
      nEps = grammar.n_rules[grammar.non_terminals.index(protNT)]
      protDic[prodName] = nEps
  return protDic

def generalMapper(genome, grammar, max_depth, minEps, maxEps, mapType=WEIGHT, compProb=0.9):
  idxGenome, idx_depth, nodes, structure = 0, 0, 0, []
  epsCPerProt = getNCEpisPerProt(getProtsDic(grammar, comp=True), compProb=compProb)

  # Min/Maxs number of sampling epitopes are the numbers not filled by compulsory
  nComp = sum([nc for nc in epsCPerProt.values()])
  minEps, maxEps = max(0, minEps - nComp), max(0, maxEps - nComp)
  epsPerProt = getNEpisPerProt(getProtsDic(grammar), minEps, maxEps)

  curEpsPerProt, curCEpsPerProt = {pName: 1 for pName in epsPerProt}, {pName: 1 for pName in epsCPerProt}

  phenotype = grammar.start_rule
  next_NT = re.search(r"\<(\w+)\>", phenotype).group()
  n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>", phenotype)])
  list_depth = [1] * n_starting_NTs  # it keeps the depth of each branch

  # Used indexes are stored to not be repeated for non-repeatable elements
  noRepNTs = grammar.noRepNTs
  repDic = {NT_label: [] for NT_label in noRepNTs} if noRepNTs else {}
  while next_NT and idxGenome < len(genome):
    NT_index = grammar.non_terminals.index(next_NT)
    repIdxs = repDic[next_NT] if next_NT in repDic else []

    prodName, index_production_chosen = next_NT[1:-1], None
    if prodName in epsPerProt:
      # Grow eps tree until the decided number of eps per prot is met
      index_production_chosen = 1 if curEpsPerProt[prodName] < epsPerProt[prodName] else 0
      curEpsPerProt[prodName] += 1

    elif prodName in epsCPerProt:
      # Grow eps tree until the decided number of comp eps per prot is met
      index_production_chosen = 1 if curCEpsPerProt[prodName] < epsCPerProt[prodName] else 0
      curCEpsPerProt[prodName] += 1

    else:
      if mapType == WEIGHT:
        weights = grammar.weightDic[next_NT] if next_NT in grammar.weightDic else {}
        index_production_chosen = getWeightedIndex(genome[idxGenome], repIdxs, grammar.n_rules[NT_index], weights)
      elif mapType == NOREP:
        index_production_chosen = genome[idxGenome] % grammar.n_rules[NT_index]
        index_production_chosen = getNonRepeatedIndex(index_production_chosen, repIdxs, grammar.n_rules[NT_index])

    if next_NT in repDic:
        if index_production_chosen == None:
          return [None] * 7
        repDic[next_NT].append(index_production_chosen)

    structure.append(index_production_chosen)
    idxGenome += 1

    phenotype = phenotype.replace(next_NT, grammar.production_rules[NT_index][index_production_chosen][0], 1)
    list_depth[idx_depth] += 1
    if list_depth[idx_depth] > max_depth:
        break
    if grammar.production_rules[NT_index][index_production_chosen][2] == 0:  # arity 0 (T)
        idx_depth += 1
        nodes += 1
    elif grammar.production_rules[NT_index][index_production_chosen][2] == 1:  # arity 1 (PR with one NT)
        pass
    else:  # it is a PR with more than one NT
        arity = grammar.production_rules[NT_index][index_production_chosen][2]
        if idx_depth == 0:
            list_depth = [list_depth[idx_depth], ] * arity + list_depth[idx_depth + 1:]
        else:
            list_depth = list_depth[0:idx_depth] + [list_depth[idx_depth], ] * arity + list_depth[idx_depth + 1:]

    next_ = re.search(r"\<(\w+)\>", phenotype)
    if next_:
        next_NT = next_.group()
    else:
        next_NT = None

  if next_NT:
      invalid = True
      used_codons = 0
  else:
      invalid = False
      used_codons = idxGenome

  depth = max(list_depth)
  phenotype = checkPhenotype(phenotype)
  return phenotype, nodes, depth, used_codons, invalid, 0, structure


def getProteinNames(bnfGrammar):
    return [pName[0][1:-5] for pName in bnfGrammar.production_rules[1]]

def getProteinPhenotype(ind, pName):
    '''Return the substring containing the epitopes and linkers of the given protein name from the phenotype
    '''
    eps = re.findall(fr'{pName}C?\[\d+\]', ind.phenotype)
    if len(eps) > 1:
        epStr = re.findall(fr'{pName}C?\[\d+\].*{pName}C?\[\d+\]', ind.phenotype)[0]
    else:
        epStr = eps[0]
    return epStr

def performCrossProtein(child0, child1, proteinCrosses):
    '''Perform the protein epitopes swapping
    proteinCrosses = {pName: True/False}'''
    for pName, doCross in proteinCrosses.items():
        if doCross:
            phen0, phen1 = getProteinPhenotype(child0, pName), getProteinPhenotype(child1, pName)
            child0.phenotype = child0.phenotype.replace(phen0, phen1)
            child1.phenotype = child1.phenotype.replace(phen1, phen0)
    return child0, child1


def getProteinEpitopes(bnfGrammar):
    '''Returns a dictionary of the form:
    protEpDic = {protName: nEps}
    '''
    protEpDic = {}
    protNames = getProteinNames(bnfGrammar)
    for pName in protNames:
        for pRule in bnfGrammar.production_rules:
            option = pRule[0]
            if re.search(fr'{pName}\[\d+\]', option[0]):
                protEpDic[pName] = len(pRule)
                break
    return protEpDic

def getProteinLinkers(bnfGrammar):
    '''Returns a dictionary of the form:
    protLinkDic = {protName: nLinkers}
    '''
    protLinkDic = {}
    protNames = getProteinNames(bnfGrammar)
    for pName in protNames:
        for pRule in bnfGrammar.production_rules:
            option = pRule[0]
            if re.search(fr'{pName}L\[\d+\]', option[0]):
                protLinkDic[pName] = len(pRule)
                break
    return protLinkDic

def performReplaceEpitope(ind, pName, nEps):
    '''Performs the epitope mutation by replacement with another available epitope
    '''
    epIdxs = re.findall(fr'{pName}\[(\d+)\]', ind.phenotype)
    if len(epIdxs) == 0:
      return ind, True
    phen = getProteinPhenotype(ind, pName)
    mutIdx = random.choice(epIdxs)
    posibleIdx = list(set(range(nEps)) - set(epIdxs))
    if len(posibleIdx) == 0:
      return ind, True

    newIdx = random.choice(posibleIdx)
    phen = phen.replace(f'{pName}[{mutIdx}]', f'{pName}[{newIdx}]')
    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, False

def performSwapEpitopes(ind, pName):
    '''Performs the epitope mutation by swapping the order of the epitopes of a protein
    '''
    tmpEp, mutated_ = 'tmpEpXxXxX', False
    eps = re.findall(fr'{pName}C?\[\d+\]', ind.phenotype)
    if len(eps) < 2:
      return ind, True

    phen = getProteinPhenotype(ind, pName)
    ep1, ep2 = np.random.choice(eps, 2, replace=False)
    phen = phen.replace(ep1, tmpEp)
    phen = phen.replace(ep2, ep1)
    phen = phen.replace(tmpEp, ep2)

    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, False

def performReplaceLinker(ind, pName, nLinks):
    '''Performs the linker epitope mutation by replacement with another epitope linker from a protein
    '''
    linkIdxs = re.findall(fr'{pName}L\[(\d+)\]', ind.phenotype)
    phen = getProteinPhenotype(ind, pName)
    if len(linkIdxs) == 0:
      return ind, True
    mutIdx = random.choice(linkIdxs)
    newIdx = random.choice(range(nLinks))
    phen = phen.replace(f'{pName}L[{mutIdx}]', f'{pName}L[{newIdx}]')

    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, False

def mutationEpitope(ind, bnfGrammar, mutEProb=0.3, mutSwapProb=None, mutLProb=None):
    """This mutation operates at phenotype level, so no further mapping is necessary.
    The mutation will swap an epitope from 1 protein with another epitope from the same protein.
    The mutation will ensure that no repeated epitopes are generated
    - ind: Individual to be mutated
    - bnfGrammar: RepGrammar
    - mutProb: possibility of the individual to be mutated by epitope replacement
    - mutSwapProb: possibility of the individual to be mutated by epitope cluster swapping
    - mutLProb: possibility of the individual to be mutated by linker replacement
    """
    mutEProb = 0.3 if mutEProb is None else mutEProb
    mutSwapProb = mutEProb if mutSwapProb is None else mutSwapProb
    mutLProb = mutEProb if mutLProb is None else mutLProb
    
    mutated_ = False
    newInd = copy.deepcopy(ind)
    protEpDic = getProteinEpitopes(bnfGrammar)
    
    # Determining which mutation(s) to be performed. At least, one must be done
    mutRandoms = [random.random() for i in range(3)]
    doMuts = [mutR < mProb for mutR, mProb in zip(mutRandoms, [mutSwapProb, mutEProb, mutLProb])]
    if not any(doMuts):
      doMuts[random.randint(0, 2)] = True

    failMut = False
    if doMuts[0]:
      mutated_ = True
      pName = random.choice(list(protEpDic.keys()))
      newInd, failMut = performSwapEpitopes(newInd, pName)

    if doMuts[1] or failMut:
      mutated_ = True
      pName = random.choice(list(protEpDic.keys()))
      newInd, failMut = performReplaceEpitope(newInd, pName, protEpDic[pName])

    if doMuts[2] or failMut:
      mutated_ = True
      protLinkDic = getProteinLinkers(bnfGrammar)
      pName = random.choice(list(protEpDic.keys()))
      newInd, failMut = performReplaceLinker(newInd, pName, protLinkDic[pName])

    if mutated_:
        del newInd.fitness.values
    return newInd,

def crossover_multiepitope(parent0, parent1, bnfGrammar):
  '''This crossover operates at phenotype level, so no further mapping is necessary.
  The crossover will swap the entire set of epitopes from 1 protein with the same protein set of the other parent.
  - parents: Individuals to be crossed
  - bnfGrammar: RepGrammar
  '''
  child0, child1 = copy.deepcopy(parent0), copy.deepcopy(parent1)
  protNames = getProteinNames(bnfGrammar)
  proteinCrosses = {pName: False for pName in protNames}

  proteinCrosses[random.choice(protNames)] = True  # Crossing only one random protein
  child0, child1 = performCrossProtein(child0, child1, proteinCrosses)
  del child0.fitness.values, child1.fitness.values
  return child0, child1

######### ALGORITM UTILS ##############


def phenVarAnd(population, toolbox, bnfGrammar, cxpb, mutpb, mutProbs):
  """Part of an evolutionary algorithm applying only the variation part
  (crossover **and** mutation). The modified individuals have their
  fitness invalidated. The individuals are cloned so returned population is
  independent of the input population.

  :param population: A list of individuals to vary.
  :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                  operators.
  :param cxpb: The probability of mating two individuals.
  :param mutpb: The probability of mutating an individual.
  :returns: A list of varied individuals that are independent of their
            parents.

  """
  offspring = [toolbox.clone(ind) for ind in population]

  # Apply crossover and mutation on the offspring
  for i in range(1, len(offspring), 2):
    if random.random() < cxpb:
      offspring[i - 1], offspring[i] = toolbox.mate(offspring[i - 1], offspring[i], bnfGrammar)

  for i in range(len(offspring)):
    if random.random() < mutpb:
      offspring[i], = toolbox.mutate(offspring[i], bnfGrammar, *mutProbs)

  return offspring

def phenVarOr(population, toolbox, bnfGrammar, lambda_, cxpb, mutpb, mutProbs):
  """Part of an evolutionary algorithm applying only the variation part
  (crossover **or** mutation). The modified individuals have their
  fitness invalidated. The individuals are cloned so returned population is
  independent of the input population.

  :param population: A list of individuals to vary.
  :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                  operators.
  :param cxpb: The probability of mating two individuals.
  :param mutpb: The probability of mutating an individual.
  :returns: A list of varied individuals that are independent of their
            parents.

  """

  offspring = []
  while len(offspring) < lambda_:
    op_choice = random.random()
    if op_choice < cxpb:  # Apply crossover
      ind1, ind2 = [toolbox.clone(i) for i in random.sample(population, 2)]
      ind1, ind2 = toolbox.mate(ind1, ind2, bnfGrammar)
      del ind1.fitness.values
      del ind2.fitness.values
      if not ind1.invalid:
        offspring.append(ind1)
      elif not ind2.invalid:
        offspring.append(ind2)
    elif op_choice < cxpb + mutpb:  # Apply mutation
      ind = toolbox.clone(random.choice(population))
      ind, = toolbox.mutate(ind, bnfGrammar, *mutProbs)
      del ind.fitness.values
      if not ind.invalid:
        offspring.append(ind)
    else:  # Apply reproduction
      offspring.append(random.choice(population))

  return offspring

def performEvaluation(p, fitDic, toolbox):
  invalidInd = [ind for ind in p if not ind.fitness.valid]
  newInds = [ind for ind in invalidInd if ind.phenotype not in fitDic]
  oldInds = [ind for ind in invalidInd if ind.phenotype in fitDic]

  if len(newInds) > 0:
    fitnesses = toolbox.evaluate(newInds)
    for i, ind in enumerate(newInds):
      fitDic[ind.phenotype] = fitnesses[i]
      ind.fitness.values = fitnesses[i]

  for ind in oldInds:
    ind.fitness.values = fitDic[ind.phenotype]

  return p, fitDic

def initPop(toolbox, bnfGrammar, nPop, minEps, maxEps, cProb):
  return toolbox.populationCreator(pop_size=nPop, bnf_grammar=bnfGrammar, minEps=minEps, maxEps=maxEps, cProb=cProb,
                                  min_init_genome_length=95, max_init_genome_length=115, max_init_depth=35,
                                  codon_size=255, codon_consumption='weighted',
                                  genome_representation='list')


def getValids(population):
  valid0 = [ind for ind in population if not ind.invalid]
  valid = [ind for ind in valid0 if not math.isnan(ind.fitness.values[0])]
  if len(valid0) != len(valid):
    warnings.warn("Warning: There are valid individuals with fitness = NaN in the population. We will avoid them.")
  return valid

def updateHOF(halloffame, population):
  if halloffame is not None:
    halloffame.update(population)
  return halloffame

def generateRecord(logbook, stats, population, valid, halloffame, ngen, nevals, verbose):
  if verbose:
    record = stats.compile(population) if stats is not None else {}
    scores = [ind.fitness.values[0] for ind in valid]
    hofScores = [ind.fitness.values[0] for ind in halloffame]
    logKwargs = {"gen": ngen, "nevals": nevals,
                 "maxScore": max(scores), "minScore": min(scores), "meanScore": sum(scores) / len(scores),
                 "hofScores": hofScores}

    logbook.record(**logKwargs, **record)
    print(logbook.stream)
  return logbook

######### ALGORITMS ##############

def ge_eaSimpleChem(population, toolbox, bnfGrammar, ngen,
                    cxpb, mutpb,
                    mutEpb=None, mutSpb=None, mutLpb=None,
                    stats=None, halloffame=None, verbose=__debug__):
  """This algorithm reproduce the simplest evolutionary algorithm as
  presented in chapter 7 of [Back2000]_, with some adaptations to run GE
  on GRAPE.
  :param population: A list of individuals.
  :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                  operators.
  :param cxpb: The probability of mating two individuals.
  :param mutpb: The probability of mutating an individual.
  :param ngen: The number of generation.
  :param elite_size: The number of best individuals to be copied to the
                  next generation.
  :params bnfGrammar, codon_size, max_tree_depth: Parameters
                  used to mapper the individuals after crossover and
                  mutation in order to check if they are valid.
  :param stats: A :class:`~deap.tools.Statistics` object that is updated
                inplace, optional.
  :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                     contain the best individuals, optional.
  :param verbose: Whether or not to log the statistics.
  :returns: The final population
  :returns: A class:`~deap.tools.Logbook` with the statistics of the
            evolution
  """
  logbook = tools.Logbook()
  logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])
  fitnessDic = {}

  # Evaluate the individuals with an invalid fitness
  population, fitnessDic = performEvaluation(population, fitnessDic, toolbox)
  valid = getValids(population)

  # Update the hall of fame with the generated individuals
  halloffame = updateHOF(halloffame, valid)

  # Begin the generational process
  for gen in range(1, ngen + 1):
    # Select the next generation individuals
    offspring = toolbox.select(valid, len(population))
    # Vary the pool of individuals
    mutProbs = [mutEpb, mutSpb, mutLpb]
    offspring = phenVarAnd(offspring, toolbox, bnfGrammar, cxpb, mutpb, mutProbs)

    # Evaluate the individuals with an invalid fitness
    invalidInd = [ind for ind in offspring if not ind.fitness.valid]
    offspring, fitnessDic = performEvaluation(offspring, fitnessDic, toolbox)

    # Update population for next generation
    population[:] = offspring
    # Include in the population the elitist individuals
    valid = getValids(population)

    # Update the hall of fame with the generated individuals
    halloffame = updateHOF(halloffame, valid)

    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=len(invalidInd), **record)
    if verbose:
      print(logbook.stream)

  return population, logbook


def ge_eaMuLambdaChem(gaType, population, toolbox, bnfGrammar, minEps, maxEps, cProb,
                      ngen, mu, lambda_, nMigr,
                      cxpb, mutpb,
                      mutEpb=None, mutSpb=None, mutLpb=None,
                      stats=None, halloffame=None, verbose=__debug__):
  """This algorithm reproduce the Mu Lambda evolutionary algorithms from DEAP (both muPlusLambda nd MuCommaLambda)
  :param population: A list of individuals.
  :param toolbox: A :class:`~deap.base.Toolbox` that contains the evolution
                  operators.
  :param cxpb: The probability of mating two individuals.
  :param mutpb: The probability of mutating an individual.
  :param ngen: The number of generation.
  :params bnfGrammar, codon_size, max_tree_depth: Parameters
                  used to mapper the individuals after crossover and
                  mutation in order to check if they are valid.
  :param stats: A :class:`~deap.tools.Statistics` object that is updated
                inplace, optional.
  :param halloffame: A :class:`~deap.tools.HallOfFame` object that will
                     contain the best individuals, optional.
  :param verbose: Whether or not to log the statistics.
  :returns: The final population
  :returns: A class:`~deap.tools.Logbook` with the statistics of the
            evolution
  """

  logbook = tools.Logbook()
  logbook.header = ['gen', 'nevals', 'maxScore', 'minScore', 'meanScore', 'hofScores'] + (stats.fields if stats else [])
  fitnessDic = {}

  # Evaluate the individuals with an invalid fitness
  population, fitnessDic = performEvaluation(population, fitnessDic, toolbox)
  valid = getValids(population)

  # Update the hall of fame with the generated individuals
  halloffame = updateHOF(halloffame, valid)

  logbook = generateRecord(logbook, stats, population, valid, halloffame, 0, len(population), verbose)

  # Begin the generational process
  for gen in range(1, ngen + 1):
    # Select the next generation individuals
    mutProbs = [mutEpb, mutSpb, mutLpb]
    offspring = phenVarOr(population, toolbox, bnfGrammar, lambda_, cxpb, mutpb, mutProbs)
    migrants = initPop(toolbox, bnfGrammar, nMigr, minEps, maxEps, cProb) if nMigr > 0 else []

    # Evaluate the individuals with an invalid fitness
    invalidInd = [ind for ind in offspring if not ind.fitness.valid]
    offspring, fitnessDic = performEvaluation(offspring + migrants, fitnessDic, toolbox)

    # Update population for next generation
    toSelect = offspring if gaType == COMMA else population + offspring
    population[:] = toolbox.select(toSelect, mu)
    valid = getValids(population)

    # Update the hall of fame with the generated individuals
    halloffame = updateHOF(halloffame, valid)

    logbook = generateRecord(logbook, stats, population, valid, halloffame, gen, len(invalidInd), verbose)

  return population, logbook

def ge_eaMuPlusLambdaChem(*args, **kwargs):
  return ge_eaMuLambdaChem(PLUS, *args, **kwargs)

def ge_eaMuCommaLambdaChem(*args, **kwargs):
  return ge_eaMuLambdaChem(COMMA, *args, **kwargs)
