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

def myEaSimple(population, toolbox, cxpb, mutpb, ngen, stats=None, halloffame=None, verbose=__debug__):
  r"""Variation of the eaSimple which allows for batch evaluation of the individuals
  """
  logbook = tools.Logbook()
  logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

  # Evaluate the individuals with an invalid fitness
  invalid_ind = [ind for ind in population if not ind.fitness.valid]
  fitnesses = toolbox.evaluate(invalid_ind)
  for i, ind in enumerate(invalid_ind):
    ind.fitness.values = fitnesses[i]

  if halloffame is not None:
    halloffame.update(population)

  record = stats.compile(population) if stats is not None else {}
  logbook.record(gen=0, nevals=len(invalid_ind), **record)
  if verbose:
    print(logbook.stream)

  # Begin the generational process
  for gen in range(1, ngen + 1):
    offspring = toolbox.select(population, len(population))

    # Vary the pool of individuals
    offspring = algorithms.varAnd(offspring, toolbox, cxpb, mutpb)

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.evaluate(invalid_ind)
    for i, ind in enumerate(invalid_ind):
      ind.fitness.values = fitnesses[i]

    # Update the hall of fame with the generated individuals
    if halloffame is not None:
      halloffame.update(offspring)

    # Replace the current population by the offspring
    population[:] = offspring

    # Append the current generation statistics to the logbook
    record = stats.compile(population) if stats else {}
    logbook.record(gen=gen, nevals=len(invalid_ind), **record)
    if verbose:
      print(logbook.stream)

  return population, logbook

def myEaMuPlusLambda(population, toolbox, mu, lambda_, cxpb, mutpb, ngen,
                     stats=None, halloffame=None, verbose=__debug__):
  r"""Variation of the eaMuPlusLambda which allows for batch evaluation of the individuals
  """
  logbook = tools.Logbook()
  logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

  # Evaluate the individuals with an invalid fitness
  invalid_ind = [ind for ind in population if not ind.fitness.valid]
  fitnesses = toolbox.evaluate(invalid_ind)
  for i, ind in enumerate(invalid_ind):
    ind.fitness.values = fitnesses[i]

  if halloffame is not None:
    halloffame.update(population)

  record = stats.compile(population) if stats is not None else {}
  logbook.record(gen=0, nevals=len(invalid_ind), **record)
  if verbose:
    print(logbook.stream)

  # Begin the generational process
  for gen in range(1, ngen + 1):
    # Vary the population
    offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.evaluate(invalid_ind)
    for i, ind in enumerate(invalid_ind):
      ind.fitness.values = fitnesses[i]

    # Update the hall of fame with the generated individuals
    if halloffame is not None:
      halloffame.update(offspring)

    # Select the next generation population
    population[:] = toolbox.select(population + offspring, mu)

    # Update the statistics with the new population
    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=len(invalid_ind), **record)
    if verbose:
      print(logbook.stream)
    print(f'Generation ({gen}/{ngen}) done')
    print(population)

  return population, logbook

def myEaMuCommaLambda(population, toolbox, mu, lambda_, cxpb, mutpb, ngen,
                      stats=None, halloffame=None, verbose=__debug__):
  r"""Variation of the eaMuPlusLambda which allows for batch evaluation of the individuals
  """
  assert lambda_ >= mu, "lambda must be greater or equal to mu."

  logbook = tools.Logbook()
  logbook.header = ['gen', 'nevals'] + (stats.fields if stats else [])

  # Evaluate the individuals with an invalid fitness
  invalid_ind = [ind for ind in population if not ind.fitness.valid]
  fitnesses = toolbox.evaluate(invalid_ind)
  for i, ind in enumerate(invalid_ind):
    ind.fitness.values = fitnesses[i]

  if halloffame is not None:
    halloffame.update(population)

  record = stats.compile(population) if stats is not None else {}
  logbook.record(gen=0, nevals=len(invalid_ind), **record)
  if verbose:
    print(logbook.stream)

  # Begin the generational process
  for gen in range(1, ngen + 1):
    # Vary the population
    offspring = algorithms.varOr(population, toolbox, lambda_, cxpb, mutpb)

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.evaluate(invalid_ind)
    for i, ind in enumerate(invalid_ind):
      ind.fitness.values = fitnesses[i]

    # Update the hall of fame with the generated individuals
    if halloffame is not None:
      halloffame.update(offspring)

    # Select the next generation population
    population[:] = toolbox.select(offspring, mu)

    # Update the statistics with the new population
    record = stats.compile(population) if stats is not None else {}
    logbook.record(gen=gen, nevals=len(invalid_ind), **record)
    if verbose:
      print(logbook.stream)
    print(f'Generation ({gen}/{ngen}) done')

  return population, logbook


# GRAPE UTILS

class RepIndividual(Individual):
    '''Individual including the option to not repeated nodes'''
    def __init__(self, genome, grammar, max_depth, codon_consumption):
        self.genome = genome
        if codon_consumption == 'lazy':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_lazy(genome, grammar, max_depth)
        elif codon_consumption == 'eager':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_eager(genome, grammar, max_depth)
        elif codon_consumption == 'noRepeated':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_noRepeated(genome, grammar, max_depth)
        elif codon_consumption == 'weighted':
            self.phenotype, self.nodes, self.depth, \
            self.used_codons, self.invalid, self.n_wraps, \
            self.structure = mapper_weighted(genome, grammar, max_depth)
        else:
            raise ValueError("Unknown mapper")

class RepGrammar(Grammar):
    def __init__(self, file_address, noRepNTs, weightDic):
        super().__init__(file_address)
        self.noRepNTs = noRepNTs
        self.weightDic = weightDic

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

def getWeightedIndex(genomeSeed, repIdxs, nRules, weights):
    if len(repIdxs) >= nRules: #if not enough rules, invalid individual
        return None
    pIdxs, pWs = [], []
    for i in range(nRules):
        if not i in repIdxs:
            pIdxs.append(i)
            if weights:
                pWs.append(weights[i])
            else:
                pWs.append(1)

    pWs = np.array(pWs) / sum(pWs)
    np.random.seed(genomeSeed)
    return np.random.choice(pIdxs, 1, p=pWs)[0]

def mapper_noRepeated(genome, grammar, max_depth):
    """
    This mapper is similar to the previous one, but it does not consume codons
    when mapping a production rule with a single option."""

    idx_genome = 0
    phenotype = grammar.start_rule
    next_NT = re.search(r"\<(\w+)\>", phenotype).group()
    n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>", phenotype)])
    list_depth = [1] * n_starting_NTs  # it keeps the depth of each branch
    idx_depth = 0
    nodes = 0
    structure = []

    # Used indexes are stored to not be repeated for non-repeatable elements
    noRepNTs = grammar.noRepNTs
    repDic = {NT_label: [] for NT_label in noRepNTs} if noRepNTs else {}
    while next_NT and idx_genome < len(genome):
        NT_index = grammar.non_terminals.index(next_NT)
        # todo: switch from genome->indexProduction to weighted productions (use genes as random seeds?)
        index_production_chosen = genome[idx_genome] % grammar.n_rules[NT_index]
        if next_NT in repDic:
            index_production_chosen = getNonRepeatedIndex(index_production_chosen, repDic[next_NT], grammar.n_rules[NT_index])
            if index_production_chosen == None:
                break
            repDic[next_NT].append(index_production_chosen)

        structure.append(index_production_chosen)
        idx_genome += 1

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
        used_codons = idx_genome

    depth = max(list_depth)

    return phenotype, nodes, depth, used_codons, invalid, 0, structure

def mapper_weighted(genome, grammar, max_depth):
    """
    This mapper uses the genome as random seeds to generate the phenotype, having the epitopes weights that are
    used to randomly choose them"""

    idx_genome = 0
    phenotype = grammar.start_rule
    next_NT = re.search(r"\<(\w+)\>", phenotype).group()
    n_starting_NTs = len([term for term in re.findall(r"\<(\w+)\>", phenotype)])
    list_depth = [1] * n_starting_NTs  # it keeps the depth of each branch
    idx_depth = 0
    nodes = 0
    structure = []

    # Used indexes are stored to not be repeated for non-repeatable elements
    noRepNTs = grammar.noRepNTs
    repDic = {NT_label: [] for NT_label in noRepNTs} if noRepNTs else {}
    while next_NT and idx_genome < len(genome):
        NT_index = grammar.non_terminals.index(next_NT)
        repIdxs = repDic[next_NT] if next_NT in repDic else []
        weights = grammar.weightDic[next_NT] if next_NT in grammar.weightDic else {}
        index_production_chosen = getWeightedIndex(genome[idx_genome], repIdxs, grammar.n_rules[NT_index], weights)
        if next_NT in repDic:
            if index_production_chosen == None:
                break
            repDic[next_NT].append(index_production_chosen)

        structure.append(index_production_chosen)
        idx_genome += 1

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
        used_codons = idx_genome

    depth = max(list_depth)

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

def crossover_multiepitope(parent0, parent1, bnfGrammar, crossProtein=0):
    '''This crossover operates at phenotype level, so no further mapping is necessary.
    The crossover will swap the entire set of epitopes from 1 protein with the same protein set of the other parent.
    - parents: Individuals to be crossed
    - bnf_grammar: RepGrammar
    - crossProtein: individual possibility per protein to be swapped
    # todo: maybe a lower level crossover for epitope subsets
    '''
    child0, child1 = copy.deepcopy(parent0), copy.deepcopy(parent1)
    if crossProtein > 0:
        protNames = getProteinNames(bnfGrammar)
        proteinCrosses = {pName: np.random.choice([True, False], 1, p=[crossProtein, 1-crossProtein])[0]
                          for pName in protNames}
        child0, child1 = performCrossProtein(child0, child1, proteinCrosses)
        del child0.fitness.values, child1.fitness.values
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

def performReplaceEpitope(ind, pName, nEps, mutProb):
    '''Performs the epitope mutation
    '''
    mutated_ = False
    epIdxs = re.findall(fr'{pName}\[(\d+)\]', ind.phenotype)
    phen = getProteinPhenotype(ind, pName)
    for i in range(len(epIdxs)):
        if mutProb > random.random():
            mutated_ = True
            posibleIdx = list(set(range(nEps)) - set(epIdxs))
            newIdx = random.choice(posibleIdx)
            phen = phen.replace(f'{pName}[{epIdxs[i]}]', f'{pName}[{newIdx}]')
            epIdxs[i] = newIdx

    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, mutated_

def performSwapEpitopes(ind, pName, mutProb):
    '''Performs the epitope mutation
    '''
    tmpEp, mutated_ = 'tmpEpXxXxX', False
    eps = re.findall(fr'{pName}C?\[\d+\]', ind.phenotype)
    phen = getProteinPhenotype(ind, pName)
    for ep in eps:
        if mutProb > random.random():
            posEps = list(set(eps) - {ep})
            if posEps:
                mutated_ = True
                otherEp = random.choice(posEps)
                phen = phen.replace(ep, tmpEp)
                phen = phen.replace(otherEp, ep)
                phen = phen.replace(tmpEp, otherEp)

    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, mutated_

def performReplaceLinker(ind, pName, nLinks, mutProb):
    '''Performs the epitope mutation
    '''
    mutated_ = False
    linkIdxs = re.findall(fr'{pName}L\[(\d+)\]', ind.phenotype)
    phen = getProteinPhenotype(ind, pName)
    for lIdx in range(len(linkIdxs)):
        if mutProb < random.random():
            mutated_ = True
            newIdx = random.choice(range(nLinks))
            phen = phen.replace(f'{pName}L[{lIdx}]', f'{pName}L[{newIdx}]')

    ind.phenotype = ind.phenotype.replace(getProteinPhenotype(ind, pName), phen)
    return ind, mutated_

def mutationEpitope(ind, bnfGrammar, mutProb, mutSwapProb, mutLProb=None):
    """This mutation operates at phenotype level, so no further mapping is necessary.
    The mutation will swap an epitope from 1 protein with another epitope from the same protein.
    The mutation will ensure that no repeated epitopes are generated
    - ind: Individual to be mutated
    - bnfGrammar: RepGrammar
    - mutProb: individual possibility per epitope to be mutated by replacement
    - mutSwapProb: individual possibility per epitope to be mutated by swapping
    """
    mutLProb = mutProb if mutLProb is None else mutLProb
    mutated_ = False
    newInd = copy.deepcopy(ind)
    protEpDic = getProteinEpitopes(bnfGrammar)
    for pName, nEps in protEpDic.items():
        newInd, mutated_ = performReplaceEpitope(newInd, pName, nEps, mutProb)

    for pName, nEps in protEpDic.items():
        newInd, mutated_ = performSwapEpitopes(newInd, pName, mutSwapProb)

    protLinkDic = getProteinLinkers(bnfGrammar)
    for pName, nLinks in protLinkDic.items():
        newInd, mutated_ = performReplaceLinker(newInd, pName, nLinks, mutLProb)

    if mutated_:
        del newInd.fitness.values
    return newInd,

######### ALGORITMS ##############


def phenVarAnd(population, toolbox, bnfGrammar, cxpb, mutpb, mutSpb, mutLpb):
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
    offspring[i - 1], offspring[i] = toolbox.mate(offspring[i - 1], offspring[i],
                                                  bnfGrammar, cxpb)

  for i in range(len(offspring)):
    offspring[i], = toolbox.mutate(offspring[i], bnfGrammar, mutpb, mutSpb, mutLpb)

  return offspring


def ge_eaSimpleWithElitismNoTrain(population, toolbox, cxpb, mutpb, mutSpb, ngen, elite_size,
                                  bnf_grammar, mutLpb=None,
                                  report_items=None, stats=None, halloffame=None, verbose=__debug__):
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
  :params bnf_grammar, codon_size, max_tree_depth: Parameters
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

  if halloffame is None:
    if elite_size != 0:
      raise ValueError("You should add a hof object to use elitism.")
    else:
      warnings.warn('You will not register results of the best individual while not using a hof object.', UserWarning)
      logbook.header = ['gen', 'invalid'] + (stats.fields if stats else []) + ['avg_length', 'avg_nodes', 'avg_depth',
                                                                               'avg_used_codons',
                                                                               'behavioural_diversity',
                                                                               'structural_diversity',
                                                                               'fitness_diversity', 'selection_time',
                                                                               'generation_time']
  else:
    if halloffame.maxsize < 1:
      raise ValueError("HALLOFFAME_SIZE should be greater or equal to 1")
    if elite_size > halloffame.maxsize:
      raise ValueError("HALLOFFAME_SIZE should be greater or equal to ELITE_SIZE")
    logbook.header = ['gen', 'invalid'] + (stats.fields if stats else []) + ['best_ind_length', 'avg_length',
                                                                               'best_ind_nodes', 'avg_nodes',
                                                                               'best_ind_depth', 'avg_depth',
                                                                               'avg_used_codons',
                                                                               'best_ind_used_codons',
                                                                               'behavioural_diversity',
                                                                               'structural_diversity',
                                                                               'fitness_diversity', 'selection_time',
                                                                               'generation_time']

  start_gen = time.time()
  # Evaluate the individuals with an invalid fitness
  invalid_ind = [ind for ind in population if not ind.fitness.valid]
  fitnesses = toolbox.evaluate(invalid_ind)
  for i, ind in enumerate(invalid_ind):
    ind.fitness.values = fitnesses[i]

  valid0 = [ind for ind in population if not ind.invalid]
  valid = [ind for ind in valid0 if not math.isnan(ind.fitness.values[0])]
  if len(valid0) != len(valid):
    warnings.warn("Warning: There are valid individuals with fitness = NaN in the population. We will avoid them.")
  invalid = len(population) - len(
    valid0)  # We use the original number of invalids in this case, because we just want to count the completely mapped individuals

  list_structures = []
  if 'fitness_diversity' in report_items:
    list_fitnesses = []
  if 'behavioural_diversity' in report_items:
    behaviours = np.zeros([len(valid), len(valid[0].fitness_each_sample)], dtype=float)

  # for ind in offspring:
  for idx, ind in enumerate(valid):
    list_structures.append(str(ind.structure))
    if 'fitness_diversity' in report_items:
      list_fitnesses.append(str(ind.fitness.values[0]))
    if 'behavioural_diversity' in report_items:
      behaviours[idx, :] = ind.fitness_each_sample

  unique_structures = np.unique(list_structures, return_counts=False)
  if 'fitness_diversity' in report_items:
    unique_fitnesses = np.unique(list_fitnesses, return_counts=False)
  if 'behavioural_diversity' in report_items:
    unique_behaviours = np.unique(behaviours, axis=0)

  structural_diversity = len(unique_structures) / len(population)
  behavioural_diversity = len(unique_behaviours) / len(population) if 'behavioural_diversity' in report_items else 0

  # Update the hall of fame with the generated individuals
  if halloffame is not None:
    halloffame.update(valid)
    best_ind_length = len(halloffame.items[0].genome)
    best_ind_nodes = halloffame.items[0].nodes
    best_ind_depth = halloffame.items[0].depth
    best_ind_used_codons = halloffame.items[0].used_codons
    if not verbose:
      print("gen =", 0, ", Best fitness =", halloffame.items[0].fitness.values)

  length = [len(ind.genome) for ind in valid]
  avg_length = sum(length) / len(length)

  nodes = [ind.nodes for ind in valid]
  avg_nodes = sum(nodes) / len(nodes)

  depth = [ind.depth for ind in valid]
  avg_depth = sum(depth) / len(depth)

  used_codons = [ind.used_codons for ind in valid]
  avg_used_codons = sum(used_codons) / len(used_codons)

  end_gen = time.time()
  generation_time = end_gen - start_gen

  selection_time = 0

  record = stats.compile(population) if stats else {}
  logbook.record(gen=0, invalid=invalid, **record,
                   best_ind_length=best_ind_length, avg_length=avg_length,
                   best_ind_nodes=best_ind_nodes,
                   avg_nodes=avg_nodes,
                   best_ind_depth=best_ind_depth,
                   avg_depth=avg_depth,
                   avg_used_codons=avg_used_codons,
                   best_ind_used_codons=best_ind_used_codons,
                   behavioural_diversity=behavioural_diversity,
                   structural_diversity=structural_diversity,
                   selection_time=selection_time,
                   generation_time=generation_time)
  if verbose:
    print(logbook.stream)

  # Begin the generational process
  for gen in range(logbook.select("gen")[-1] + 1, ngen + 1):
    start_gen = time.time()

    # Select the next generation individuals
    start = time.time()
    offspring = toolbox.select(valid, len(population) - elite_size)
    end = time.time()
    selection_time = end - start
    # Vary the pool of individuals
    offspring = phenVarAnd(offspring, toolbox, bnf_grammar, cxpb, mutpb, mutSpb, mutLpb)

    # Evaluate the individuals with an invalid fitness
    invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
    fitnesses = toolbox.evaluate(invalid_ind)
    print('Fitnesses: ', fitnesses)
    for i, ind in enumerate(invalid_ind):
      ind.fitness.values = fitnesses[i]

    # Update population for next generation
    population[:] = offspring
    # Include in the population the elitist individuals
    for i in range(elite_size):
      population.append(halloffame.items[i])

    valid0 = [ind for ind in population if not ind.invalid]
    valid = [ind for ind in valid0 if not math.isnan(ind.fitness.values[0])]
    if len(valid0) != len(valid):
      warnings.warn(
        "Warning: There are valid individuals with fitness = NaN in the population. We will avoid in the statistics.")
    invalid = len(population) - len(
      valid0)  # We use the original number of invalids in this case, because we just want to count the completely mapped individuals

    list_structures = []
    if 'fitness_diversity' in report_items:
      list_fitnesses = []
    if 'behavioural_diversity' in report_items:
      behaviours = np.zeros([len(valid), len(valid[0].fitness_each_sample)], dtype=float)

    for idx, ind in enumerate(valid):
      list_structures.append(str(ind.structure))
      if 'fitness_diversity' in report_items:
        list_fitnesses.append(str(ind.fitness.values[0]))
      if 'behavioural_diversity' in report_items:
        behaviours[idx, :] = ind.fitness_each_sample

    unique_structures = np.unique(list_structures, return_counts=False)
    if 'fitness_diversity' in report_items:
      unique_fitnesses = np.unique(list_fitnesses, return_counts=False)
    if 'behavioural_diversity' in report_items:
      unique_behaviours = np.unique(behaviours, axis=0)

    structural_diversity = len(unique_structures) / len(population)
    behavioural_diversity = len(unique_behaviours) / len(population) if 'behavioural_diversity' in report_items else 0

    # Update the hall of fame with the generated individuals
    if halloffame is not None:
      halloffame.update(valid)
      best_ind_length = len(halloffame.items[0].genome)
      best_ind_nodes = halloffame.items[0].nodes
      best_ind_depth = halloffame.items[0].depth
      best_ind_used_codons = halloffame.items[0].used_codons
      if not verbose:
        print("gen =", gen, ", Best fitness =", halloffame.items[0].fitness.values, ", Number of invalids =", invalid)

    length = [len(ind.genome) for ind in valid]
    avg_length = sum(length) / len(length)

    nodes = [ind.nodes for ind in valid]
    avg_nodes = sum(nodes) / len(nodes)

    depth = [ind.depth for ind in valid]
    avg_depth = sum(depth) / len(depth)

    used_codons = [ind.used_codons for ind in valid]
    avg_used_codons = sum(used_codons) / len(used_codons)

    end_gen = time.time()
    generation_time = end_gen - start_gen

    # Append the current generation statistics to the logbook
    record = stats.compile(population) if stats else {}
    logbook.record(gen=gen, invalid=invalid, **record,
                   best_ind_length=best_ind_length, avg_length=avg_length,
                   best_ind_nodes=best_ind_nodes,
                   avg_nodes=avg_nodes,
                   best_ind_depth=best_ind_depth,
                   avg_depth=avg_depth,
                   avg_used_codons=avg_used_codons,
                   best_ind_used_codons=best_ind_used_codons,
                   behavioural_diversity=behavioural_diversity,
                   structural_diversity=structural_diversity,
                   selection_time=selection_time,
                   generation_time=generation_time)

    if verbose:
      print(logbook.stream)

  return population, logbook



