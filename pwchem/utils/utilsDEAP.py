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

from deap import base, creator, tools, algorithms

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