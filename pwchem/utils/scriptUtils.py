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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

'''Utils used in scripts. Imports containing Bio (and other modules) cannot be performed, so these functions
must be placed separatedly. These are also available in scipion-chem.'''

import threading

def divideMultiPDB(file):
  '''Divides a PDB into the different models (start with MODEL finish with ENDMDL) and return the string
  with the text for each of them'''
  pdbStrs = []

  with open(file) as fIn:
    fileText = fIn.read()
    endSplit = 'ENDMDL' if 'ENDMDL' in fileText else 'END'
    pdbs = fileText.split(f'\n{endSplit}\n')
    if len(pdbs) > 1:
      pdbs = pdbs[:-1]

  names = []
  for i, pdb in enumerate(pdbs):
    molName = "molecule_%s" % (i + 1)
    for line in pdb.strip().split('\n'):
      if line.startswith('COMPND'):
        preName = '_'.join(line.split()[1:])
        if not preName in names:
          molName = preName
          names.append(molName)
          break

    pdbStrs.append(pdb)
  return pdbStrs

def makeSubsets(oriSet, nt, cloneItem=True):
  '''Returns a list of subsets, given a set and the number of subsets'''
  if nt > len(oriSet):
      nt = len(oriSet)

  subsets = [[] for _ in range(nt)]
  for i, obj in enumerate(oriSet):
    nObj = obj.clone() if cloneItem else obj
    subsets[i%nt].append(nObj)
  return subsets

def performBatchThreading(task, inSet, nt, cloneItem=True, **kwargs):
  '''Uses threading to divide a task over an input set among a number of threads
  Input:
      -task: function. Task to perform
      -inSet: iterable. Will be divided and inputed for the task
      -nt: int. Number of threads to use

      -cloneItem: whether to perform "clone" function over the objects in inSet (necessary for Scipion objects)
  '''
  threads, outLists = [], [[] for i in range(nt)]
  subsets = makeSubsets(inSet, nt, cloneItem=cloneItem)

  for it, curSet in enumerate(subsets):
    t = threading.Thread(target=task, args=(curSet, outLists, it), kwargs=kwargs, daemon=False)
    threads.append(t)
    t.start()

  for t in threads:
    t.join()

  return [item for sublist in outLists for item in sublist]