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
This protocol is used to cluster a set of molecules and extract a representatives for each cluster

"""
import os

from pyworkflow.utils import Message
from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC
from pwchem import Plugin
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_filter import ProtocolBaseLibraryToSetOfMols

FINGERPRINTS = ['Morgan', 'RDKit', 'AtomPair', 'TopologicalTorsion', 'MACCS', 'Pattern', 'Layered', 'Avalon']
CLUSTERING = ['Butina', 'DBSCAN', 'HDBSCAN', 'KMedoids', 'Birch', 'BitBirch', 'LeaderPick', 'MaxMinPick']
DISTANCES = ['Tanimoto', 'Dice', 'Cosine']

scriptName = 'cluster_molecules.py'
RDKIT_FINGER, ATOMPAIR_FINGER = 'finger == 1', 'finger == 2'


class ProtClusterMolecules(ProtocolBaseLibraryToSetOfMols):
  """
  Cluster set of small molecules
  """
  _label = 'Cluster molecules'
  _OUTNAME = 'outputAtomStruct'
  _possibleOutputs = {_OUTNAME: SetOfSmallMolecules}

  # -------------------------- DEFINE param functions ----------------------
  def _defineParams(self, form):
    """ """
    form.addSection(label=Message.LABEL_INPUT)
    group = form.addGroup('Input')
    group = self.addInputParams(group)
    group.addParam('outputOnlyReps', params.BooleanParam, label="Output only representatives: ", default=True,
                   help='Output only the molecules representative of each cluster. If False, the input set with the '
                        'cluster index for each molecule will also be output')

    group = form.addGroup('Fingerprints')
    group.addParam('finger', params.EnumParam, label="Fingerprints to use: ", choices=FINGERPRINTS, default=0,
                   help='Fingerprints to use for defining the molecule features')
    group.addParam('fingerSize', params.IntParam, label="Fingerprint size: ", default=2048,
                   condition='finger < 4', help='Size of the fingerprint to generate')
    group.addParam('useChiralty', params.BooleanParam, label="Use chiralty: ", default=False,
                   expertLevel=params.LEVEL_ADVANCED, condition='finger in [0, 2, 3]',
                   help='Whether to use or not the chiralty of the molecules')
    group.addParam('radius', params.IntParam, label="Circular radius: ", default=2, condition='finger == 0',
                   expertLevel=params.LEVEL_ADVANCED, help='Radius of the circular Morgan fingerprint')
    line = group.addLine('Path distances: ', condition=RDKIT_FINGER,
                         help='Minimum and maximum path distance for the RDKit fingerprints')
    line.addParam('minPath', params.IntParam, label="Minimum: ", default=1)
    line.addParam('maxPath', params.IntParam, label="Maximum: ", default=7)
    group.addParam('useHs', params.BooleanParam, label="Use H: ", default=True, condition=RDKIT_FINGER,
                   expertLevel=params.LEVEL_ADVANCED,
                   help='Whether to use or not the Hydrogen info of the molecules')
    line = group.addLine('Atom pair distances: ', condition=ATOMPAIR_FINGER,
                         help='Minimum and maximum distance for the AtomPair fingerprints. Maximum not set if -1.')
    line.addParam('minDistance', params.IntParam, label="Minimum: ", default=1)
    line.addParam('maxDistance', params.IntParam, label="Maximum: ", default=-1)
    group.addParam('use2D', params.BooleanParam, label="Use 2D: ", default=True, condition=ATOMPAIR_FINGER,
                   expertLevel=params.LEVEL_ADVANCED, help='Whether to use or not the 2D info of the molecules')

    group = form.addGroup('Clustering')
    group.addParam('cluster', params.EnumParam, label="Clustering to use: ", choices=CLUSTERING, default=0,
                   help='Clustering method to use for grouping the molecules')
    group.addParam('distance', params.EnumParam, label="Distance to use: ", choices=DISTANCES, default=0,
                   help='Distance method to use in the clustering')

    line = group.addLine('Cutoffs: ', condition='cluster in [0, 1, 4, 5, 6]',
                         help='Cutoff to use in the clustering method to check if two items are clustered together.'
                              'In BitBirch, you can set a min, max and step values to optimize the threshold in order '
                              'to obtain the desired number of clusters.')
    line.addParam('cutoff', params.FloatParam, label="Cutoff: ", default=0.7)
    line.addParam('cutoffLow', params.FloatParam, label="Lowest: ", default=0.2, condition='cluster in [5]')
    line.addParam('cutoffStep', params.FloatParam, label="Step: ", default=0.05, condition='cluster in [5]')
    
    group.addParam('nClusters', params.IntParam, label="Number of clusters: ", default=10,
                   condition='cluster in [3, 4, 5, 7]',
                   help='Number of clusters (and therefore, output molecules) to generate.')

    group.addParam('minSamples', params.IntParam, label="Min. samples: ", default=5, condition='finger in [1, 2]',
                   help='The number of samples (or total weight) in a neighborhood for a point to be considered as '
                        'a core point')
    group.addParam('minClusterSize', params.IntParam, label="Min. cluster size: ", default=5,
                   condition='cluster in [1, 2]',
                   help='The minimum number of samples in a group for that group to be considered a cluster')
    group.addParam('branchingFactor', params.IntParam, label="Branching factor: ", default=50,
                   condition='cluster in [4, 5]', expertLevel=params.LEVEL_ADVANCED,
                   help='Maximum number of CF subclusters in each node. If a new samples enters such that the '
                        'number of subclusters exceed the branching_factor then that node is split into two nodes '
                        'with the subclusters redistributed in each')

  # --------------------------- STEPS functions ------------------------------
  def _insertAllSteps(self):
    # Insert processing steps
    self._insertFunctionStep('clusteringStep')
    self._insertFunctionStep('defineOutputStep')

  def clusteringStep(self):
    molFiles = self.getInputMolFiles()
    self.describeClustering(molFiles)

  def defineOutputStep(self):
    molDic, reps = {}, []
    with open(self._getPath('results.tsv')) as f:
        for line in f:
            sline = line.strip().split('\t')
            if not self.outputOnlyReps.get():
              molDic[sline[0]] = sline[1]

            if sline[2] == 'True':
              reps.append(sline[0])

    if self.useLibrary.get():
      repMols, outputLib = self.defineLibraryOutput(molDic, reps)
      if outputLib:
        self._defineOutputs(outputLibrary=outputLib)

    else:
        repMols, outMols = self.defineMolsOutput(molDic, reps)
        if not self.outputOnlyReps.get():
          self._defineOutputs(outputSmallMolecules=outMols)

    self._defineOutputs(outputRepSmallMolecules=repMols)

  # --------------------------- STEPS subfunctions ------------------------------
  def defineLibraryOutput(self, molDic, reps):
    repMols = SetOfSmallMolecules().create(prefix='representatives', outputPath=self._getPath())
    outputLib = None
    inLib = self.inputLibrary.get()
    mapDic = inLib.getLibraryMap(inverted=True, fullLine=True)

    if not self.outputOnlyReps.get():
      oLibFile = self._getPath('outputLibrary.smi')
      with open(oLibFile, 'w') as f:
        for smiName, line in mapDic.items():
          molFile = os.path.abspath(self._getTmpPath(f'{smiName}.smi'))
          if molFile in molDic:
            f.write(f'{line}\t{molDic[molFile]}\n')

      prevHeaders = inLib.getHeaders()
      outputLib = inLib.clone()
      outputLib.setFileName(oLibFile)
      outputLib.setHeaders(prevHeaders + ['Cluster'])

    for repName in reps:
      repName = getBaseName(repName)

      newSmallMol = SmallMolecule()
      newSmallMol = self.addLibAttributes(newSmallMol, mapDic[repName])
      molFile = os.path.abspath(self._getTmpPath(f'{repName}.smi'))
      newFile = self._getPath(f'{repName}.smi')

      tryMove = 0
      while tryMove < 5:
        try:
          os.rename(molFile, newFile)
          tryMove = 10
        except:
          tryMove += 1

      newSmallMol.setFileName(newFile)
      newSmallMol.guessMolName()
      repMols.append(newSmallMol)
    return repMols, outputLib

  def defineMolsOutput(self, molDic, reps):
    repMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(),
                                             prefix='representatives', copyInfo=True)
    outMols = None
    if not self.outputOnlyReps.get():
      outMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(),
                                               prefix='clusters', copyInfo=True)

    for mol in self.inputSmallMolecules.get():
      molFile = os.path.abspath(mol.getFileName())

      if not self.outputOnlyReps.get() and molFile in molDic:
          mol.cluster = pwobj.String(molDic[molFile])
          mol.guessMolName()
          outMols.append(mol)

      if molFile in reps:
          mol.guessMolName()
          repMols.append(mol)
    return repMols, outMols


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

  # --------------------------- UTILS functions -----------------------------------
  def getInputMolFiles(self):
    if self.useLibrary.get():
      inDir = os.path.abspath(self._getTmpPath())
      ligFiles = self.inputLibrary.get().splitInFiles(inDir)
    else:
      ligFiles = [os.path.abspath(mol.getFileName()) for mol in self.inputSmallMolecules.get()]
    return ligFiles

  def describeClustering(self, molFiles):
    paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
    self.writeParamsFile(paramsPath, molFiles)

    envDic = RDKIT_DIC if self.getEnumText("cluster") != 'BitBirch' else OPENBABEL_DIC
    Plugin.runScript(self, scriptName, paramsPath, env=envDic, cwd=self._getPath())

  def writeParamsFile(self, paramsFile, molFiles):
    molFiles = molsPDBQT2PDB(self, molFiles, self._getTmpPath())

    with open(paramsFile, 'w') as f:
      f.write('outputPath: results.tsv\n')
      f.write('molFiles: {}\n'.format(' '.join(molFiles)))

      f.write(self.getFingerParamsString())
      f.write(self.getClusterParamsString())

    return paramsFile

  def getFingerParamsString(self):
    s = f'finger: {self.getEnumText("finger")}\n'
    if self.finger.get() < 4:
      s += f'fingerSize: {self.fingerSize.get()}\n'
    if self.finger.get() in [0, 2, 3]:
      s += f'useChiralty: {self.useChiralty.get()}\n'

    if self.finger.get() == 0:
      s += f'radius: {self.radius.get()}\n'
    elif self.finger.get() == 1:
      s += f'minPath: {self.minPath.get()}\n'
      s += f'maxPath: {self.maxPath.get()}\n'
      s += f'useHs: {self.useHs.get()}\n'
    elif self.finger.get() == 2:
      s += f'minDistance: {self.minDistance.get()}\n'
      s += f'maxDistance: {self.maxDistance.get()}\n'
      s += f'use2D: {self.use2D.get()}\n'

    return s

  def getClusterParamsString(self):
    s = f'cluster: {self.getEnumText("cluster")}\n'
    s += f'distance: {self.getEnumText("distance")}\n'
    if self.cluster.get() in [4, 5]:
      s += f'branchingFactor: {self.branchingFactor.get()}\n'

    if self.cluster.get() in [0, 1, 4, 5, 6]:
      s += f'cutoff: {self.cutoff.get()}\n'
      if self.cluster.get() == 5:
          s += f'cutoffLow: {self.cutoffLow.get()}\n'
          s += f'cutoffStep: {self.cutoffStep.get()}\n'
    if self.cluster.get() in [3, 4, 5, 7]:
      s += f'nClusters: {self.nClusters.get()}\n'
    if self.cluster.get() in [1, 2]:
      s += f'minSamples: {self.minSamples.get()}\n'
      s += f'minClusterSize: {self.minClusterSize.get()}\n'

    return s
