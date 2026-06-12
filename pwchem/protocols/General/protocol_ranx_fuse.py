# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:      Carlos Oscar Sorzano (coss@cnb.csic.es)
# *               Natalia del Rey
# *               Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Natl. Center of Biotechnology CSIC
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

"""
This protocol is used to merge different lists scores generating a combined ranking.
"""

import os, json

from pyworkflow.constants import BETA
import pyworkflow.protocol.params as params
from pyworkflow.utils import Message
from pyworkflow.object import Object, Float, String, Integer
from pwem.protocols import EMProtocol
from pwem.objects.data import SetOfStats

from pwchem.constants import NORM_STRATEGY, SCORE_BASED_METHODS, RANK_BASED_METHODS, RANX_DIC
from pwchem import Plugin as pwchemPlugin


class ProtocolRANXFuse(EMProtocol):
  """
  AI Generated:

    Protocol to fuse multiple Sets of Statistics into a single combined ranking using Ranx-based
    score fusion and rank aggregation methods.

    This protocol integrates multiple independent scoring sources (e.g. ΔΔG predictions,
    mutation scores, docking scores, or other SetOfStats outputs) into a unified consensus
    score and ranking.

    It supports both score-based and rank-based fusion strategies using the Ranx framework.

    Overview
    --------
    The protocol performs three main operations:

    1. Input parsing and alignment
       - Reads multiple input SetOfStats
       - Extracts ID attributes and score attributes
       - Builds per-set dictionaries of scores

    2. Score fusion (external execution)
       - Prepares a configuration file describing:
         * normalization method
         * fusion algorithm
         * input score dictionaries
       - Executes external script (ranx_fusion.py) using Ranx library

    3. Output reconstruction
       - Reads fused results file
       - Computes final ranking
       - Rebuilds Scipion SetOfStats with:
         * fused score
         * rank
         * optionally expanded per-set attributes

    Input
    -----
    inputSets:
        List of EMSet objects (typically SetOfStats) containing scored items.

    inSetID:
        Name/identifier of each input set used in configuration of fusion.

    inAttrName:
        Attribute used as unique identifier for each item (e.g. mutation ID).

    inAttrVal:
        Attribute containing the score values to be fused across sets.

    high:
        Boolean flag indicating whether higher scores are better.
        If False, scores are inverted internally for correct ranking.

    outName:
        Name of the output fused score attribute.

    inAttrs:
        Text configuration defining:
        - which sets contribute
        - which attributes are used
        - mapping between IDs and values

    Fusion configuration
    --------------------
    typeAlgorithm:
        Select fusion family:
        - Score-based methods
        - Rank-based methods

    scoreAlgorithm:
        Score-based fusion methods (e.g. median, average).

    rankAlgorithm:
        Rank-based fusion methods:
        - ISR
        - Log_ISR
        - LogN_ISR
        - RRF
        - RBC

    normStrategy:
        Normalization strategy applied before fusion (Ranx normalization options).

    sigma:
        Parameter for LogN_ISR method (controls frequency smoothing).

    kparameter:
        Parameter for RRF method (controls rank damping).

    phi:
        Parameter for RBC method (controls rank persistence).

    Workflow
    --------
    1. Parse input sets and build per-set ID → score mappings.
    2. Apply direction normalization (higher/lower is better).
    3. Generate run dictionaries for each dataset.
    4. Write parameter file for external Ranx execution.
    5. Execute ranx_fusion.py via plugin environment.
    6. Read fused output file (ranked scores).
    7. Build final ranked dictionary:
       - score
       - rank
    8. Clone original items and attach fused attributes.
    9. Produce final output SetOfStats.

    Output
    ------
    outputSet:
        EMSet containing original items enriched with:
        - fused score (outName)
        - fused rank (RanxRank)
        - optionally expanded per-set attributes

    Internal data structures
    -------------------------
    runDics:
        Dictionary mapping input set indices and attributes to score dictionaries.

    perSetAttrs:
        Expanded representation of all attributes per item across sets.

    ranked:
        Final mapping:
        ID → {score, rank}

    Validation
    ----------
    - Ensures sigma ∈ [0,1] for LogN_ISR
    - Ensures k ∈ [10,100] for RRF
    - Ensures phi ∈ [0,1] for RBC

    Summary
    -------
    This protocol enables meta-analysis of multiple scoring systems by:
    - normalizing heterogeneous score sources
    - applying state-of-the-art fusion algorithms (Ranx)
    - producing a unified ranking of biological entities

    Notes
    -----
    - Relies on external script (ranx_fusion.py) for fusion computation.
    - Designed for large-scale comparative scoring datasets.
    - Particularly useful for integrating multiple prediction models.
    - Final ranking is deterministic based on selected fusion method.
  """
  _label = 'Ranx Score Fusion'
  _devStatus = BETA

  # -------------------------- DEFINE param functions ----------------------
  def _addSetOfStatsForm(self, form):
    group = form.addGroup('Input sets')
    group.addParam('inputSets', params.MultiPointerParam, pointerClass="EMSet",
                   label='Set of Stats', allowsNull=False,
                   help='Select the SetOfStats to perform the fusion function.')
    group = form.addGroup('Define attribute ID')
    group.addParam('inSetID', params.StringParam, label="Input set name: ",
                   help='Select the set from the input sets to define its attributes to fuse')
    group.addParam('inAttrName', params.StringParam, label="Input ID attribute: ",
                   help='Select the attribute to use as ID (only1 per input Set!) for the score combination')
    group = form.addGroup('Define attribute values')
    group.addParam('inAttrVal', params.StringParam, label="Input scores attribute: ",
                   help='Select the attribute in the sets where the value of each element is found. '
                        'e.g: mutation score.\nThey are the values that will be fused')
    group.addParam('high', params.BooleanParam, label="Higher is better: ", default=True,
                   help='Determines whether the score is ascending or decresing, so if a higher value of the score is '
                        'a "better" value, in order to properly correlate inverse scores')
    group.addParam('addAttr', params.LabelParam, label='Add score: ',
                   help='Add ID-score defined')

    group = form.addGroup('Summary')
    group.addParam('outName', params.StringParam, label="Output name for the fused score: ", default='RanxScore',
                   expertLevel=params.LEVEL_ADVANCED, help='Output name for fused scores')
    group.addParam('inAttrs', params.TextParam, width=70, default='', label='List of inputs to fuse: ',
                   help='Inputs to define the attribute fusion.\n'
                        'This contains the input list of attributes that will be fused, with their respective input '
                        'sets, attribute IDs and value columns.\nThe output will combine the input sets, stacking '
                        'elements with the same attribute ID. If other attributes are shared, the ones in earliest sets'
                        ' will be saved.')

  def addFusionAlgorithmForm(self, form):
    form.addParam('typeAlgorithm', params.EnumParam, default=0,
                  choices=['Score-based Methods', 'Rank-based Methods'],
                  label="Type of fusion method", allowsNull=False,
                  help='Select the type of fusion method:\n'
                       '1. Score-based methods use the value of the ΔΔG as a basis for combining '
                       'different SetOfStats.\n'
                       '2. Rank-based methods consider the inverse of the rank of each mutation in '
                       'each of the SetOfStats to combine them.'
                       '\n\nYou can find more information at https://amenra.github.io/ranx/fusion/')

    form.addParam('scoreAlgorithm', params.EnumParam, condition='typeAlgorithm==0', default=0,
                  choices=SCORE_BASED_METHODS, label='Score-based Methods', allowsNull=False,
                  help='Select the fusion algorithm:\n'
                       '1. MED: calculates the median.\n'
                       '2. ANZ: calculates the average, ignoring lists in which a mutation does not appear.')

    form.addParam('normStrategy', params.EnumParam, condition='typeAlgorithm==0', default=4,
                  choices=NORM_STRATEGY, label="Normalization Strategy", allowsNull=False,
                  help='Select the normalization strategy you want to apply to the data.'
                       '\nMore information about each strategy can be found at '
                       'https://amenra.github.io/ranx/normalization/.')

    form.addParam('rankAlgorithm', params.EnumParam, condition='typeAlgorithm==1', default=0,
                  choices=RANK_BASED_METHODS, label='Rank-based Methods', allowsNull=False,
                  help='Select the merging algorithm:\n'
                       '1. ISR (Inverted Square Ranks): uses the inverse of the rank of each mutation as a '
                       'score, and combines it with the frequency of appearance in the lists.\n'
                       '2. Log_ISR: similar to ISR, but uses the logarithm of the frequency of appearance of '
                       'each mutation.\n'
                       '3. LogN_ISR: adjusted version of Log_ISR in which a sigma parameter is added to the '
                       'frequency of appearance, so that mutations that appear in only one SetOfStats receive'
                       'a small score, avoiding it from being zero.\n'
                       '4. RRF (Reciprocal Rank Fusion): uses the inverse of the rank of each mutation as a '
                       'score, and a parameter k that prevents mutations with high ddg values in one SetOfStats '
                       'but not in the other from getting too high a final score.\n'
                       '5. RBC (Rank-Biased Centroids): based on information on the ranks given to the mutations. '
                       "It includes a persistence parameter phi that controls the importance of the mutation's position "
                       'in each list for the combined ranking.')

    form.addParam('sigma', params.FloatParam, condition='typeAlgorithm==1 and rankAlgorithm==2',
                  label='Sigma value for LogN_ISR:', allowsNull=False,
                  help='Specify the parameter sigma that will regulate the frequency of appearance of the mutations. It '
                       'must be between 0-1.')
    form.addParam('kparameter', params.FloatParam, condition='typeAlgorithm==1 and rankAlgorithm==3',
                  label='k value for RRF:', allowsNull=False,
                  help='Specify the parameter k that will regulate the final score of the mutations. It must be'
                       'between 10-100.')
    form.addParam('phi', params.FloatParam, condition='typeAlgorithm==1 and rankAlgorithm==4',
                  label='phi value for RBC:', allowsNull=False,
                  help="Specify the persistence parameter phi that controls the importance of the mutation's position. "
                       'It must be between 0-1.')

  def _defineParams(self, form):
    """ Define the input parameters that will be used.
    Params:
        form: this is the form to be populated with sections and params.
    """
    form.addSection(label=Message.LABEL_INPUT)
    self._addSetOfStatsForm(form)

    form.addSection(label='Fusion parameters')
    group = form.addGroup('Fusion algorithms')
    self.addFusionAlgorithmForm(group)

  # --------------------------- STEPS functions ------------------------------
  def _insertAllSteps(self):
    self._insertFunctionStep('performFuseStep')
    self._insertFunctionStep('createOutputStep')

  def performFuseStep(self):
    inSets = self.getInputSets()
    inAttrDic = self.getInputAttrsDic()
    normMethod, fusionMethod, kwargs = self.getArgs()

    runDics = {}
    for key, attrVals in inAttrDic.items():
      inPointIdx, attrID = key.split('-')
      inSet = inSets[int(inPointIdx)]
      for attrVal in attrVals:
        runDics[f"{inPointIdx}-{attrVal[0]}"] = self.createAttrDic(inSet, attrID, attrVal)

    paramsFile = self.writeParamsFile(runDics, normMethod, fusionMethod, kwargs)
    pwchemPlugin.runScript(self, 'ranx_fusion.py', paramsFile, RANX_DIC)

  def createOutputStep(self):
    inSet = self.inputSets[0].get()
    outAttrName = self.getOutAttrName()
    outData = self.parseOutputFile()
    inAttrDic = self.getInputAttrsDic()
    originalAttrs = self.getOriginalAttrs(inAttrDic)
    perSetAttrs = self.buildPerSetAttributeDictionary(inAttrDic)
    outSet = inSet.createCopy(self._getPath(), copyInfo=True)
    for item in inSet:
      inID = str(item.getAttributeValue(outAttrName))
      # Clone item to avoid modifying original object
      newItem = item.clone()

      # Remove original attributes dynamically
      for attrName in originalAttrs:
        if hasattr(newItem, attrName):
          delattr(newItem, attrName)

      # Add renamed attributes in deterministic order
      if inID in perSetAttrs:

        orderedAttrs = perSetAttrs[inID]

        for attrName, value in orderedAttrs.items():
          try:
            setattr(newItem, attrName, Float(float(value)))
          except:
            setattr(newItem, attrName, String(str(value)))

      # Add fusion outputs at the end
      scoreComb = outData[inID]["score"]
      rankComb = outData[inID]["rank"]

      setattr(newItem, self.outName.get(), Float(scoreComb))
      setattr(newItem, "RanxRank", Integer(rankComb))
      outSet.append(newItem)

    self._defineOutputs(outputSet=outSet)

  #  -------------------------- MAIN FUNCTIONS -----------------------------------

  def getInputSets(self):
    return [inPoint.get() for inPoint in self.inputSets]

  def getInputAttrsDic(self):
    '''Returns: {pointerIdx-AttrID: [[AttrValue, highBetter=-1|1], ...]}.
    If lower is better, highBetter=-1 and the values will be multiplied by -1, changing the order
    If a pointer idx has several AttrID, only the first will be used as ID
    for that set.
    '''
    inAttrsDic, idxsKey = {}, {}
    for inLine in self.inAttrs.get().split('\n'):
      if inLine:
        inLineDic = json.loads(inLine.strip())
        pointIdx = int(inLineDic['Set Idx'])
        lineID = inLineDic['ID']
        key = f'{pointIdx}-{lineID}'
        if pointIdx not in idxsKey or idxsKey[pointIdx] == key:
          if key not in inAttrsDic:
            inAttrsDic[key] = []
          highMult = 1 if inLineDic["Higher_is_better"] == 'True' else -1
          inAttrsDic[key].append([inLineDic['Values'], highMult])
          idxsKey[pointIdx] = key
    return inAttrsDic

  def getOutAttrName(self):
    attrDic = self.getInputAttrsDic()
    for k in attrDic:
      if k.split('-')[0] == '0':
        return k.split('-')[1]


  def getArgs(self):
    if self.normStrategy.get():
      normMethod = NORM_STRATEGY[self.normStrategy.get()]
    else:
      normMethod = None

    if self.typeAlgorithm.get() == 0:
      fusionMethod = SCORE_BASED_METHODS[self.scoreAlgorithm.get()]
    else:
      fusionMethod = RANK_BASED_METHODS[self.rankAlgorithm.get()]

    if fusionMethod == "logn_isr":
      kwargs = {"sigma": self.sigma.get()}
    elif fusionMethod == "rrf":
      kwargs = {"k": self.kparameter.get()}
    elif fusionMethod == "rbc":
      kwargs = {"phi": self.phi.get()}
    else:
      kwargs = {}
    return normMethod, fusionMethod, kwargs

  def createAttrDic(self, inSet, attrName, attrValue):
    '''Returns a dictionary {ID: Value} with the IDs in the ID attribute and value attributes columns respectively
    '''
    mutations = {}
    for obj in inSet:
      value, mult = float(obj.getAttributeValue(attrValue[0])), attrValue[1]
      mutations[str(obj.getAttributeValue(attrName))] = value * mult

    return mutations

  def buildPerSetAttributeDictionary(self, inAttrDic):
    """ Builds a dictionary where the values for each attribute of each input set are stored.
    Returns:
      {
        itemID: {
          "attr1_setIdx_0": value,
          "attr2_setIdx_0": value,
          "attr1_setIdx_1": value,
          ...
        }
      }
    """
    inSets = self.getInputSets()
    data = {}

    for key, attrVals in inAttrDic.items():
      inPointIdx, attrID = key.split('-')
      inSet = inSets[int(inPointIdx)]
      for obj in inSet:
        itemID = str(obj.getAttributeValue(attrID))
        if itemID not in data:
          data[itemID] = {}
        for attrVal in attrVals:
          attrName = attrVal[0]
          newAttrName = f"{attrName}_setIdx_{inPointIdx}"
          value = obj.getAttributeValue(attrName)
          data[itemID][newAttrName] = value
        
    orderedData = {}
    for itemID, attrs in data.items():
      orderedAttrs = dict(sorted(attrs.items(), key=lambda x: (int(x[0].split('_setIdx_')[-1]), x[0].split('_setIdx_')[0])))
      orderedData[itemID] = orderedAttrs
    return orderedData

  def getOriginalAttrs(self, inAttrDic):
    originalAttrs = set()
    for _, attrVals in inAttrDic.items():
      for attrVal in attrVals:
        originalAttrs.add(attrVal[0])
    return originalAttrs

  def getOutFile(self):
    return self._getExtraPath("rankAggregation.tsv")

  def writeParamsFile(self, runDics, norm, fus, kwargs):
    paramsFile = self._getExtraPath('inputParams.txt')
    with open(paramsFile, 'w') as f:
      f.write(f'outputPath:: {self.getOutFile()}\n')
      f.write(f'normMethod:: {norm}\n')
      f.write(f'fusionMethod:: {fus}\n')
      f.write(f'kwargs:: {kwargs}\n')

      dicLine = f'runDics:: {runDics}\n'
      f.write(dicLine.replace("'", '"'))

    return paramsFile

  def parseOutputFile(self):
    scores = {}
    with open(self.getOutFile()) as f:
      f.readline()
      for line in f:
        _, mutName, scoreComb = line.strip().split('\t')
        scores[mutName] = float(scoreComb)

    sortedItems = sorted(scores.items(), key=lambda x: x[1], reverse=True)

    ranked = {}
    currentRank = 1
    for idx, (mut, score) in enumerate(sortedItems):
      currentRank = idx + 1
      ranked[mut] = {"score": score, "rank": currentRank}

    with open(self.getOutFile(), 'w') as f:
      f.write("Rank\tMut\tScore\n")
      for mut, data in ranked.items():
        f.write(f"{data['rank']}\t{mut}\t{data['score']}\n")
    return ranked

  # --------------------------- INFO functions -----------------------------------
  def _validate(self):
    errors = []

    if self.typeAlgorithm.get() == 1:
      if self.rankAlgorithm.get() == 2:
        if float(self.sigma.get()) > 1 or float(self.sigma.get()) < 0:
          errors.append('The sigma parameter must be a float between 0 and 1.')
      elif self.rankAlgorithm.get() == 3:
        if float(self.kparameter.get()) > 100 or float(self.kparameter.get()) < 10:
          errors.append('The k parameter must be a float between 10 and 100.')
      elif self.rankAlgorithm.get() == 4:
        if float(self.phi.get()) > 1 or float(self.phi.get()) < 0:
          errors.append('The phi parameter must be a float between 0 and 1.')

    return errors

  def _summary(self):
    summary = []
    outputFile = self._getExtraPath('RankAggregation.tsv')
    if os.path.exists(outputFile):
      with open(outputFile) as f:
        summary.append(f.read())
    return summary

  def _methods(self):
    methods = []
    methods.append("Combination of the free energy change (ΔΔG) in protein-protein interactions "
                   "due to a mutation calculated by different protocols into a single value. "
                   "A combined ranking is also generated.")
    return methods

  def _citations(self):
    return ['BassaniR22']

  # --------------------------- WIZARD functions -----------------------------------
  def createElementLine(self):
    setStr, attrName, attrValue, high = self.inSetID.get().strip(), self.inAttrName.get().strip(), \
                                        self.inAttrVal.get().strip(), self.high.get()
    setIdx = setStr.split('//')[0]

    roiDef = f'{{"Set Idx": "{setIdx}", "ID": "{attrName}", "Values": "{attrValue}", "Higher_is_better": "{high}"}}\n'
    return roiDef
