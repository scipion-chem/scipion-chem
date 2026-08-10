# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Judith Maestro Ciria
# *
# * Biocomputing Unit, CNB-CSIC
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
Viewer for the RANX Fuse protocol. Shows how often each "NativePosition" appears among the
best-ranked mutations, to spot the native positions most favorable to change.
"""

import re
import operator
import statistics
from collections import Counter, defaultdict

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.gui.dialog import showError
import pyworkflow.protocol.params as params

from pwem.viewers import TableView

from pwchem.protocols import ProtocolRANXFuse

PER_SET_ATTR_PATTERN = re.compile(r'^(.+)_setIdx_(\d+)$')

CONSENSUS_OPS = ['<', '<=', '>', '>=']
CONSENSUS_OP_FUNCS = {'<': operator.lt, '<=': operator.le, '>': operator.gt, '>=': operator.ge}


class RANXFuseViewer(ProtocolViewer):
  """ Viewer for the RANX Fuse protocol: NativePosition frequency among the top ranked mutations. """
  _label = 'Viewer RANX Fuse'
  _targets = [ProtocolRANXFuse]
  _environments = [DESKTOP_TKINTER]

  def _defineParams(self, form):
    form.addSection(label='Native position frequency')
    form.addParam('topN', params.IntParam, default=20,
                  label='Number of top RanxRank mutations to consider: ',
                  help='Only the mutations with the best (lowest) "RanxRank", up to this number, are used '
                       'to count how often each "NativePosition" appears.\n'
                       'Requires the output set to have the "NativePosition" attribute, i.e. the protocol '
                       'must have been run with "Extract native position from mutation ID" enabled.')
    form.addParam('displayNativePosFreq', params.LabelParam,
                  label='Show NativePosition frequency table')

    form.addSection(label='Predictor score distribution')
    self._statAttrChoices = self.getAvailableStatAttrs()
    form.addParam('statAttr', params.EnumParam, choices=self._statAttrChoices, default=0,
                  label='Attribute to summarize: ',
                  help='Choose which of the per-predictor score attributes fused by RANX Fuse '
                       '(e.g. "zscore", "ddg"...) to summarize below.')
    form.addParam('displayPredictorStats', params.LabelParam,
                  label='Show predictor score distribution table',
                  help='For each input SetOfStats fused, shows the '
                       'distribution of the attribute selected above across all the mutations: '
                       'mean, standard deviation and percentage of positive/negative values.\n'
                       'This table is independent from the output set: it summarizes each '
                       'predictor as a whole, it does not add columns to the per-mutation output '
                       'table.')

    form.addParam('consensusOperator', params.EnumParam, choices=CONSENSUS_OPS, default=0,
                  label='Comparison: ',
                  help='Comparison applied, for each predictor, between its value of the '
                       'attribute selected above and the threshold below.\n')
    form.addParam('consensusThreshold', params.FloatParam, default=0.0,
                  label='Threshold: ',
                  help='A predictor "agrees" on a mutation when "value <comparison> threshold" '
                       'holds for it, where <comparison> is chosen above.')
    form.addParam('displayConsensusMutations', params.LabelParam,
                  label='Show values in all predictors',
                  help='Lists the mutations whose value of the attribute selected above (e.g. '
                       '"zscore" or "ddg") satisfies the chosen comparison/threshold in EVERY '
                       'fused input set (predictor) at once, i.e. mutations agreed on by all '
                       'predictors, not just one.\n')

  def getAvailableStatAttrs(self):
    try:
      names = set()
      for item in self.protocol.outputSet:
        for attrName, attrObj in item.getAttributes():
          match = PER_SET_ATTR_PATTERN.match(attrName)
          if match and attrObj.get() is not None:
            names.add(match.group(1))
      if names:
        return sorted(names, key=lambda n: (n.lower() != 'zscore', n.lower()))
    except Exception:
      pass
    return ['zscore']

  def _getVisualizeDict(self):
    return {'displayNativePosFreq': self._viewNativePosFreq,
            'displayPredictorStats': self._viewPredictorStats,
            'displayConsensusMutations': self._viewConsensusMutations}

  def _viewNativePosFreq(self, paramName=None):
    outSet = self.protocol.outputSet
    topN = self.topN.get()

    firstItem = outSet.getFirstItem()
    if not firstItem.hasAttribute('NativePosition'):
      showError('Missing attribute',
               'The output set does not have a "NativePosition" attribute.\nRe-run the protocol with '
               '"Extract native position from mutation ID" enabled.', self.getTkRoot())
      return

    records = [(item.getAttributeValue('RanxRank'), item.getAttributeValue('NativePosition'))
               for item in outSet]
    records.sort(key=lambda r: r[0])
    topRecords = records[:topN]

    counts = Counter(nativePos for _, nativePos in topRecords)
    dataList = counts.most_common()

    if not dataList:
      showError('No data', 'No data available to display.', self.getTkRoot())
      return

    TableView(headerList=['NativePosition', 'Frequency'],
             dataList=dataList,
             mesg='Frequency of each native position among the %d best ranked mutations (RanxRank)'
                  % len(topRecords),
             title='NativePosition frequency (top %d)' % topN,
             height=min(len(dataList), 25), width=350)

  def getSelectedAttr(self):
    choices = getattr(self, '_statAttrChoices', None) or self.getAvailableStatAttrs()
    return choices[self.statAttr.get()]

  def getPerItemValues(self, outSet, attrUsed):
    idAttrName = self.protocol.getOutAttrName()
    perItem = {}
    for item in outSet:
      vals = {}
      for attrName, attrObj in item.getAttributes():
        match = PER_SET_ATTR_PATTERN.match(attrName)
        if match and match.group(1) == attrUsed and attrObj.get() is not None:
          vals[int(match.group(2))] = float(attrObj.get())
      if vals:
        perItem[str(item.getAttributeValue(idAttrName))] = vals
    return perItem

  def _viewPredictorStats(self, paramName=None):
    outSet = self.protocol.outputSet
    attrUsed = self.getSelectedAttr()

    valuesBySet = defaultdict(list)
    for vals in self.getPerItemValues(outSet, attrUsed).values():
      for setIdx, value in vals.items():
        valuesBySet[setIdx].append(value)

    if not valuesBySet:
      showError('No data', 'No "%s_setIdx_<n>" attribute was found in the output set.' % attrUsed,
                self.getTkRoot())
      return

    dataList = []
    for setIdx in sorted(valuesBySet.keys()):
      values = valuesBySet[setIdx]
      n = len(values)
      mean = statistics.mean(values)
      std = statistics.pstdev(values) if n > 1 else 0.0
      pctPos = 100 * sum(1 for v in values if v > 0) / n
      pctNeg = 100 * sum(1 for v in values if v < 0) / n
      dataList.append((self.getSetLabel(setIdx), n, round(mean, 3), round(std, 3),
                       round(pctPos, 1), round(pctNeg, 1)))

    TableView(headerList=['Predictor', 'N mutations', 'Mean %s' % attrUsed,
                         'Std %s' % attrUsed, '% positive', '% negative'],
             dataList=dataList,
             mesg='Distribution of the "%s" attribute per input predictor, across all its mutations.'
                  % attrUsed,
             title='Predictor score distribution',
             height=min(len(dataList), 25), width=550)

  def _viewConsensusMutations(self, paramName=None):
    outSet = self.protocol.outputSet
    attrUsed = self.getSelectedAttr()
    threshold = self.consensusThreshold.get()
    opStr = CONSENSUS_OPS[self.consensusOperator.get()]
    opFunc = CONSENSUS_OP_FUNCS[opStr]

    perItem = self.getPerItemValues(outSet, attrUsed)
    if not perItem:
      showError('No data', 'No "%s_setIdx_<n>" attribute was found in the output set.' % attrUsed,
                self.getTkRoot())
      return

    allSetIdxs = set()
    for vals in perItem.values():
      allSetIdxs.update(vals.keys())
    setIdxs = sorted(allSetIdxs)

    rows = []
    for mutID, vals in perItem.items():
      if set(vals.keys()) == allSetIdxs and all(opFunc(vals[i], threshold) for i in setIdxs):
        rows.append((mutID,) + tuple(round(vals[i], 3) for i in setIdxs))
    rows.sort(key=lambda r: sum(r[1:]) / len(setIdxs), reverse=opStr in ('>', '>='))

    if not rows:
      showError('No matching values',
                'No mutation has "%s" %s %.3g in all %d predictors at once.'
                % (attrUsed, opStr, threshold, len(setIdxs)), self.getTkRoot())
      return

    headerList = ['Mutation'] + [self.getSetLabel(i) for i in setIdxs]
    TableView(headerList=headerList, dataList=rows,
             mesg='Mutations with "%s" %s %.3g in all %d fused predictors at once (%d found), '
                  'sorted by average value.'
                  % (attrUsed, opStr, threshold, len(setIdxs), len(rows)),
             title='Values in all predictors',
             height=min(len(rows), 30), width=min(1200, 220 + 160 * len(setIdxs)))

  def getSetLabel(self, setIdx):
    try:
      inSet = self.protocol.getInputSets()[setIdx]
    except Exception:
      return 'Set %d' % setIdx

    try:
      protocol = self.protocol.getProject().getObject(inSet.getObjParentId())
      if protocol is not None and protocol.getRunName():
        return protocol.getRunName()
    except Exception:
      pass

    try:
      if inSet is not None and inSet.getObjLabel():
        return inSet.getObjLabel()
    except Exception:
      pass
    return 'Set %d' % setIdx
