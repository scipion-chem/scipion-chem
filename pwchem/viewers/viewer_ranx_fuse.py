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

from collections import Counter

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.gui.dialog import showError
import pyworkflow.protocol.params as params

from pwem.viewers import TableView

from pwchem.protocols import ProtocolRANXFuse


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

  def _getVisualizeDict(self):
    return {'displayNativePosFreq': self._viewNativePosFreq}

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
