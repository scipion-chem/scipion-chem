# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

from pyworkflow.protocol.params import LabelParam
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pwem.viewers import TableView

from pwchem.protocols import ProtPicard
from .viewer_alignment import AlignmentViewer


class PicardViewer(ProtocolViewer):
    """Open the final BAM in IGV and the optional MarkDuplicates report."""

    _label = 'Picard results viewer'
    _targets = [ProtPicard]
    _environments = [DESKTOP_TKINTER]

    def _defineParams(self, form):
        form.addSection(label='Picard results')
        form.addParam('showAlignment', LabelParam, label='Open final BAM in IGV')
        if self.protocol.markDuplicates.get():
            form.addParam('showMetricsFile', LabelParam, label='Show duplication metrics table',
                          help='Picard marks duplicate reads without removing them from the BAM.')

    def _getVisualizeDict(self):
        actions = {'showAlignment': self._showAlignment}
        if self.protocol.markDuplicates.get():
            actions['showMetricsFile'] = self._showMetricsFile
        return actions

    def _getOutputAlignment(self):
        output = getattr(self.protocol, 'outputAlignment', None)
        if output is None:
            raise RuntimeError('Picard has not published an output alignment yet.')
        AlignmentViewer._checkFile(output.getFileName(), 'Picard output BAM')
        return output

    def _getMetricsPath(self):
        if not self.protocol.markDuplicates.get():
            raise RuntimeError('MarkDuplicates was not selected for this protocol.')
        path = self.protocol._getMarkDuplicatesMetricsFile()
        AlignmentViewer._checkFile(path, 'Picard duplication metrics')
        return path

    def _showMetricsFile(self, paramName=None):
        try:
            metrics = self.protocol._readDuplicationMetrics(self._getMetricsPath())
        except (OSError, RuntimeError, ValueError) as error:
            return [self.errorMessage(str(error), 'Picard duplication metrics')]
        labels = {
            'LIBRARY': 'Library',
            'UNPAIRED_READS_EXAMINED': 'Unpaired reads examined',
            'READ_PAIRS_EXAMINED': 'Pairs examined',
            'SECONDARY_OR_SUPPLEMENTARY_RDS': 'Secondary/supplementary reads',
            'UNMAPPED_READS': 'Unmapped reads',
            'UNPAIRED_READ_DUPLICATES': 'Duplicate unpaired reads',
            'READ_PAIR_DUPLICATES': 'Duplicate pairs',
            'READ_PAIR_OPTICAL_DUPLICATES': 'Optical duplicate pairs',
            'PERCENT_DUPLICATION': 'Duplication (%)',
            'ESTIMATED_LIBRARY_SIZE': 'Estimated library size (molecules)',
        }
        columns = list(metrics[0])
        rows = [tuple(self._formatMetric(column, row[column]) for column in columns)
                for row in metrics]
        TableView(headerList=[labels.get(column, column) for column in columns], dataList=rows,
                  title='Picard duplication metrics',
                  mesg='One row per library. Duplicates are marked, not removed. '
                       'Optical duplicate pairs are included in duplicate pairs. '
                       'Read counts and pair counts have different units.',
                  height=min(max(len(rows), 3), 20), width=1350)
        return []

    @staticmethod
    def _formatMetric(column, value):
        if value in ('', None):
            return 'N/A'
        if column == 'PERCENT_DUPLICATION':
            return '{:.2f}'.format(100 * value)
        if column != 'LIBRARY':
            try:
                return '{:,}'.format(int(value))
            except ValueError:
                pass
        return value

    def _showAlignment(self, paramName=None):
        try:
            output = self._getOutputAlignment()
            viewer = AlignmentViewer(project=self.getProject(), parent=self.formWindow)
            return viewer._visualize(output)
        except (OSError, RuntimeError) as error:
            return [self.errorMessage(str(error), 'Picard alignment')]
