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

import os
import tempfile
import subprocess

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import AlignmentFile


class AlignmentViewer(Viewer):
    """Open a SAM/BAM/CRAM alignment file in IGV."""

    _label = 'Alignment IGV viewer'
    _targets = [AlignmentFile]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **kwargs):

        alignmentFile = obj.getFileName()
        alignmentFormat = (obj.getFormat() or '').upper()

        indexFile = (
            obj.getIndexFile()
            if obj.hasIndexFile()
            else None
        )

        fastaFile = (
            obj.getReferenceFasta()
            if obj.hasReferenceFasta()
            else None
        )

        gtfFile = (
            obj.getReferenceGtf()
            if obj.hasReferenceGtf()
            else None
        )

        self._checkFile(
            alignmentFile,
            '{} alignment'.format(
                alignmentFormat or 'Sequencing'
            )
        )

        if alignmentFormat in ('BAM', 'CRAM'):
            self._checkFile(
                indexFile,
                '{} index'.format(alignmentFormat)
            )

        if fastaFile:
            self._checkFile(
                fastaFile,
                'Reference FASTA'
            )

        if gtfFile:
            self._checkFile(
                gtfFile,
                'Reference GTF'
            )

        batchFile = self._createIgvBatchFile(
            alignmentFile=alignmentFile,
            fastaFile=fastaFile,
            gtfFile=gtfFile
        )

        command = (
            Plugin.getEnvActivationCommand(RNASEQ_DIC) +
            ' && igv -b "{}"'.format(batchFile)
        )

        subprocess.Popen(
            command,
            shell=True
        )

        return []

    def _createIgvBatchFile(
            self,
            alignmentFile,
            fastaFile=None,
            gtfFile=None):
        """Create an IGV batch script."""

        region = self._getFirstMappedRegion(
            alignmentFile
        )

        lines = ['new']

        if fastaFile:
            lines.append(
                'genome "{}"'.format(
                    os.path.abspath(fastaFile)
                )
            )

        if gtfFile:
            lines.append(
                'load "{}"'.format(
                    os.path.abspath(gtfFile)
                )
            )

        lines.append(
            'load "{}"'.format(
                os.path.abspath(alignmentFile)
            )
        )

        if region:
            lines.append(
                'goto {}'.format(region)
            )

        lines.append(
            'sort position'
        )

        lines.append(
            'snapshotDirectory .'
        )

        lines.append(
            'collapse'
        )

        fd, batchFile = tempfile.mkstemp(
            prefix='scipion_igv_',
            suffix='.batch',
            text=True
        )

        with os.fdopen(fd, 'w') as outputFile:
            outputFile.write(
                '\n'.join(lines) + '\n'
            )

        return batchFile

    def _getFirstMappedRegion(
            self,
            alignmentFile):
        """Return a region around the first mapped read."""

        command = (
            Plugin.getEnvActivationCommand(RNASEQ_DIC) +
            ' && samtools view -@ 2 -F 4 "{}" | head -n 1'
            .format(alignmentFile)
        )

        try:
            output = subprocess.check_output(
                command,
                shell=True,
                text=True,
                stderr=subprocess.DEVNULL
            ).strip()

        except (
                subprocess.CalledProcessError,
                OSError):
            return None

        if not output:
            return None

        fields = output.split('\t')

        if len(fields) < 4:
            return None

        chromosome = fields[2]

        try:
            position = int(fields[3])

        except ValueError:
            return None

        start = max(1, position - 500)
        end = position + 1500

        return '{}:{}-{}'.format(
            chromosome,
            start,
            end
        )

    @staticmethod
    def _checkFile(
            filePath,
            label):
        """Check that a file exists."""

        if not filePath:
            raise RuntimeError(
                '{} file was not defined.'.format(label)
            )

        if not os.path.exists(filePath):
            raise RuntimeError(
                '{} file does not exist:\n{}'
                .format(label, filePath)
            )