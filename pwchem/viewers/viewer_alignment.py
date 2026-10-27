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
import shlex
import threading
import subprocess

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import AlignmentFile


class AlignmentViewer(Viewer):
    """Open a SAM/BAM/CRAM alignment file in IGV.

    The batch script uses the stored reference, annotation and explicit index
    paths, including indexes kept outside the alignment directory. A bounded
    samtools query selects an initial locus when possible; failure to obtain
    that locus does not prevent viewing. Source files are not modified.
    """

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
            gtfFile=gtfFile,
            indexFile=indexFile
        )

        command = (
            Plugin.getEnvActivationCommand(RNASEQ_DIC) +
            ' && exec igv -b {}'.format(shlex.quote(batchFile))
        )

        try:
            process = subprocess.Popen(command, shell=True, env=Plugin.getEnviron())
        except OSError:
            os.remove(batchFile)
            raise
        threading.Thread(target=self._removeBatchAfterExit, args=(process, batchFile), daemon=True).start()
        return []

    @staticmethod
    def _removeBatchAfterExit(process, batchFile):
        process.wait()
        try:
            os.remove(batchFile)
        except FileNotFoundError:
            pass

    @staticmethod
    def _batchPath(path):
        path = os.path.abspath(path)
        if any(character in path for character in ('"', '\r', '\n')):
            raise RuntimeError('IGV batch paths cannot contain double quotes or line breaks: {}'.format(path))
        return '"{}"'.format(path)

    def _createIgvBatchFile(
            self,
            alignmentFile,
            fastaFile=None,
            gtfFile=None,
            indexFile=None):
        """Create an IGV batch script."""

        region = self._getFirstMappedRegion(
            alignmentFile, fastaFile
        )

        lines = ['new']

        if fastaFile:
            lines.append(
                'genome {}'.format(self._batchPath(fastaFile))
            )

        if gtfFile:
            lines.append(
                'load {}'.format(self._batchPath(gtfFile))
            )

        loadCommand = 'load ' + self._batchPath(alignmentFile)
        if indexFile:
            loadCommand += ' index=' + self._batchPath(indexFile)
        lines.append(loadCommand)

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
            alignmentFile, fastaFile=None):
        """Return a region around the first mapped read; tolerate unavailable samtools."""

        arguments = ['samtools', 'view', '-@', '2', '-F', '4']
        if fastaFile:
            arguments.extend(['-T', fastaFile])
        arguments.append(alignmentFile)
        command = (Plugin.getEnvActivationCommand(RNASEQ_DIC) + ' && ' +
                   ' '.join(shlex.quote(value) for value in arguments) + ' | head -n 1')

        try:
            output = subprocess.check_output(
                command,
                shell=True,
                text=True,
                stderr=subprocess.DEVNULL,
                env=Plugin.getEnviron(),
                timeout=10
            ).strip()

        except (
                subprocess.CalledProcessError,
                subprocess.TimeoutExpired,
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

        if chromosome == '*' or position < 1:
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

        if not os.path.isfile(filePath) or os.path.getsize(filePath) == 0:
            raise RuntimeError(
                '{} file is missing or empty:\n{}'
                .format(label, filePath)
            )
