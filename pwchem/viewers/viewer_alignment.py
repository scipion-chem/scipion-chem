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
import hashlib
import fcntl
import logging
import queue
import signal
import tkinter as tk
from tkinter import ttk

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
        job = {'cancel': threading.Event(), 'messages': queue.Queue()}
        window = tk.Toplevel(self.getTkRoot())
        window.title('Preparing IGV')
        status = ttk.Label(window, text='Preparing the alignment viewer...', padding=15)
        status.pack(fill='x')
        progress = ttk.Progressbar(window, mode='indeterminate', length=360)
        progress.pack(padx=15, pady=5)
        progress.start()
        ttk.Button(window, text='Cancel', command=window.destroy).pack(pady=10)
        window.bind('<Destroy>', lambda event: job['cancel'].set() if event.widget == window else None)
        worker = threading.Thread(target=self._prepareAndLaunch,
                                  args=(job, alignmentFile, fastaFile, gtfFile, indexFile), daemon=True)
        worker.start()
        window.after(100, self._pollPreparation, window, status, job)
        return []

    def _pollPreparation(self, window, status, job):
        """Handle UI updates only on Tk's main thread."""
        if job['cancel'].is_set():
            return
        while True:
            try:
                kind, message = job['messages'].get_nowait()
            except queue.Empty:
                break
            if kind == 'status':
                status.configure(text=message)
            else:
                window.destroy()
                if kind == 'error':
                    self.errorMessage(message, 'IGV viewer').show()
                return
        window.after(100, self._pollPreparation, window, status, job)

    def _prepareAndLaunch(self, job, alignmentFile, fastaFile, gtfFile, indexFile):
        """Prepare files and start IGV without blocking the Scipion event loop."""
        batchFile = None
        try:
            if gtfFile:
                gtfFile = self._prepareAnnotation(gtfFile, job)
            if job['cancel'].is_set():
                return
            job['messages'].put(('status', 'Finding the initial locus...'))
            batchFile = self._createIgvBatchFile(alignmentFile, fastaFile, gtfFile, indexFile)
            if job['cancel'].is_set():
                os.remove(batchFile)
                return
            job['messages'].put(('status', 'Starting IGV...'))
            self._launchIgv(batchFile, job)
            job['messages'].put(('done', ''))
        except Exception as error:
            if batchFile and os.path.isfile(batchFile):
                os.remove(batchFile)
            if not job['cancel'].is_set():
                logging.getLogger(__name__).exception('Could not open the IGV viewer')
                job['messages'].put(('error', str(error)))

    def _launchIgv(self, batchFile, job):
        command = (
            Plugin.getEnvActivationCommand(RNASEQ_DIC) +
            ' && exec igv -b {}'.format(shlex.quote(batchFile))
        )

        with tempfile.NamedTemporaryFile(prefix='scipion_igv_', suffix='.log', delete=False) as log:
            logPath = log.name
            process = subprocess.Popen(command, shell=True, executable='/bin/bash',
                                       env=Plugin.getEnviron(), stdout=log, stderr=subprocess.STDOUT,
                                       start_new_session=True)
        try:
            process.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass
        if job['cancel'].is_set():
            self._stopProcess(process)
            return
        if process.poll() not in (None, 0):
            with open(logPath) as log:
                detail = log.read()[-6000:]
            raise RuntimeError('IGV failed to start. Log: {}\n{}'.format(logPath, detail))
        logging.getLogger(__name__).info('IGV log: %s', logPath)
        threading.Thread(target=self._removeBatchAfterExit, args=(process, batchFile), daemon=True).start()

    @staticmethod
    def _stopProcess(process):
        """Stop only a viewer command's own process group, never protocols."""
        if process.poll() is not None:
            return
        try:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
        except ProcessLookupError:
            pass

    def _runPreparationCommand(self, command, job):
        process = subprocess.Popen(command, shell=True, executable='/bin/bash',
                                   env=Plugin.getEnviron(), stdout=subprocess.PIPE,
                                   stderr=subprocess.PIPE, text=True, start_new_session=True)
        while True:
            if job['cancel'].is_set():
                self._stopProcess(process)
                process.communicate()
                raise RuntimeError('IGV preparation cancelled.')
            try:
                output, errors = process.communicate(timeout=0.2)
                break
            except subprocess.TimeoutExpired:
                continue
        if process.returncode:
            raise subprocess.CalledProcessError(process.returncode, command, output, errors)

    def _prepareAnnotation(self, gtfFile, job):
        """Prepare a BGZF/Tabix annotation without loading the whole GTF in Java."""
        source = os.path.realpath(gtfFile)
        sourceStat = os.stat(source)
        signature = repr(('tabix-v1', source, sourceStat.st_size, sourceStat.st_mtime_ns))
        key = hashlib.sha256(signature.encode()).hexdigest()
        cacheRoot = os.path.join(tempfile.gettempdir(), 'scipion_igv_annotations_{}'.format(os.getuid()))
        os.makedirs(cacheRoot, mode=0o700, exist_ok=True)
        annotation = os.path.join(cacheRoot, key, 'annotation.gtf.gz')
        index = annotation + '.tbi'

        # Different viewers may prepare the same annotation concurrently.
        with open(os.path.join(cacheRoot, key + '.lock'), 'a') as lock:
            while True:
                if job['cancel'].is_set():
                    raise RuntimeError('IGV preparation cancelled.')
                try:
                    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    break
                except BlockingIOError:
                    job['messages'].put(('status', 'Waiting for another viewer to finish indexing...'))
                    job['cancel'].wait(0.5)
            if os.path.isfile(annotation) and os.path.isfile(index) and os.path.getsize(index) > 0:
                return annotation
            os.makedirs(os.path.dirname(annotation), exist_ok=True)
            activation = Plugin.getEnvActivationCommand(RNASEQ_DIC)
            try:
                with tempfile.TemporaryDirectory(dir=cacheRoot) as staging:
                    sortedGtf = os.path.join(staging, 'annotation.gtf.gz')
                    job['messages'].put(('status', 'Sorting and compressing reference annotation...'))
                    # External sorting bounds RAM usage even for multi-GB GTFs.
                    sort = "LC_ALL=C sort -t $'\\t' -k1,1 -k4,4n -S 256M -T " + shlex.quote(staging)
                    with open(source, 'rb') as handle:
                        compressed = handle.read(2) == b'\x1f\x8b'
                    if compressed:
                        sort = 'gzip -dc -- ' + shlex.quote(source) + ' | ' + sort
                    else:
                        sort += ' -- ' + shlex.quote(source)
                    command = (activation + ' && set -o pipefail && ' + sort +
                               ' | bgzip -@ 2 -c > ' + shlex.quote(sortedGtf))
                    self._runPreparationCommand(command, job)
                    self._checkFile(sortedGtf, 'Compressed reference GTF')
                    job['messages'].put(('status', 'Creating the reference annotation Tabix index...'))
                    command = activation + ' && tabix -p gff ' + shlex.quote(sortedGtf)
                    self._runPreparationCommand(command, job)
                    self._checkFile(sortedGtf + '.tbi', 'Reference GTF Tabix index')
                    os.replace(sortedGtf, annotation)
                    os.replace(sortedGtf + '.tbi', index)
            except subprocess.CalledProcessError as error:
                raise RuntimeError('Could not prepare the reference GTF for IGV:\n{}'.format(
                    error.stderr or error.stdout or str(error))) from error
        return annotation

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
                'load {} index={}'.format(self._batchPath(gtfFile), self._batchPath(gtfFile + '.tbi'))
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
                executable='/bin/bash',
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
