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

import gzip
import os
import re
import shlex
import shutil
import stat
import struct
import tempfile
import zlib

from pwchem import Plugin
from pwchem.constants import CUFFLINKS_DIC, RNASEQ_DIC
from pwchem.utils import assertHandle


HTML_EXTENSION = '.html'


def openFastq(fn):
    """Open a plain FASTQ or a .gz FASTQ as text; accept path-like objects."""
    fn = os.fspath(fn)
    return gzip.open(fn, 'rt') if fn.endswith('.gz') else open(fn, 'r')


def getFastqStats(fn):
    """Count FASTQ records and calculate their rounded mean sequence length.

    Files are streamed without loading all reads into memory. Standard
    four-line records are required: an @ header, sequence, + separator and
    quality string of the same length. Truncated or malformed records raise
    RuntimeError with their record number. Wrapped FASTQ is not supported.
    An empty file returns (0, 0), including an empty gzip stream.
    """
    numReads = 0
    totalLength = 0
    try:
        with openFastq(fn) as handle:
            while True:
                header = handle.readline()
                if not header:
                    break
                record = numReads + 1
                sequence = handle.readline()
                separator = handle.readline()
                quality = handle.readline()
                if not header.startswith('@') or not separator.startswith('+'):
                    raise RuntimeError('Invalid FASTQ header or separator in {} at record {}.'.format(fn, record))
                if not sequence or not quality:
                    raise RuntimeError('Truncated FASTQ record {} in {}.'.format(record, fn))
                length = len(sequence.rstrip('\r\n'))
                if length != len(quality.rstrip('\r\n')):
                    raise RuntimeError('Sequence and quality lengths differ in {} at record {}.'.format(fn, record))
                numReads += 1
                totalLength += length
    except (OSError, EOFError, UnicodeError, zlib.error) as error:
        raise RuntimeError('Could not read FASTQ {}: {}'.format(fn, error)) from error
    return numReads, int(round(totalLength / numReads)) if numReads else 0


def runFastqc(protocol, fastqFiles):
    """Run FastQC and return non-empty HTML reports in input order.

    Parameters
    ----------
    protocol
        Scipion protocol executing FastQC.
    fastqFiles : iterable of paths
        Non-empty FASTQ files. Names must produce distinct report basenames,
        even when inputs are located in different directories.

    Returns
    -------
    list
        Final HTML report paths under the protocol fastqc directory.

    Notes
    -----
    Reports are generated in a temporary directory. All HTML files are checked
    before publication, so reports from an earlier run cannot hide a missing
    result. Generated ZIP reports are retained beside their HTML files.
    Input files are not changed. An empty input list is an error.
    """
    fastqFiles = [os.fspath(path) for path in fastqFiles]
    reportNames = getFastqcReportNames(fastqFiles)
    outDir = protocol._getExtraPath('fastqc')
    os.makedirs(outDir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='run_', dir=outDir) as tmpDir:
        arguments = ' '.join(shlex.quote(path) for path in ['-o', tmpDir, *fastqFiles])
        Plugin.runCondaCommand(protocol, arguments, RNASEQ_DIC, 'fastqc')
        publishFastqcReports(reportNames, tmpDir, outDir)
    return [os.path.join(outDir, name + HTML_EXTENSION) for name in reportNames]


def getFastqcReportNames(fastqFiles):
    if not fastqFiles:
        raise ValueError('At least one FASTQ file is required for FastQC.')
    reportNames = []
    for path in fastqFiles:
        if not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise RuntimeError('FastQC input is missing or empty: {}'.format(path))
        baseName = os.path.basename(path)
        for extension in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if baseName.endswith(extension):
                baseName = baseName[:-len(extension)]
                break
        reportNames.append(baseName + '_fastqc')
    if len(set(reportNames)) != len(reportNames):
        raise ValueError('FASTQ inputs produce duplicate FastQC report names. Use distinct filenames.')
    return reportNames


def publishFastqcReports(reportNames, tmpDir, outDir):
    for name in reportNames:
        report = os.path.join(tmpDir, name + HTML_EXTENSION)
        if not os.path.isfile(report) or os.path.getsize(report) == 0:
            raise RuntimeError('FastQC HTML report is missing or empty: {}'.format(name))
    for name in reportNames:
        for extension in (HTML_EXTENSION, '.zip'):
            source = os.path.join(tmpDir, name + extension)
            target = os.path.join(outDir, name + extension)
            if os.path.isfile(source) and os.path.getsize(source) > 0:
                os.replace(source, target)
            elif extension == '.zip' and os.path.isfile(target):
                os.remove(target)


def assertOutputExists(test, protocol, output):
    """Check that the protocol generated an output object."""
    assertHandle(
        test.assertIsNotNone,
        output,
        cwd=protocol.getWorkingDir()
    )


def assertGenomeFiles(test, protocol, genome, expectGtf=True):
    """Check that required genome files are regular, non-empty files."""
    paths = [genome.getFastaFile()]
    if expectGtf:
        assertHandle(test.assertTrue, genome.hasGtfFile(), cwd=protocol.getWorkingDir())
        paths.append(genome.getGtfFile())
    for path in paths:
        assertHandle(test.assertTrue, bool(path) and os.path.isfile(path), cwd=protocol.getWorkingDir())
        assertHandle(test.assertGreater, os.path.getsize(path), 0, cwd=protocol.getWorkingDir())


class GATKPreparationMixin:
    """Prepare local input links and compatible indexes without changing source files.

    BAM and reference dictionaries are compared by sequence names, lengths and
    order. Resource VCF headers are inspected before indexing. Copies of existing
    indexes are published atomically and both BAI naming conventions are handled.
    Signature and timestamp checks are deliberately lightweight: they do not prove
    that every index entry or alignment record is valid, or compare reference MD5s.
    Subclasses supply protocol paths, thread parameters and the command-log path.
    """

    def prepareReference(self, reference):
        self.checkNonempty(reference, 'Reference FASTA')
        with open(reference, 'rb') as handle:
            if not handle.readline().lstrip().startswith(b'>'):
                raise RuntimeError('Reference must be an uncompressed FASTA: {}'.format(reference))

        preparedReference = self.getPreparedReference()
        self.stageFile(reference, preparedReference)
        self.reuseIndex(reference + '.fai', preparedReference + '.fai', reference,
                         self.validFastaIndex)
        self.reuseIndex(self.getSequenceDictionaryPath(reference),
                         self.getSequenceDictionaryPath(preparedReference), reference,
                         self.validSequenceDictionary)
        self.ensureFastaIndex(preparedReference)
        self.ensureSequenceDictionary(preparedReference)
        contigs = self.readFastaIndex(preparedReference + '.fai')
        dictionary = self.readSequenceDictionary(self.getSequenceDictionaryPath(preparedReference))
        if list(contigs.items()) != list(dictionary.items()):
            self.warning('Reference auxiliary files disagree; rebuilding the local FASTA index and dictionary.')
            self.removeIndex(preparedReference + '.fai')
            self.removeIndex(self.getSequenceDictionaryPath(preparedReference))
            self.ensureFastaIndex(preparedReference)
            self.ensureSequenceDictionary(preparedReference)
            contigs = self.readFastaIndex(preparedReference + '.fai')
            dictionary = self.readSequenceDictionary(self.getSequenceDictionaryPath(preparedReference))
            if list(contigs.items()) != list(dictionary.items()):
                raise RuntimeError('Generated FASTA index and dictionary disagree in sequence names, '
                                   'lengths or order.')
        return contigs

    def prepareAlignment(self, alignment, target, referenceContigs, label, requireReadGroups=False):
        self.stageFile(alignment.getFileName(), target)
        self.runTool('samtools', ['quickcheck', '-v', target])
        header = target + '.header.sam'
        self.runTool('samtools', ['view', '-H', '-o', header, target])
        contigs, groups, sortOrder = self.readSamHeader(header)
        self.checkContigs(referenceContigs, contigs, label, requireAll=True)
        if sortOrder in ('queryname', 'unsorted'):
            raise RuntimeError('{} declares {} sorting; coordinate sorting is required.'.format(label, sortOrder))
        if requireReadGroups and (not groups or any(
                not group.get('ID', '').strip() or not group.get('SM', '').strip() for group in groups)):
            raise RuntimeError('{} requires read groups with non-empty ID and SM tags. '
                               'Prepare read groups with Picard first.'.format(label))
        ids = [group.get('ID') for group in groups]
        if len(ids) != len(set(ids)):
            raise RuntimeError('{} contains duplicate read-group IDs.'.format(label))
        for index in (alignment.getIndexFile(), alignment.getFileName() + '.bai',
                      os.path.splitext(alignment.getFileName())[0] + '.bai'):
            self.reuseIndex(index, target + '.bai', alignment.getFileName(),
                             lambda path: self.hasIndexSignature(path, b'BAI\x01'))
        self.ensureBamIndex(target)
        return {group['SM'] for group in groups if group.get('SM')}

    def getLocalVCFPath(self, source, stem):
        self.checkNonempty(source, 'VCF')
        with open(source, 'rb') as handle:
            compressed = handle.read(2) == b'\x1f\x8b'
        return self._getExtraPath('inputs', stem + ('.vcf.gz' if compressed else '.vcf'))

    def prepareVCF(self, vcf, target, referenceContigs, label, validateHeader=None):
        source = vcf.getFileName()
        self.stageFile(source, target)
        contigs = self.readVCFContigs(target)
        if contigs:
            self.checkContigs(referenceContigs, contigs, label)
        else:
            self.warning('{} has no contig declarations; reference compatibility could not be checked.'.format(label))
        if validateHeader is not None:
            validateHeader(target)
        extension = '.tbi' if target.endswith('.gz') else '.idx'
        for index in (vcf.getIndexFile(), source + extension):
            self.reuseIndex(index, target + extension, source,
                             lambda path: self.validFeatureIndex(path, target))
        self.ensureFeatureIndex(target)

    def runTool(self, program, arguments):
        command = ' '.join(self.quote(value) for value in arguments)
        environment = CUFFLINKS_DIC if program in ('cufflinks', 'cuffdiff') else RNASEQ_DIC
        Plugin.runCondaCommand(self, command, environment, program)
        self.appendCommand('{} {}'.format(program, command))

    def appendCommand(self, command):
        with open(self._getCommandFile(), 'a') as outputFile:
            outputFile.write(command + '\n')

    def readCommandHistory(self):
        commandFile = self._getCommandFile()

        if not os.path.isfile(commandFile):
            return ''

        with open(commandFile) as inputFile:
            commands = [line.strip() for line in inputFile if line.strip()]

        return ' && '.join(commands)

    def getPreparedReference(self):
        return self._getExtraPath('reference', 'reference.fa')

    def stageFile(self, source, target):
        self.checkNonempty(source, 'Input file')
        source = os.path.abspath(source)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        if os.path.lexists(target):
            if os.path.realpath(target) == os.path.realpath(source):
                return
            os.remove(target)
            for path in (target + '.fai', target + '.bai', target + '.idx', target + '.tbi',
                         target + '.csi', target + '.gzi', os.path.splitext(target)[0] + '.bai',
                         self.getSequenceDictionaryPath(target)):
                if os.path.isfile(path):
                    os.remove(path)
        os.symlink(source, target)

    @staticmethod
    def isCurrentIndex(index, dataFile):
        if not index:
            return False
        try:
            indexInfo = os.stat(index)
            dataInfo = os.stat(dataFile)
        except OSError:
            return False
        return (stat.S_ISREG(indexInfo.st_mode) and indexInfo.st_size > 0
                and indexInfo.st_mtime_ns >= dataInfo.st_mtime_ns)

    @staticmethod
    def copyIndex(source, target):
        """Replace a local index only after its copy has finished."""
        partialFile = target + '.part'
        try:
            shutil.copy2(source, partialFile)
            os.replace(partialFile, target)
        finally:
            if os.path.exists(partialFile):
                os.remove(partialFile)

    def reuseIndex(self, source, target, dataFile, validator):
        if self.isCurrentIndex(target, dataFile) and validator(target):
            return
        if self.isCurrentIndex(source, dataFile) and validator(source):
            if os.path.abspath(source) != os.path.abspath(target):
                self.copyIndex(source, target)

    @staticmethod
    def hasIndexSignature(path, signature):
        try:
            with open(path, 'rb') as handle:
                compressed = handle.read(2) == b'\x1f\x8b'
            opener = gzip.open if compressed else open
            with opener(path, 'rb') as handle:
                return handle.read(4) == signature and len(handle.read(4)) == 4
        except (OSError, EOFError, zlib.error):
            return False

    @staticmethod
    def readFastaIndex(path):
        contigs = {}
        fasta = path[:-4]
        fastaSize = os.path.getsize(fasta) if os.path.isfile(fasta) else None
        with open(path) as handle:
            for line in handle:
                fields = line.rstrip().split('\t')
                if len(fields) != 5:
                    raise ValueError('Invalid FASTA index record: {}'.format(path))
                name = fields[0]
                length, offset, bases, width = map(int, fields[1:])
                if (not name or name in contigs or length <= 0 or offset < 0
                        or bases <= 0 or width < bases):
                    raise ValueError('Invalid FASTA index values: {}'.format(path))
                if fastaSize is not None:
                    lastBase = offset + ((length - 1) // bases) * width + (length - 1) % bases
                    if lastBase >= fastaSize:
                        raise ValueError('FASTA index points beyond the reference file: {}'.format(path))
                contigs[name] = length
        if not contigs:
            raise ValueError('Empty FASTA index: {}'.format(path))
        return contigs

    @classmethod
    def readSamHeader(cls, path):
        """Read sequence lengths, read groups and sort order in a single pass."""
        contigs, groups, sortOrder = {}, [], None
        with open(path) as handle:
            for line in handle:
                record, *fields = line.rstrip('\r\n').split('\t')
                if record not in ('@SQ', '@RG', '@HD'):
                    continue
                try:
                    tags = cls.parseSamTags(fields)
                    if record == '@SQ':
                        name, length = tags['SN'], int(tags['LN'])
                        cls.checkSamSequence(name, length, contigs)
                        contigs[name] = length
                    elif record == '@RG':
                        groups.append(tags)
                    else:
                        sortOrder = tags.get('SO')
                except (ValueError, KeyError) as error:
                    raise ValueError('Invalid {} record in {}'.format(record, path)) from error
        if not contigs:
            raise ValueError('Sequence dictionary has no @SQ records: {}'.format(path))
        return contigs, groups, sortOrder

    @staticmethod
    def parseSamTags(fields):
        pairs = [field.split(':', 1) for field in fields]
        tags = dict(pairs)
        if len(tags) != len(pairs):
            raise ValueError('Duplicate SAM tag')
        return tags

    @staticmethod
    def checkSamSequence(name, length, contigs):
        if not name or name in contigs or length <= 0:
            raise ValueError('Empty/duplicate sequence name or invalid length')

    @classmethod
    def readSequenceDictionary(cls, path):
        return cls.readSamHeader(path)[0]

    def validFastaIndex(self, path):
        try:
            self.readFastaIndex(path)
            return True
        except (OSError, ValueError):
            return False

    def validSequenceDictionary(self, path):
        try:
            self.readSequenceDictionary(path)
            return True
        except (OSError, ValueError, KeyError):
            return False

    @staticmethod
    def removeIndex(path):
        if os.path.lexists(path):
            os.remove(path)

    def ensureFastaIndex(self, referenceFasta):
        index = referenceFasta + '.fai'
        if self.isCurrentIndex(index, referenceFasta) and self.validFastaIndex(index):
            return
        self.removeIndex(index)
        self.runTool('samtools', ['faidx', referenceFasta])
        if not self.validFastaIndex(index):
            raise RuntimeError('samtools did not create a valid FASTA index: {}'.format(index))

    def ensureSequenceDictionary(self, referenceFasta):
        dictionary = self.getSequenceDictionaryPath(referenceFasta)
        if (self.isCurrentIndex(dictionary, referenceFasta)
                and self.validSequenceDictionary(dictionary)):
            return
        self.removeIndex(dictionary)
        self.runTool('gatk', ['CreateSequenceDictionary', '-R', referenceFasta, '-O', dictionary])
        if not self.validSequenceDictionary(dictionary):
            raise RuntimeError('GATK did not create a valid sequence dictionary: {}'.format(dictionary))

    def validFeatureIndex(self, index, featureFile):
        signature = b'TBI\x01' if featureFile.endswith('.gz') else b'TIDX'
        return self.hasIndexSignature(index, signature)

    def ensureFeatureIndex(self, featureFile):
        index = featureFile + ('.tbi' if featureFile.endswith('.gz') else '.idx')
        if self.isCurrentIndex(index, featureFile) and self.validFeatureIndex(index, featureFile):
            return
        self.removeIndex(index)
        self.runTool('gatk', ['IndexFeatureFile', '-I', featureFile])
        if not self.validFeatureIndex(index, featureFile):
            raise RuntimeError('GATK did not create a compatible VCF index: {}'.format(index))

    def clearBamOutput(self, path):
        """Remove a generated BAM and its two standard BAI names before rerunning a step."""
        for output in (path, path + '.bai', os.path.splitext(path)[0] + '.bai'):
            self.removeIndex(output)


    def ensureBamIndex(self, bamFile):
        self.validateBamFile(bamFile)
        index = bamFile + '.bai'
        if self.isCurrentIndex(index, bamFile) and self.hasIndexSignature(index, b'BAI\x01'):
            return
        alternate = os.path.splitext(bamFile)[0] + '.bai'
        if self.isCurrentIndex(alternate, bamFile) and self.hasIndexSignature(alternate, b'BAI\x01'):
            self.copyIndex(alternate, index)
            return
        self.removeIndex(index)
        self.runTool('samtools', ['quickcheck', '-v', bamFile])
        self.runTool('samtools', ['index', '-@', self.numberOfThreads.get(), bamFile, index])
        if not self.hasIndexSignature(index, b'BAI\x01'):
            raise RuntimeError('samtools did not create a valid BAM index: {}'.format(index))

    @classmethod
    def readVCFContigs(cls, path, compressed=None):
        """Validate VCF headers; allow explicit compression for temporary download paths."""
        path = os.fspath(path)
        if compressed is None:
            compressed = path.endswith('.gz')
        if compressed:
            cls.checkBgzfHeader(path)
        opener = gzip.open if compressed else open
        contigs = {}
        with opener(path, 'rt') as handle:
            if not handle.readline().startswith('##fileformat=VCFv'):
                raise RuntimeError('Input is not a VCF: {}'.format(path))
            for line in handle:
                if line.startswith('##contig=<'):
                    name, length = cls.parseVcfContig(line, path, contigs)
                    contigs[name] = length
                if line.startswith('#CHROM\t'):
                    cls.checkVcfColumns(line, path)
                    return contigs
                if not line.startswith('##'):
                    break
        raise RuntimeError('VCF has no column header: {}'.format(path))

    @classmethod
    def checkBgzfHeader(cls, path):
        with open(path, 'rb') as handle:
            header = handle.read(12)
            if len(header) != 12 or header[:3] != b'\x1f\x8b\x08' or not header[3] & 4:
                raise RuntimeError('VCF must be BGZF-compressed, not ordinary gzip: {}'.format(path))
            extraSize = struct.unpack('<H', header[10:12])[0]
            extra = handle.read(extraSize)
        if len(extra) != extraSize:
            raise RuntimeError('Truncated BGZF header: {}'.format(path))
        blockSize = cls.getBgzfBlockSize(extra, path)
        if blockSize is None or blockSize < 12 + extraSize + 8:
            raise RuntimeError('VCF must use valid BGZF blocks for tabix indexing: {}'.format(path))
        if blockSize > os.path.getsize(path):
            raise RuntimeError('Truncated first BGZF block: {}'.format(path))

    @staticmethod
    def getBgzfBlockSize(extra, path):
        offset, blockSize = 0, None
        while offset < len(extra):
            if offset + 4 > len(extra):
                raise RuntimeError('Malformed BGZF extra fields: {}'.format(path))
            size = struct.unpack('<H', extra[offset + 2:offset + 4])[0]
            end = offset + 4 + size
            if end > len(extra):
                raise RuntimeError('Truncated BGZF extra field: {}'.format(path))
            if extra[offset:offset + 2] == b'BC':
                if size != 2 or blockSize is not None:
                    raise RuntimeError('Invalid BGZF BC field: {}'.format(path))
                blockSize = struct.unpack('<H', extra[offset + 4:end])[0] + 1
            offset = end
        return blockSize

    @classmethod
    def parseVcfContig(cls, line, path, contigs):
        declaration = line.rstrip('\r\n')
        if not declaration.endswith('>'):
            raise RuntimeError('Unterminated VCF contig declaration: {}'.format(path))
        fields = {}
        pattern = r'([^=,<>\s]+)=("(?:\\.|[^"\\])*"|[^,>]*)'
        for match in re.finditer(pattern, declaration[len('##contig=<'):-1]):
            key, value = match.group(1).lower(), match.group(2).strip('"')
            if key in fields:
                raise RuntimeError('Duplicate VCF contig field {} in {}'.format(key, path))
            fields[key] = value
        name = fields.get('id')
        if not name or name in contigs:
            raise RuntimeError('Empty or duplicate VCF contig ID in {}'.format(path))
        return name, cls.getVcfContigLength(fields, name)

    @staticmethod
    def getVcfContigLength(fields, name):
        if 'length' not in fields:
            return None
        try:
            length = int(fields['length'])
            if length <= 0:
                raise ValueError
        except ValueError as error:
            raise RuntimeError('Invalid VCF contig length for {}'.format(name)) from error
        return length

    @staticmethod
    def checkVcfColumns(line, path):
        columns = line.rstrip('\r\n').split('\t')
        if columns[:8] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
            raise RuntimeError('VCF has an invalid column header: {}'.format(path))
        if len(columns) == 8:
            return
        samples = columns[9:]
        if (len(columns) < 10 or columns[8] != 'FORMAT'
                or any(not name for name in samples) or len(set(samples)) != len(samples)):
            raise RuntimeError('VCF has invalid or duplicate sample columns: {}'.format(path))

    def checkContigs(self, reference, contigs, label, requireAll=False):
        shared = reference.keys() & contigs.keys()
        if not shared:
            raise RuntimeError('{} and reference FASTA have no sequence names in common.'.format(label))
        for name, length in contigs.items():
            if name not in reference:
                if requireAll:
                    raise RuntimeError('{} sequence {} is absent from the reference FASTA.'.format(label, name))
                continue
            if length is not None and length != reference[name]:
                raise RuntimeError('{} sequence {} has length {}, but the reference has {}.'.format(
                    label, name, length, reference[name]
                ))
        if requireAll and list(contigs) != [name for name in reference if name in contigs]:
            raise RuntimeError('{} sequence order differs from the reference. Reorder the BAM '
                               'with Picard ReorderSam first.'.format(label))
        if not requireAll and contigs.keys() - reference.keys():
            self.warning('{} contains extra contigs absent from the reference; GATK will perform '
                         'its own sequence dictionary validation.'.format(label))

    @staticmethod
    def getReferenceFasta(alignment):
        return alignment.getReferenceFasta() if alignment is not None else None

    @staticmethod
    def getSequenceDictionaryPath(referenceFasta):
        root, extension = os.path.splitext(referenceFasta)

        if extension.lower() in ('.fa', '.fasta', '.fna'):
            return root + '.dict'

        return referenceFasta + '.dict'

    @classmethod
    def validateBamFile(cls, bamFile):
        cls.checkNonempty(bamFile, 'BAM')

    @staticmethod
    def quote(value):
        return shlex.quote(str(value))

    @staticmethod
    def checkNonempty(path, label):
        if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise RuntimeError('{} is missing or empty: {}'.format(label, path))


def iterVepRows(path):
    """Yield native VEP tab rows, expanding Extra without discarding its text.

    Accept comment-only metadata and header-only tables, but require the
    standard VEP fields. Both Extra ENSP=... and explicit ENSP columns work.
    Malformed rows are errors rather than silently truncated annotations.
    """
    columns = None
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('##'):
                continue
            if line.startswith('#'):
                if columns is not None:
                    raise RuntimeError('Repeated VEP header at line {}.'.format(number))
                columns = parseVepColumns(line)
                continue
            if columns is None:
                raise RuntimeError('VEP data appears before the column header.')
            yield parseVepRow(line, columns, number)
    if columns is None:
        raise RuntimeError('VEP column header not found.')


def parseVepColumns(line):
    required = {'Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
                'Feature_type', 'Consequence', 'Protein_position', 'Amino_acids'}
    columns = line.rstrip('\r\n').lstrip('#').split('\t')
    if len(columns) != len(set(columns)) or not required.issubset(columns):
        raise RuntimeError('Missing or duplicate required VEP columns.')
    return columns


def parseVepRow(line, columns, number):
    values = line.rstrip('\r\n').split('\t')
    if len(values) != len(columns):
        raise RuntimeError('VEP line {} has the wrong number of columns.'.format(number))
    row = dict(zip(columns, values))
    for item in row.get('Extra', '').split(';'):
        if '=' in item:
            key, value = item.split('=', 1)
            if key not in row:
                row[key] = value
    return row
