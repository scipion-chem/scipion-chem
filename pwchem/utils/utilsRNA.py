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
from pwchem.constants import RNASEQ_DIC
from pwchem.utils import assertHandle


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

    outDir = protocol._getExtraPath('fastqc')
    os.makedirs(outDir, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='run_', dir=outDir) as tmpDir:
        arguments = ' '.join(shlex.quote(path) for path in ['-o', tmpDir, *fastqFiles])
        Plugin.runCondaCommand(protocol, arguments, RNASEQ_DIC, 'fastqc')
        for name in reportNames:
            report = os.path.join(tmpDir, name + '.html')
            if not os.path.isfile(report) or os.path.getsize(report) == 0:
                raise RuntimeError('FastQC HTML report is missing or empty: {}'.format(name))
        for name in reportNames:
            for extension in ('.html', '.zip'):
                source = os.path.join(tmpDir, name + extension)
                target = os.path.join(outDir, name + extension)
                if os.path.isfile(source) and os.path.getsize(source) > 0:
                    os.replace(source, target)
                elif extension == '.zip' and os.path.isfile(target):
                    os.remove(target)
    return [os.path.join(outDir, name + '.html') for name in reportNames]


# ============================================================================
# Species utilities
# ============================================================================
COMMON_SPECIES = {
    "Homo sapiens": {
        "ensembl": "homo_sapiens",
        "ncbi": "Homo sapiens",
        "assembly": "GRCh38",
    },
    "Mus musculus": {
        "ensembl": "mus_musculus",
        "ncbi": "Mus musculus",
        "assembly": "GRCm39",
    },
    "Rattus norvegicus": {
        "ensembl": "rattus_norvegicus",
        "ncbi": "Rattus norvegicus",
        "assembly": "mRatBN7.2",
    },
    "Danio rerio": {
        "ensembl": "danio_rerio",
        "ncbi": "Danio rerio",
        "assembly": "GRCz11",
    },
    "Drosophila melanogaster": {
        "ensembl": "drosophila_melanogaster",
        "ncbi": "Drosophila melanogaster",
        "assembly": "BDGP6.46",
    },
    "Caenorhabditis elegans": {
        "ensembl": "caenorhabditis_elegans",
        "ncbi": "Caenorhabditis elegans",
        "assembly": "WBcel235",
    },
    "Saccharomyces cerevisiae": {
        "ensembl": "saccharomyces_cerevisiae",
        "ncbi": "Saccharomyces cerevisiae",
        "assembly": "R64-1-1",
    },
    "Arabidopsis thaliana": {
        "ensembl": "arabidopsis_thaliana",
        "ncbi": "Arabidopsis thaliana",
        "assembly": "TAIR10",
    },
}

def getCommonSpecies():
    """Return the scientific names of the predefined common species.

    Returns
    -------
    list
        Scientific species names available as common genomic resources.
    """
    return list(COMMON_SPECIES.keys())

def getProviderSpeciesName(species, provider):
    """Return the species identifier expected by a genomic data provider.

    Parameters
    ----------
    species : str
        Scientific species name, for example ``"Homo sapiens"``.
    provider : str
        Genomic data provider. Currently supported values are
        ``"ensembl"`` and ``"ncbi"``.

    Returns
    -------
    str
        Species identifier expected by the selected provider.

    Notes
    -----
    If the species is not included in ``COMMON_SPECIES``, a generic
    conversion is used for Ensembl and the original scientific name is
    preserved for NCBI. This allows custom species to be used without
    requiring them to be predefined. Whitespace is normalized and common
    species names are recognized regardless of capitalization. Empty names
    and unsupported providers raise ValueError.
    """
    provider = (provider or '').strip().lower()
    species = ' '.join((species or '').split())
    if not species:
        raise ValueError('A scientific species name is required.')

    if provider not in ("ensembl", "ncbi"):
        raise ValueError(
            f"Unsupported genomic data provider: {provider}"
        )

    commonName = next((name for name in COMMON_SPECIES if name.casefold() == species.casefold()), None)
    speciesInfo = COMMON_SPECIES.get(commonName)

    if speciesInfo is not None:
        return speciesInfo[provider]

    if provider == "ensembl":
        return species.lower().replace(" ", "_")

    return species


def parseCustomSpecies(speciesText):
    """Parse a comma- or semicolon-separated list of custom species.

    Parameters
    ----------
    speciesText : str
        Species entered by the user, separated by commas or semicolons.
        For example::

            Homo sapiens;Mus musculus;Canis lupus familiaris

    Returns
    -------
    list
        Cleaned species names. Empty entries are ignored and duplicate
        species are removed case-insensitively while preserving their order.
        Repeated whitespace is collapsed and common species use their
        canonical scientific-name capitalization.
    """
    if not speciesText:
        return []

    commonNames = {name.casefold(): name for name in COMMON_SPECIES}
    selected = {}
    for item in re.split(r'[;,]', speciesText):
        name = ' '.join(item.split())
        if name:
            key = name.casefold()
            selected.setdefault(key, commonNames.get(key, name))
    return list(selected.values())


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

    def _prepareReference(self, reference):
        self._checkNonempty(reference, 'Reference FASTA')
        with open(reference, 'rb') as handle:
            if not handle.readline().lstrip().startswith(b'>'):
                raise RuntimeError('Reference must be an uncompressed FASTA: {}'.format(reference))

        preparedReference = self._getPreparedReference()
        self._stageFile(reference, preparedReference)
        self._reuseIndex(reference + '.fai', preparedReference + '.fai', reference,
                         self._validFastaIndex)
        self._reuseIndex(self._getSequenceDictionaryPath(reference),
                         self._getSequenceDictionaryPath(preparedReference), reference,
                         self._validSequenceDictionary)
        self._ensureFastaIndex(preparedReference)
        self._ensureSequenceDictionary(preparedReference)
        contigs = self._readFastaIndex(preparedReference + '.fai')
        dictionary = self._readSequenceDictionary(self._getSequenceDictionaryPath(preparedReference))
        if list(contigs.items()) != list(dictionary.items()):
            self.warning('Reference auxiliary files disagree; rebuilding the local FASTA index and dictionary.')
            self._removeIndex(preparedReference + '.fai')
            self._removeIndex(self._getSequenceDictionaryPath(preparedReference))
            self._ensureFastaIndex(preparedReference)
            self._ensureSequenceDictionary(preparedReference)
            contigs = self._readFastaIndex(preparedReference + '.fai')
            dictionary = self._readSequenceDictionary(self._getSequenceDictionaryPath(preparedReference))
            if list(contigs.items()) != list(dictionary.items()):
                raise RuntimeError('Generated FASTA index and dictionary disagree in sequence names, '
                                   'lengths or order.')
        return contigs

    def _prepareAlignment(self, alignment, target, referenceContigs, label, requireReadGroups=False):
        self._stageFile(alignment.getFileName(), target)
        self._runTool('samtools', ['quickcheck', '-v', target])
        header = target + '.header.sam'
        self._runTool('samtools', ['view', '-H', '-o', header, target])
        contigs, groups, sortOrder = self._readSamHeader(header)
        self._checkContigs(referenceContigs, contigs, label, requireAll=True)
        if sortOrder in ('queryname', 'unsorted'):
            raise RuntimeError('{} declares {} sorting; coordinate sorting is required.'.format(label, sortOrder))
        if requireReadGroups:
            if not groups or any(not group.get('ID', '').strip() or not group.get('SM', '').strip()
                                 for group in groups):
                raise RuntimeError('{} requires read groups with non-empty ID and SM tags. '
                                   'Prepare read groups with Picard first.'.format(label))
        ids = [group.get('ID') for group in groups]
        if len(ids) != len(set(ids)):
            raise RuntimeError('{} contains duplicate read-group IDs.'.format(label))
        for index in (alignment.getIndexFile(), alignment.getFileName() + '.bai',
                      os.path.splitext(alignment.getFileName())[0] + '.bai'):
            self._reuseIndex(index, target + '.bai', alignment.getFileName(),
                             lambda path: self._hasIndexSignature(path, b'BAI\x01'))
        self._ensureBamIndex(target)
        return {group['SM'] for group in groups if group.get('SM')}

    def _getLocalVCFPath(self, source, stem):
        self._checkNonempty(source, 'VCF')
        with open(source, 'rb') as handle:
            compressed = handle.read(2) == b'\x1f\x8b'
        return self._getExtraPath('inputs', stem + ('.vcf.gz' if compressed else '.vcf'))

    def _prepareVCF(self, vcf, target, referenceContigs, label, validateHeader=None):
        source = vcf.getFileName()
        self._stageFile(source, target)
        contigs = self._readVCFContigs(target)
        if contigs:
            self._checkContigs(referenceContigs, contigs, label)
        else:
            self.warning('{} has no contig declarations; reference compatibility could not be checked.'.format(label))
        if validateHeader is not None:
            validateHeader(target)
        extension = '.tbi' if target.endswith('.gz') else '.idx'
        for index in (vcf.getIndexFile(), source + extension):
            self._reuseIndex(index, target + extension, source,
                             lambda path: self._validFeatureIndex(path, target))
        self._ensureFeatureIndex(target)

    def _runTool(self, program, arguments):
        command = ' '.join(self._quote(value) for value in arguments)
        Plugin.runCondaCommand(self, command, RNASEQ_DIC, program)
        self._appendCommand('{} {}'.format(program, command))

    def _appendCommand(self, command):
        with open(self._getCommandFile(), 'a') as outputFile:
            outputFile.write(command + '\n')

    def _readCommandHistory(self):
        commandFile = self._getCommandFile()

        if not os.path.isfile(commandFile):
            return ''

        with open(commandFile) as inputFile:
            commands = [line.strip() for line in inputFile if line.strip()]

        return ' && '.join(commands)

    def _getPreparedReference(self):
        return self._getExtraPath('reference', 'reference.fa')

    def _stageFile(self, source, target):
        self._checkNonempty(source, 'Input file')
        source = os.path.abspath(source)
        os.makedirs(os.path.dirname(target), exist_ok=True)
        if os.path.lexists(target):
            if os.path.realpath(target) == os.path.realpath(source):
                return
            os.remove(target)
            for path in (target + '.fai', target + '.bai', target + '.idx', target + '.tbi',
                         target + '.csi', target + '.gzi', os.path.splitext(target)[0] + '.bai',
                         self._getSequenceDictionaryPath(target)):
                if os.path.isfile(path):
                    os.remove(path)
        os.symlink(source, target)

    @staticmethod
    def _isCurrentIndex(index, dataFile):
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
    def _copyIndex(source, target):
        """Replace a local index only after its copy has finished."""
        partialFile = target + '.part'
        try:
            shutil.copy2(source, partialFile)
            os.replace(partialFile, target)
        finally:
            if os.path.exists(partialFile):
                os.remove(partialFile)

    def _reuseIndex(self, source, target, dataFile, validator):
        if self._isCurrentIndex(target, dataFile) and validator(target):
            return
        if self._isCurrentIndex(source, dataFile) and validator(source):
            if os.path.abspath(source) != os.path.abspath(target):
                self._copyIndex(source, target)

    @staticmethod
    def _hasIndexSignature(path, signature):
        try:
            with open(path, 'rb') as handle:
                compressed = handle.read(2) == b'\x1f\x8b'
            opener = gzip.open if compressed else open
            with opener(path, 'rb') as handle:
                return handle.read(4) == signature and len(handle.read(4)) == 4
        except (OSError, EOFError, zlib.error):
            return False

    @staticmethod
    def _readFastaIndex(path):
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

    @staticmethod
    def _readSamHeader(path):
        """Read sequence lengths, read groups and sort order in a single pass."""
        contigs, groups, sortOrder = {}, [], None
        with open(path) as handle:
            for line in handle:
                record, *fields = line.rstrip('\r\n').split('\t')
                if record not in ('@SQ', '@RG', '@HD'):
                    continue
                try:
                    pairs = [field.split(':', 1) for field in fields]
                    tags = dict(pairs)
                    if len(tags) != len(pairs):
                        raise ValueError('Duplicate SAM tag')
                    if record == '@SQ':
                        name, length = tags['SN'], int(tags['LN'])
                        if not name or name in contigs or length <= 0:
                            raise ValueError('Empty/duplicate sequence name or invalid length')
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

    @classmethod
    def _readSequenceDictionary(cls, path):
        return cls._readSamHeader(path)[0]

    def _validFastaIndex(self, path):
        try:
            self._readFastaIndex(path)
            return True
        except (OSError, ValueError):
            return False

    def _validSequenceDictionary(self, path):
        try:
            self._readSequenceDictionary(path)
            return True
        except (OSError, ValueError, KeyError):
            return False

    @staticmethod
    def _removeIndex(path):
        if os.path.lexists(path):
            os.remove(path)

    def _ensureFastaIndex(self, referenceFasta):
        index = referenceFasta + '.fai'
        if self._isCurrentIndex(index, referenceFasta) and self._validFastaIndex(index):
            return
        self._removeIndex(index)
        self._runTool('samtools', ['faidx', referenceFasta])
        if not self._validFastaIndex(index):
            raise RuntimeError('samtools did not create a valid FASTA index: {}'.format(index))

    def _ensureSequenceDictionary(self, referenceFasta):
        dictionary = self._getSequenceDictionaryPath(referenceFasta)
        if (self._isCurrentIndex(dictionary, referenceFasta)
                and self._validSequenceDictionary(dictionary)):
            return
        self._removeIndex(dictionary)
        self._runTool('gatk', ['CreateSequenceDictionary', '-R', referenceFasta, '-O', dictionary])
        if not self._validSequenceDictionary(dictionary):
            raise RuntimeError('GATK did not create a valid sequence dictionary: {}'.format(dictionary))

    def _validFeatureIndex(self, index, featureFile):
        signature = b'TBI\x01' if featureFile.endswith('.gz') else b'TIDX'
        return self._hasIndexSignature(index, signature)

    def _ensureFeatureIndex(self, featureFile):
        index = featureFile + ('.tbi' if featureFile.endswith('.gz') else '.idx')
        if self._isCurrentIndex(index, featureFile) and self._validFeatureIndex(index, featureFile):
            return
        self._removeIndex(index)
        self._runTool('gatk', ['IndexFeatureFile', '-I', featureFile])
        if not self._validFeatureIndex(index, featureFile):
            raise RuntimeError('GATK did not create a compatible VCF index: {}'.format(index))

    def _clearBamOutput(self, path):
        """Remove a generated BAM and its two standard BAI names before rerunning a step."""
        for output in (path, path + '.bai', os.path.splitext(path)[0] + '.bai'):
            self._removeIndex(output)

    def _clearVCFOutput(self, path):
        """Remove a generated VCF and companion indexes before rerunning a step."""
        for output in (path, path + '.tbi', path + '.csi', path + '.idx'):
            self._removeIndex(output)

    def _ensureBamIndex(self, bamFile):
        self._validateBamFile(bamFile)
        index = bamFile + '.bai'
        if self._isCurrentIndex(index, bamFile) and self._hasIndexSignature(index, b'BAI\x01'):
            return
        alternate = os.path.splitext(bamFile)[0] + '.bai'
        if self._isCurrentIndex(alternate, bamFile) and self._hasIndexSignature(alternate, b'BAI\x01'):
            self._copyIndex(alternate, index)
            return
        self._removeIndex(index)
        self._runTool('samtools', ['quickcheck', '-v', bamFile])
        self._runTool('samtools', ['index', '-@', self.numberOfThreads.get(), bamFile, index])
        if not self._hasIndexSignature(index, b'BAI\x01'):
            raise RuntimeError('samtools did not create a valid BAM index: {}'.format(index))

    @staticmethod
    def _readVCFContigs(path, compressed=None):
        """Validate VCF headers; allow explicit compression for temporary download paths."""
        path = os.fspath(path)
        if compressed is None:
            compressed = path.endswith('.gz')
        if compressed:
            with open(path, 'rb') as handle:
                header = handle.read(12)
                if len(header) != 12 or header[:3] != b'\x1f\x8b\x08' or not header[3] & 4:
                    raise RuntimeError('VCF must be BGZF-compressed, not ordinary gzip: {}'.format(path))
                extraSize = struct.unpack('<H', header[10:12])[0]
                extra = handle.read(extraSize)
            if len(extra) != extraSize:
                raise RuntimeError('Truncated BGZF header: {}'.format(path))
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
            if blockSize is None or blockSize < 12 + extraSize + 8:
                raise RuntimeError('VCF must use valid BGZF blocks for tabix indexing: {}'.format(path))
            if blockSize > os.path.getsize(path):
                raise RuntimeError('Truncated first BGZF block: {}'.format(path))
        opener = gzip.open if compressed else open
        contigs = {}
        with opener(path, 'rt') as handle:
            if not handle.readline().startswith('##fileformat=VCFv'):
                raise RuntimeError('Input is not a VCF: {}'.format(path))
            for line in handle:
                if line.startswith('##contig=<'):
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
                    length = None
                    if 'length' in fields:
                        try:
                            length = int(fields['length'])
                            if length <= 0:
                                raise ValueError
                        except ValueError as error:
                            raise RuntimeError('Invalid VCF contig length for {}'.format(name)) from error
                    contigs[name] = length
                if line.startswith('#CHROM\t'):
                    columns = line.rstrip('\r\n').split('\t')
                    if columns[:8] != ['#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO']:
                        raise RuntimeError('VCF has an invalid column header: {}'.format(path))
                    if len(columns) != 8:
                        samples = columns[9:]
                        if (len(columns) < 10 or columns[8] != 'FORMAT'
                                or any(not name for name in samples) or len(set(samples)) != len(samples)):
                            raise RuntimeError('VCF has invalid or duplicate sample columns: {}'.format(path))
                    return contigs
                if not line.startswith('##'):
                    break
        raise RuntimeError('VCF has no column header: {}'.format(path))

    def _checkContigs(self, reference, contigs, label, requireAll=False):
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
    def _getReferenceFasta(alignment):
        return alignment.getReferenceFasta() if alignment is not None else None

    @staticmethod
    def _getSequenceDictionaryPath(referenceFasta):
        root, extension = os.path.splitext(referenceFasta)

        if extension.lower() in ('.fa', '.fasta', '.fna'):
            return root + '.dict'

        return referenceFasta + '.dict'

    @classmethod
    def _validateBamFile(cls, bamFile):
        cls._checkNonempty(bamFile, 'BAM')

    @staticmethod
    def _quote(value):
        return shlex.quote(str(value))

    @staticmethod
    def _checkNonempty(path, label):
        if not path or not os.path.isfile(path) or os.path.getsize(path) == 0:
            raise RuntimeError('{} is missing or empty: {}'.format(label, path))


def iterVepRows(path):
    """Yield native VEP tab rows, expanding Extra without discarding its text.

    Accept comment-only metadata and header-only tables, but require the
    standard VEP fields. Both Extra ENSP=... and explicit ENSP columns work.
    Malformed rows are errors rather than silently truncated annotations.
    """
    columns = None
    required = {'Uploaded_variation', 'Location', 'Allele', 'Gene', 'Feature',
                'Feature_type', 'Consequence', 'Protein_position', 'Amino_acids'}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('##'):
                continue
            if line.startswith('#'):
                if columns is not None:
                    raise RuntimeError('Repeated VEP header at line {}.'.format(number))
                columns = line.rstrip('\r\n').lstrip('#').split('\t')
                if len(columns) != len(set(columns)) or not required.issubset(columns):
                    raise RuntimeError('Missing or duplicate required VEP columns.')
                continue
            if columns is None:
                raise RuntimeError('VEP data appears before the column header.')
            values = line.rstrip('\r\n').split('\t')
            if len(values) != len(columns):
                raise RuntimeError('VEP line {} has the wrong number of columns.'.format(number))
            row = dict(zip(columns, values))
            for item in row.get('Extra', '').split(';'):
                if '=' in item:
                    key, value = item.split('=', 1)
                    if key not in row:
                        row[key] = value
            yield row
    if columns is None:
        raise RuntimeError('VEP column header not found.')


def _readSymbolExpression(path, column=''):
    """Read first FPKM per exact gene symbol, matching the reference script.

    The source is isoforms.fpkm_tracking; retaining the first symbol is not
    a gene-level sum. Cuffdiff needs an explicit condition FPKM column.
    """
    import csv
    import math
    values = {}
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames or []
        column = column or ('FPKM' if 'FPKM' in fields else '')
        if len(fields) != len(set(fields)) or 'gene_short_name' not in fields:
            raise RuntimeError('Expression table needs unique columns and gene_short_name.')
        if not column or column not in fields or not (column == 'FPKM' or column.endswith('_FPKM')):
            raise RuntimeError('Select an existing FPKM column (Cuffdiff requires a condition).')
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise RuntimeError('Malformed expression table row.')
            symbol = row['gene_short_name']
            if symbol in ('', '-') or symbol in values:
                continue
            try:
                value = float(row[column])
                if not math.isfinite(value) or value < 0:
                    raise ValueError
            except ValueError as error:
                raise RuntimeError('FPKM must be finite and nonnegative for {}.'.format(symbol)) from error
            values[symbol] = value
    return values


def _getUniParcId(value):
    """Return one unambiguous UniParc accession, never choose among isoforms."""
    accessions = set((value or '').split(','))
    if len(accessions) == 1:
        accession = next(iter(accessions))
        if re.fullmatch(r'UPI[0-9A-F]{10,}', accession):
            return accession
    return None


def _readUniParcFasta(path, requested):
    """Read requested UPI sequences from plain/gzip FASTA, rejecting duplicates."""
    from Bio import SeqIO
    result = {}
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        found = False
        for record in SeqIO.parse(handle, 'fasta'):
            found = True
            matches = set(re.findall(r'\bUPI[0-9A-F]{10,}\b', record.id))
            if len(matches) != 1:
                continue
            accession = next(iter(matches))
            if accession in requested:
                if accession in result:
                    raise RuntimeError('Duplicate UniParc FASTA identifier: {}.'.format(accession))
                sequence = str(record.seq).upper()
                if not sequence:
                    raise RuntimeError('Empty UniParc protein: {}.'.format(accession))
                result[accession] = sequence
        if not found:
            raise RuntimeError('UniParc FASTA has no sequences.')
    return result


def fetchUniParcSequence(accession, cacheDirectory, session, attempts=3):
    """Fetch one exact UPI sequence, with timeout, bounded retries and cache.

    Queries the same UniProt search endpoint as the supplied script. Responses
    containing several FASTA entries or another accession are errors. A 404
    returns None; exhausted network/server failures stop the run rather than
    silently dropping mutations. Successfully validated responses are cached.
    """
    import io
    import time
    import requests
    from Bio import SeqIO

    if _getUniParcId(accession) != accession:
        raise ValueError('Invalid UniParc accession.')
    os.makedirs(cacheDirectory, exist_ok=True)
    path = os.path.join(cacheDirectory, accession + '.fasta')
    if os.path.isfile(path):
        cached = _readUniParcFasta(path, {accession})
        if accession not in cached:
            raise RuntimeError('Cached UniParc file does not contain {}.'.format(accession))
        return cached[accession]
    if attempts < 1:
        raise ValueError('At least one download attempt is required.')
    for attempt in range(attempts):
        try:
            response = session.get('https://rest.uniprot.org/uniparc/search',
                                   params={'query': accession, 'format': 'fasta'}, timeout=(10, 60))
        except requests.RequestException as error:
            if attempt + 1 == attempts:
                raise RuntimeError('UniParc download failed for {}.'.format(accession)) from error
        else:
            try:
                if response.status_code == 404:
                    return None
                if response.status_code == 429 or response.status_code >= 500:
                    if attempt + 1 == attempts:
                        raise RuntimeError('UniParc returned HTTP {} for {}.'.format(response.status_code, accession))
                else:
                    response.raise_for_status()
                    records = list(SeqIO.parse(io.StringIO(response.text), 'fasta'))
                    if not records and not response.text.strip():
                        return None
                    if len(records) != 1 or accession not in re.findall(r'\bUPI[0-9A-F]{10,}\b', records[0].id):
                        raise RuntimeError('UniParc response is not a single matching protein for {}.'.format(accession))
                    sequence = str(records[0].seq).upper()
                    if not sequence:
                        raise RuntimeError('Empty UniParc sequence for {}.'.format(accession))
                    with tempfile.NamedTemporaryFile(mode='w', dir=cacheDirectory, delete=False) as handle:
                        temporary = handle.name
                        handle.write('>{}\n{}\n'.format(accession, sequence))
                    try:
                        os.replace(temporary, path)
                    finally:
                        if os.path.exists(temporary):
                            os.remove(temporary)
                    return sequence
            except requests.HTTPError as error:
                raise RuntimeError('UniParc rejected {} (HTTP {}).'.format(accession, response.status_code)) from error
            finally:
                response.close()
        time.sleep(attempt + 1)


def buildMutationReport(annotations, proteinFasta, directory, expression,
                        expressionColumn='', minimumFPKM=0, peptideLength=25,
                        cacheDirectory=None):
    """Write reference-compatible candidates.csv and detailed mutation records.

    UniParc proteins, exact SYMBOL/gene_short_name matching, first FPKM per
    symbol, zero for absent symbols, and first unique mutant sequence follow
    the supplied script. Windows now contain exactly peptideLength residues,
    centred on the mutation (25 by default). Terminal, ambiguous and mismatched
    effects are reported instead of sliced incorrectly or fabricated.

    A local UniParc FASTA permits offline reproduction. Otherwise use cached
    UniProt REST responses; all used sequences are saved with the result.
    Unlike the original script, failed requests cannot reuse a prior nseq.
    """
    import csv
    import math
    import requests

    if not isinstance(peptideLength, int) or isinstance(peptideLength, bool) or peptideLength < 1 or peptideLength % 2 == 0:
        raise ValueError('Peptide length must be a positive odd integer.')
    if not math.isfinite(minimumFPKM) or minimumFPKM < 0:
        raise ValueError('Minimum FPKM must be finite and nonnegative.')
    if not expression:
        raise ValueError('Cufflinks/Cuffdiff expression is required.')
    values = _readSymbolExpression(expression, expressionColumn)
    os.makedirs(directory, exist_ok=True)
    requested = {_getUniParcId(row.get('UNIPARC')) for row in iterVepRows(annotations)
                 if row.get('SYMBOL') not in (None, '', '-')
                 and 'missense_variant' in row['Consequence'].replace('&', ',').split(',')}
    requested.discard(None)
    if proteinFasta:
        proteins = _readUniParcFasta(proteinFasta, requested)
    else:
        cacheDirectory = cacheDirectory or os.path.join(directory, 'uniparc_cache')
        with requests.Session() as session:
            proteins = {accession: sequence for accession in sorted(requested)
                        for sequence in [fetchUniParcSequence(accession, cacheDirectory, session)]
                        if sequence is not None}
    with open(os.path.join(directory, 'reference_sequences.fasta'), 'w') as handle:
        for accession, sequence in sorted(proteins.items()):
            handle.write('>{}\n{}\n'.format(accession, sequence))
    fields = ['effect_id', 'variant', 'location', 'allele', 'gene_id', 'symbol', 'transcript_id',
              'protein_id', 'consequence', 'protein_position', 'reference_aa', 'mutant_aa',
              'FPKM', 'window_start', 'window_end', 'mutation_offset', 'reference_peptide', 'mutant_peptide']
    aminoAcids = set('ACDEFGHIKLMNPQRSTVWY')
    flank = peptideLength // 2
    with open(os.path.join(directory, 'mutations.tsv'), 'w') as table, \
            open(os.path.join(directory, 'peptides.fasta'), 'w') as fasta, \
            open(os.path.join(directory, 'rejected.tsv'), 'w') as rejected, \
            open(os.path.join(directory, 'candidates.csv'), 'w') as candidates:
        writer = csv.writer(table, delimiter='\t', lineterminator='\n')
        rejectWriter = csv.writer(rejected, delimiter='\t', lineterminator='\n')
        candidateWriter = csv.writer(candidates, lineterminator='\n')
        writer.writerow(fields)
        rejectWriter.writerow(['variant', 'transcript_id', 'protein_id', 'reason'])
        candidateWriter.writerow(['seq', 'symbol', 'FPKM'])
        accepted, seen = 0, set()
        for row in iterVepRows(annotations):
            reason = ''
            consequence, transcript = row['Consequence'], row['Feature']
            protein, symbol = _getUniParcId(row.get('UNIPARC')), row.get('SYMBOL', '')
            change, posText = row['Amino_acids'].split('/'), row['Protein_position']
            if 'missense_variant' not in consequence.replace('&', ',').split(','):
                reason = 'not_missense'
            elif row['Feature_type'] != 'Transcript':
                reason = 'not_transcript'
            elif not protein:
                reason = 'missing_or_ambiguous_uniparc'
            elif not symbol or symbol == '-':
                reason = 'missing_symbol'
            elif len(change) != 2 or any(aa not in aminoAcids for aa in change) or change[0] == change[1]:
                reason = 'unsupported_amino_acid_change'
            elif not posText.isdigit() or int(posText) < 1:
                reason = 'unsupported_protein_position'
            sequence = proteins.get(protein, '')
            position = int(posText) - 1 if posText.isdigit() else -1
            if not reason and not sequence:
                reason = 'missing_uniparc_protein'
            if not reason and (position >= len(sequence) or sequence[position] != change[0]):
                reason = 'reference_residue_mismatch'
            start, end = position - flank, position + flank + 1
            if not reason and (start < 0 or end > len(sequence)):
                reason = 'incomplete_peptide_window'
            reference = sequence[max(0, start):end]
            if not reason and not set(reference).issubset(aminoAcids):
                reason = 'unsupported_residue_in_window'
            fpkm = values.get(symbol, 0.0)
            if not reason and fpkm < minimumFPKM:
                reason = 'below_expression_threshold'
            if reason:
                rejectWriter.writerow([row['Uploaded_variation'], transcript, row.get('UNIPARC', ''), reason])
                continue
            accepted += 1
            effect = 'effect_{:06d}'.format(accepted)
            mutant = reference[:flank] + change[1] + reference[flank + 1:]
            writer.writerow([effect, row['Uploaded_variation'], row['Location'], row['Allele'],
                row['Gene'], symbol, transcript, protein, consequence, position + 1,
                change[0], change[1], fpkm, start + 1, end, flank + 1, reference, mutant])
            fasta.write('>{}_WT\n{}\n>{}_MUT\n{}\n'.format(effect, reference, effect, mutant))
            if mutant not in seen:
                candidateWriter.writerow([mutant, symbol, fpkm])
                seen.add(mutant)
    return accepted


def checkGtfAnnotation(path, contigs, featureTypes=('exon', 'CDS', 'UTR', 'five_prime_utr', 'three_prime_utr'),
                       attribute='transcript_id'):
    """Validate GTF coordinates and feature identifiers against a BAM dictionary.

    Return annotation-only contigs for the caller to report. Cufflinks
    requires transcript-bearing features with transcript_id attributes.
    """
    sharedFeatures = 0
    extraContigs = set()
    pattern = re.compile(r'(?:^|;)\s*' + re.escape(attribute) + r'\s+"[^"\n]+"')
    with open(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.rstrip('\r\n').split('\t')
            if len(fields) != 9:
                raise RuntimeError('GTF line {} must contain nine tab-separated columns.'.format(number))
            try:
                start, end = int(fields[3]), int(fields[4])
                if start < 1 or end < start:
                    raise ValueError
            except ValueError as error:
                raise RuntimeError('Invalid GTF coordinates at line {}.'.format(number)) from error
            name = fields[0]
            if name not in contigs:
                extraContigs.add(name)
                continue
            if end > contigs[name]:
                raise RuntimeError('GTF line {} exceeds BAM sequence {} length.'.format(number, name))
            if fields[2] in featureTypes:
                if not pattern.search(fields[8]):
                    raise RuntimeError('GTF feature lacks {} at line {}.'.format(attribute, number))
                sharedFeatures += 1
    if not sharedFeatures:
        raise RuntimeError('GTF has no transcript features on the BAM reference sequences for the selected feature types.')
    return extraContigs
