# ***************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# ***************************************************************************

import csv
import json
import os
import shlex
from itertools import islice

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.objects import RNASeqExpression
from pwchem.protocols import (
    ProtCufflinks, ProtGATK, ProtImportFastq, ProtImportGenomes,
    ProtImportVCF, ProtPicard, ProtRNASeqAlignment
)
from pwchem.utils.utilsRNA import getFastqStats, openFastq


N_TEST_READS = 1000


class TestCufflinks(BaseTest):
    """Exercise Cufflinks and Cuffdiff with the mouse genomics dataset.

    Two non-overlapping blocks from mouseFastq are prepared as separate BAMs.
    These are technical subsets of one sample, not independent biological
    samples, replicates or a matched tumor-normal pair. Cuffdiff checks verify
    execution, condition order and output files, not differential-expression
    accuracy or significance.

    STAR, Picard, GATK, samtools and Cufflinks/Cuffdiff must be installed in
    the RNA-seq environment. Reference imports and prepared BAMs are reused
    between test cases; each subset follows the reference preprocessing chain.
    """

    SAMPLE_NAMES = ('mouse_subset_A', 'mouse_subset_B')
    CONDITION_LABELS = ('subset_A', 'subset_B')

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('genomics')
        cls._genomes = None
        cls._knownSites = None
        cls._alignments = {}

    def _launchWithOutput(self, protocol, name):
        self.launchProtocol(protocol)
        self.assertTrue(hasattr(protocol, name), 'Missing {}'.format(name))
        return getattr(protocol, name)

    def _assertFile(self, path, allowEmpty=False):
        self.assertTrue(path and os.path.isfile(path), 'Missing file: {}'.format(path))
        if not allowEmpty:
            self.assertGreater(os.path.getsize(path), 0, 'Empty file: {}'.format(path))

    def _getGenomes(self):
        cls = type(self)
        if cls._genomes is None:
            protocol = self.newProtocol(
                ProtImportGenomes, objLabel='Import Cufflinks test reference',
                genomesData=json.dumps([{
                    'scientificName': 'Mus musculus', 'assembly': 'GRCm39', 'release': '116',
                    'fastaFile': self.dataset.getFile('mouseGenomeFile'),
                    'gtfFile': self.dataset.getFile('mouseGtfFile')
                }])
            )
            genomes = self._launchWithOutput(protocol, 'referenceGenomes')
            self.assertEqual(genomes.getSize(), 1)
            cls._genomes = genomes
        return cls._genomes

    def _getKnownSites(self):
        cls = type(self)
        if cls._knownSites is None:
            protocol = self.newProtocol(
                ProtImportVCF, objLabel='Import Cufflinks BQSR known sites',
                vcfsData=json.dumps([{
                    'scientificName': 'Mus musculus', 'assembly': 'GRCm39',
                    'vcfFile': self.dataset.getFile('mouseVcfFile'),
                    'indexFile': self.dataset.getFile('mouseVcfIndex')
                }])
            )
            sites = self._launchWithOutput(protocol, 'outputVCFs')
            self.assertEqual(sites.getSize(), 1)
            cls._knownSites = sites
        return cls._knownSites

    def _createFastqSubset(self, subset):
        filename = self.getOutputPath(self.SAMPLE_NAMES[subset] + '.fastq')
        start = subset * N_TEST_READS * 4
        stop = start + N_TEST_READS * 4
        with openFastq(self.dataset.getFile('mouseFastq')) as source, open(filename, 'w') as target:
            target.writelines(islice(source, start, stop))
        reads, _ = getFastqStats(filename)
        self.assertEqual(reads, N_TEST_READS,
                         'mouseFastq must contain at least {} reads for subset {}.'
                         .format((subset + 1) * N_TEST_READS, subset))
        return filename

    def _getAlignment(self, subset=0):
        cls = type(self)
        if subset in cls._alignments:
            return cls._alignments[subset]

        sample = self.SAMPLE_NAMES[subset]
        importReads = self.newProtocol(
            ProtImportFastq, objLabel='Import ' + sample,
            inputFastq1=self._createFastqSubset(subset), sampleName=sample,
            isPaired=False, runFastqc=False
        )
        reads = self._launchWithOutput(importReads, 'outputFastq')
        star = self.newProtocol(
            ProtRNASeqAlignment, objLabel='STAR ' + sample,
            genomeIndex=0, aligner=ProtRNASeqAlignment.ALIGN_STAR,
            rnaStrandness=0, keepIntermediateFiles=False,
            numberOfThreads=4, numberOfMpi=1
        )
        star.inputFastq.set(reads)
        star.inputGenomes.set(self._getGenomes())
        aligned = self._launchWithOutput(star, 'outputAlignment')
        picard = self.newProtocol(
            ProtPicard, objLabel='Read groups and duplicates ' + sample,
            addReadGroups=True, markDuplicates=True, reorderBam=False
        )
        picard.inputAlignment.set(aligned)
        prepared = self._launchWithOutput(picard, 'outputAlignment')
        gatk = self.newProtocol(
            ProtGATK, objLabel='SplitNCigarReads and BQSR ' + sample,
            splitNCigarReads=True, baseRecalibrator=True, applyBQSR=True,
            vcfIndex=0, numberOfThreads=4, numberOfMpi=1
        )
        gatk.inputAlignment.set(prepared)
        gatk.knownSites.set(self._getKnownSites())
        recalibrated = self._launchWithOutput(gatk, 'outputAlignment')
        reorder = self.newProtocol(
            ProtPicard, objLabel='ReorderSam ' + sample,
            addReadGroups=False, markDuplicates=False, reorderBam=True
        )
        reorder.inputAlignment.set(recalibrated)
        result = self._launchWithOutput(reorder, 'outputAlignment')
        self._assertFile(result.getFileName())
        self._assertFile(result.getIndexFile())
        self.assertEqual(result.getSampleName(), sample)
        self.assertTrue(result.isSorted())
        self.assertTrue(result.isIndexed())
        cls._alignments[subset] = result
        return result

    def _newCufflinks(self, **kwargs):
        options = dict(numberOfThreads=2, numberOfMpi=1)
        options.update(kwargs)
        protocol = self.newProtocol(ProtCufflinks, **options)
        protocol.inputAlignment.set(self._getAlignment())
        return protocol

    def _checkTable(self, path, required):
        self._assertFile(path)
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            columns = reader.fieldnames
            self.assertIsNotNone(columns, 'Missing table header: {}'.format(path))
            self.assertEqual(len(columns), len(set(columns)))
            self.assertTrue(set(required).issubset(columns),
                            'Missing columns in {}: {}'.format(path, set(required) - set(columns)))
            for row in reader:
                self.assertNotIn(None, row, 'Extra table fields in {}'.format(path))
                self.assertTrue(all(value is not None for value in row.values()),
                                'Missing table fields in {}'.format(path))
        return columns

    def _checkOutput(self, protocol):
        self.assertTrue(hasattr(protocol, 'outputExpression'))
        output = protocol.outputExpression
        self.assertIsInstance(output, RNASeqExpression)
        first = protocol.inputAlignment.get()
        self.assertEqual(output.getSampleName(), first.getSampleName())
        self.assertEqual(output.getScientificName(), 'Mus musculus')
        self.assertEqual(output.getAssembly(), 'GRCm39')
        self.assertTrue(os.path.samefile(output.getReferenceFasta(), first.getReferenceFasta()))
        self.assertTrue(os.path.samefile(output.getReferenceGtf(), protocol._getAnnotation()))
        self.assertEqual(output.getOutputDirectory(), protocol._getOutputDir())
        self.assertEqual(os.path.basename(output.getFileName()), 'isoforms.fpkm_tracking')
        self.assertEqual(os.path.basename(output.getGenesFile()), 'genes.fpkm_tracking')

        headers = [self._checkTable(path, ['tracking_id', 'gene_id', 'locus'])
                   for path in (output.getFileName(), output.getGenesFile())]
        if protocol.mode.get() == ProtCufflinks.CUFFLINKS:
            self.assertEqual(output.getProgram(), 'cufflinks')
            self.assertEqual(output.getSampleLabel(), first.getSampleName())
            self.assertIsNone(output.getComparisonSample())
            self.assertIsNone(output.getComparisonLabel())
            self.assertIsNone(output.getIsoformDiffFile())
            self.assertIsNone(output.getGeneDiffFile())
            self._assertFile(output.getTranscriptsFile(), allowEmpty=True)
            for columns in headers:
                self.assertIn('FPKM', columns)
        else:
            self.assertEqual(output.getProgram(), 'cuffdiff')
            self.assertEqual(output.getComparisonSample(), protocol.inputComparison.get().getSampleName())
            self.assertEqual((output.getSampleLabel(), output.getComparisonLabel()), self.CONDITION_LABELS)
            self.assertIsNone(output.getTranscriptsFile())
            for columns in headers:
                self.assertEqual(len([column for column in columns if column.endswith('_FPKM')]), 2)
            required = ['test_id', 'sample_1', 'sample_2', 'status', 'p_value', 'q_value', 'significant']
            for path in (output.getIsoformDiffFile(), output.getGeneDiffFile()):
                self._checkTable(path, required)

        for path in output.getFiles():
            self._assertFile(path, allowEmpty=path == output.getTranscriptsFile())
        self.assertTrue(output.getCommand())
        with open(protocol._getCommandFile()) as handle:
            commands = [shlex.split(line) for line in handle if line.strip()]
        command = next(tokens for tokens in commands if tokens[0] == output.getProgram())
        self.assertEqual(command[command.index('--num-threads') + 1], '2')
        self.assertEqual(command[command.index('--library-type') + 1], 'fr-unstranded')
        return command

    def testCufflinks(self):
        """Assemble one prepared subset using the alignment's annotation."""
        protocol = self._newCufflinks(objLabel='Cufflinks single sample')
        self.launchProtocol(protocol)
        command = self._checkOutput(protocol)
        self.assertIn('-g', command)
        self.assertNotIn('-G', command)
        self.assertNotIn('-L', command)
        self.assertEqual(command[-1], protocol._getAlignmentPath(0))

    def testCuffdiff(self):
        """Compare two disjoint technical subsets without asserting biological differences."""
        protocol = self._newCufflinks(
            objLabel='Cuffdiff technical subsets', mode=ProtCufflinks.CUFFDIFF,
            firstLabel=self.CONDITION_LABELS[0], secondLabel=self.CONDITION_LABELS[1]
        )
        second = self._getAlignment(1)
        self.assertFalse(os.path.samefile(protocol.inputAlignment.get().getFileName(), second.getFileName()))
        protocol.inputComparison.set(second)
        self.launchProtocol(protocol)
        command = self._checkOutput(protocol)
        self.assertEqual(command[command.index('-L') + 1], ','.join(self.CONDITION_LABELS))
        self.assertEqual(command[-3:], [protocol._getAnnotationPath(),
                                        protocol._getAlignmentPath(0), protocol._getAlignmentPath(1)])

    def testAnnotationOverride(self):
        """Use an explicitly supplied GTF rather than the metadata path."""
        annotation = self.dataset.getFile('mouseGtfFile')
        protocol = self._newCufflinks(objLabel='Cufflinks explicit annotation', annotationFile=annotation)
        self.launchProtocol(protocol)
        self._checkOutput(protocol)
        self.assertTrue(os.path.samefile(protocol.outputExpression.getReferenceGtf(), annotation))

    def testInvalidInputs(self):
        """Reject missing comparisons, identical BAMs, duplicate labels and invalid options."""
        protocol = self._newCufflinks(mode=ProtCufflinks.CUFFDIFF)
        self.assertIn('Sample 2 AlignmentFile is required.', protocol._validate())
        protocol.inputComparison.set(protocol.inputAlignment.get())
        self.assertTrue(any('different BAM files' in error for error in protocol._validate()))
        protocol.secondLabel.set(protocol.firstLabel.get())
        self.assertTrue(any('labels must be distinct' in error for error in protocol._validate()))
        protocol.mode.set(ProtCufflinks.CUFFLINKS)
        protocol.libraryType.set(3)
        self.assertIn('Select a valid library type.', protocol._validate())
        protocol.libraryType.set(0)
        protocol.numberOfThreads.set(0)
        self.assertIn('At least one thread is required.', protocol._validate())
        protocol.numberOfThreads.set(2)
        protocol.annotationFile.set(self.getOutputPath('missing_annotation.gtf'))
        self.assertTrue(any('Annotation GTF' in error for error in protocol._validate()))
