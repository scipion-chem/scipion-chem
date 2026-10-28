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

import gzip
import json
import os
import shlex
from itertools import islice

from pyworkflow.tests import BaseTest, DataSet, setupTestProject

from pwchem.protocols import (
    ProtGATK, ProtGATKMutect2, ProtImportFastq, ProtImportGenomes,
    ProtImportVCF, ProtPicard, ProtRNASeqAlignment
)
from pwchem.utils.utilsRNA import getFastqStats, openFastq


N_TEST_READS = 1000


class TestGATKMutect2(BaseTest):
    """Test Mutect2 with the mouse RNA-seq files in the genomics dataset.

    A reduced single-end FASTQ is aligned with STAR, processed with Picard
    and prepared with SplitNCigarReads and BQSR. The prepared BAM is reused
    across calling tests to avoid rebuilding the STAR index for each case.
    The dataset VCF is used only as BQSR known sites, not as a germline
    population resource or a panel of normals.

    This dataset does not identify matched tumor-normal samples or expected
    somatic variants. Calling tests therefore exercise tumor-only mode and
    allow valid indexed VCFs without variant records. They check workflow
    integration, output metadata and variant selection, not somatic accuracy.
    Tumor-normal calling requires an additional paired dataset.
    """

    SAMPLE_NAME = 'mouse_mutect2_test'

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)
        cls.dataset = DataSet.getDataSet('genomics')
        cls._inputAlignment = None

    def _launchWithOutput(self, protocol, outputName):
        self.launchProtocol(protocol)
        self.assertTrue(hasattr(protocol, outputName), 'Missing {}'.format(outputName))
        return getattr(protocol, outputName)

    def _getInputAlignment(self):
        """Prepare the dataset BAM once and retain its Scipion provenance."""
        cls = type(self)
        if cls._inputAlignment is not None:
            return cls._inputAlignment

        fastq = self.getOutputPath('mouse_mutect2_test.fastq')
        with openFastq(self.dataset.getFile('mouseFastq')) as source, open(fastq, 'w') as target:
            target.writelines(islice(source, N_TEST_READS * 4))
        numReads, _ = getFastqStats(fastq)
        self.assertGreater(numReads, 0, 'Reduced FASTQ is empty.')
        self.assertLessEqual(numReads, N_TEST_READS)

        importFastq = self.newProtocol(
            ProtImportFastq, objLabel='Import Mutect2 test FASTQ', inputFastq1=fastq,
            sampleName=self.SAMPLE_NAME, isPaired=False, runFastqc=False
        )
        reads = self._launchWithOutput(importFastq, 'outputFastq')
        importGenome = self.newProtocol(
            ProtImportGenomes, objLabel='Import Mutect2 test genome',
            genomesData=json.dumps([{
                'scientificName': 'Mus musculus', 'assembly': 'GRCm39', 'release': '116',
                'fastaFile': self.dataset.getFile('mouseGenomeFile'),
                'gtfFile': self.dataset.getFile('mouseGtfFile')
            }])
        )
        genomes = self._launchWithOutput(importGenome, 'referenceGenomes')
        alignment = self.newProtocol(
            ProtRNASeqAlignment, objLabel='STAR alignment for Mutect2',
            genomeIndex=0, aligner=ProtRNASeqAlignment.ALIGN_STAR, rnaStrandness=0,
            keepIntermediateFiles=False, numberOfThreads=4, numberOfMpi=1
        )
        alignment.inputFastq.set(reads)
        alignment.inputGenomes.set(genomes)
        starBam = self._launchWithOutput(alignment, 'outputAlignment')
        picard = self.newProtocol(
            ProtPicard, objLabel='Picard preparation for Mutect2',
            addReadGroups=True, markDuplicates=True, reorderBam=False
        )
        picard.inputAlignment.set(starBam)
        picardBam = self._launchWithOutput(picard, 'outputAlignment')
        importSites = self.newProtocol(
            ProtImportVCF, objLabel='Import Mutect2 BQSR known sites',
            vcfsData=json.dumps([{
                'scientificName': 'Mus musculus', 'assembly': 'GRCm39',
                'vcfFile': self.dataset.getFile('mouseVcfFile'),
                'indexFile': self.dataset.getFile('mouseVcfIndex')
            }])
        )
        knownSites = self._launchWithOutput(importSites, 'outputVCFs')
        self.assertEqual(knownSites.getSize(), 1)
        gatk = self.newProtocol(
            ProtGATK, objLabel='RNA-seq BAM preparation for Mutect2',
            splitNCigarReads=True, baseRecalibrator=True, applyBQSR=True,
            vcfIndex=0, numberOfThreads=4, numberOfMpi=1
        )
        gatk.inputAlignment.set(picardBam)
        gatk.knownSites.set(knownSites)
        prepared = self._launchWithOutput(gatk, 'outputAlignment')
        self._assertNonemptyFile(prepared.getFileName())
        self._assertNonemptyFile(prepared.getIndexFile())
        cls._inputAlignment = prepared
        return prepared

    def _newMutect2(self, **kwargs):
        options = dict(numberOfThreads=2, numberOfMpi=1)
        options.update(kwargs)
        protocol = self.newProtocol(ProtGATKMutect2, **options)
        protocol.inputTumor.set(self._getInputAlignment())
        return protocol

    def _assertNonemptyFile(self, path):
        self.assertTrue(path and os.path.isfile(path), 'Missing file: {}'.format(path))
        self.assertGreater(os.path.getsize(path), 0, 'Empty file: {}'.format(path))

    def _checkOutputs(self, protocol):
        expectedSamples = {self.SAMPLE_NAME}
        alignment = protocol.inputTumor.get()
        records = {}
        for name in ('outputUnfilteredVCF', 'outputVCF'):
            self.assertTrue(hasattr(protocol, name), 'Missing {}'.format(name))
            output = getattr(protocol, name)
            self._assertNonemptyFile(output.getFileName())
            self._assertNonemptyFile(output.getIndexFile())
            self.assertEqual(output.getIndexFile(), output.getFileName() + '.tbi')
            with gzip.open(output.getIndexFile(), 'rb') as handle:
                self.assertEqual(handle.read(4), b'TBI\x01')
            self.assertTrue(output.isCompressed())
            self.assertEqual(output.getTumorSample(), self.SAMPLE_NAME)
            self.assertEqual(output.getNormalSample(), None)
            self.assertEqual(output.getScientificName(), alignment.getReferenceScientificName())
            self.assertEqual(output.getAssembly(), alignment.getReferenceAssembly())
            self.assertEqual(output.getSource(), 'GATK')
            self.assertEqual(output.getDatabase(), 'Mutect2')
            self.assertEqual(output.getVariantType(),
                             ['somatic candidates', 'SNP', 'INDEL'][protocol.variantType.get()])
            self.assertEqual(output.getFilteringApplied(),
                             name == 'outputVCF' and protocol.filterCalls.get())
            self.assertTrue(os.path.samefile(output.getReferenceFasta(), alignment.getReferenceFasta()))
            header = None
            records[name] = []
            with gzip.open(output.getFileName(), 'rt') as handle:
                for line in handle:
                    if line.startswith('#CHROM\t'):
                        header = line.rstrip().split('\t')
                    elif not line.startswith('#'):
                        record = line.rstrip().split('\t')
                        self.assertGreaterEqual(len(record), 10)
                        self.assertEqual(len(record), 9 + len(expectedSamples))
                        records[name].append(record)
            self.assertIsNotNone(header, 'VCF sample header is missing.')
            self.assertEqual(header[9:], sorted(expectedSamples))
        if protocol.filterCalls.get():
            self.assertNotEqual(protocol.outputVCF.getFileName(), protocol.outputUnfilteredVCF.getFileName())
            self._assertNonemptyFile(protocol._getStatsFile())
            self.assertEqual([record[:5] for record in records['outputVCF']],
                             [record[:5] for record in records['outputUnfilteredVCF']])
        else:
            self.assertEqual(protocol.outputVCF.getFileName(), protocol.outputUnfilteredVCF.getFileName())
        return records

    def _readCommands(self, protocol):
        with open(protocol._getCommandFile()) as handle:
            return [shlex.split(line) for line in handle if line.strip()]

    def testTumorOnly(self):
        """Run the reference-script calling options on the prepared dataset BAM."""
        protocol = self._newMutect2(objLabel='Mutect2 tumor only')
        self.assertFalse(protocol.filterCalls.get())
        self.launchProtocol(protocol)
        self._checkOutputs(protocol)
        commands = self._readCommands(protocol)
        command = next(tokens for tokens in commands if tokens[:2] == ['gatk', 'Mutect2'])
        self.assertEqual(command.count('-I'), 1)
        self.assertNotIn('-normal', command)
        self.assertIn('--dont-use-soft-clipped-bases', command)
        self.assertIn('MateOnSameContigOrNoMappedMateReadFilter', command)
        self.assertFalse(any(tokens[:2] == ['gatk', 'FilterMutectCalls'] for tokens in commands))

    def testFilteredSNVs(self):
        """Run optional filtering and select SNP records, including non-PASS calls."""
        protocol = self._newMutect2(objLabel='Mutect2 filtered SNVs', filterCalls=True, variantType=1)
        self.launchProtocol(protocol)
        records = self._checkOutputs(protocol)
        for variants in records.values():
            self.assertTrue(all(len(record[3]) == 1 and all(len(alt) == 1 for alt in record[4].split(','))
                                for record in variants))
        commands = self._readCommands(protocol)
        tools = [tokens[1] for tokens in commands if tokens[0] == 'gatk']
        self.assertLess(tools.index('Mutect2'), tools.index('FilterMutectCalls'))
        self.assertLess(tools.index('FilterMutectCalls'), tools.index('SelectVariants'))

    def testIndelsOnly(self):
        """Select indels; allow an indexed VCF without variant records."""
        protocol = self._newMutect2(objLabel='Mutect2 indels only', variantType=2,
                                    disableMateOnSameContigFilter=False, dontUseSoftClippedBases=False)
        self.launchProtocol(protocol)
        records = self._checkOutputs(protocol)
        for variants in records.values():
            self.assertTrue(all(all(len(alt) != len(record[3]) for alt in record[4].split(','))
                                for record in variants))
        commands = self._readCommands(protocol)
        command = next(tokens for tokens in commands if tokens[:2] == ['gatk', 'Mutect2'])
        self.assertNotIn('--disable-read-filter', command)
        self.assertNotIn('--dont-use-soft-clipped-bases', command)
        selection = next(tokens for tokens in commands if tokens[:2] == ['gatk', 'SelectVariants'])
        self.assertEqual(selection[selection.index('--select-type-to-include') + 1], 'INDEL')

    def testInvalidInputs(self):
        """Reject missing normals, a BAM used as its own normal and invalid options."""
        protocol = self._newMutect2(useNormal=True)
        self.assertIn('A normal AlignmentFile is required.', protocol._validate())
        protocol.inputNormal.set(protocol.inputTumor.get())
        self.assertTrue(any('different BAM files' in error for error in protocol._validate()))
        protocol.useNormal.set(False)
        protocol.numberOfThreads.set(0)
        self.assertIn('At least one thread is required.', protocol._validate())
        protocol.numberOfThreads.set(2)
        protocol.variantType.set(3)
        self.assertIn('Select a valid variant type.', protocol._validate())
        protocol.variantType.set(0)
        protocol.useGermlineResource.set(True)
        self.assertTrue(any('germline resource:' in error for error in protocol._validate()))
        protocol.useGermlineResource.set(False)
        protocol.usePanelOfNormals.set(True)
        self.assertTrue(any('pon resource:' in error for error in protocol._validate()))
