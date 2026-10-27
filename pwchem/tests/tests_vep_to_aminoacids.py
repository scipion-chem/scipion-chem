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
import os
import tempfile
import unittest
from unittest.mock import Mock, patch

import requests

from pwchem.objects import RNASeqExpression, VariantAnnotation
from pwchem.protocols import ProtVEPToAminoacids
from pwchem.utils.utilsRNA import buildMutationReport, fetchUniParcSequence


HEADER = ('#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t'
          'Consequence\tProtein_position\tAmino_acids\tExtra\n')
PROTEIN = 'ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY'
ACCESSION = 'UPI0000000001'


class TestVEPToAminoacids(unittest.TestCase):
    """Check 25-aa windows and the supplied script's CSV/expression semantics.

    Conversion runs on actual local FASTA and VEP/Cufflinks tables. UniProt
    responses are mocked only for HTTP tests; no live service is required.
    Expected peptide sequences and coordinates are explicit. Gene-symbol
    matching deliberately tests the first isoform FPKM rather than a sum.
    """

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(prefix='test_peptides_')
        self.addCleanup(self.work.cleanup)
        self.vep = os.path.join(self.work.name, 'vep.tsv')
        self.fasta = os.path.join(self.work.name, 'uniparc.fasta')
        self.expression = os.path.join(self.work.name, 'isoforms.fpkm_tracking')
        with open(self.fasta, 'w') as handle:
            handle.write('>{} UniParc protein\n{}\n'.format(ACCESSION, PROTEIN))
        with open(self.expression, 'w') as handle:
            handle.write('tracking_id\tgene_short_name\tFPKM\nT1\tGene1\t2.5\n')
        self._writeRows([('v1', 'T1', '20', 'Y/A', ACCESSION, 'Gene1')])

    def _writeRows(self, rows, consequence='missense_variant'):
        with open(self.vep, 'w') as handle:
            handle.write(HEADER)
            for variant, transcript, position, change, protein, symbol in rows:
                handle.write('{}\t1:10\tC\tG1\t{}\tTranscript\t{}\t{}\t{}\tUNIPARC={};SYMBOL={}\n'.format(
                    variant, transcript, consequence, position, change, protein, symbol))

    def _run(self, **kwargs):
        directory = os.path.join(self.work.name, 'report')
        kwargs.setdefault('expression', self.expression)
        buildMutationReport(self.vep, self.fasta, directory, **kwargs)
        tables = []
        for filename, delimiter in [('mutations.tsv', '\t'), ('rejected.tsv', '\t'), ('candidates.csv', ',')]:
            with open(os.path.join(directory, filename)) as handle:
                tables.append(list(csv.DictReader(handle, delimiter=delimiter)))
        return tables

    def testExactly25ResiduesAndCentralMutation(self):
        rows, rejected, candidates = self._run()
        self.assertEqual(rejected, [])
        row = rows[0]
        self.assertEqual(row['reference_peptide'], 'IKLMNPQRSTVWYACDEFGHIKLMN')
        self.assertEqual(row['mutant_peptide'], 'IKLMNPQRSTVWAACDEFGHIKLMN')
        self.assertEqual(len(row['mutant_peptide']), 25)
        self.assertEqual(row['mutant_peptide'][12], 'A')
        self.assertEqual([i for i, pair in enumerate(zip(row['reference_peptide'], row['mutant_peptide']))
                          if pair[0] != pair[1]], [12])
        self.assertEqual((row['window_start'], row['window_end'], row['mutation_offset']), ('8', '32', '13'))
        self.assertEqual(candidates, [{'seq': 'IKLMNPQRSTVWAACDEFGHIKLMN', 'symbol': 'Gene1', 'FPKM': '2.5'}])

    def testRejectBothTerminiInsteadOfShortOrPaddedPeptides(self):
        self._writeRows([('v1', 'T1', '1', 'A/Y', ACCESSION, 'Gene1'),
                         ('v2', 'T1', '40', 'Y/A', ACCESSION, 'Gene1')])
        rows, rejected, candidates = self._run()
        self.assertFalse(rows)
        self.assertFalse(candidates)
        self.assertEqual([row['reason'] for row in rejected], ['incomplete_peptide_window'] * 2)

    def testFirstAndLastCompleteWindows(self):
        self._writeRows([('v1', 'T1', '13', 'P/Y', ACCESSION, 'Gene1'),
                         ('v2', 'T1', '28', 'I/A', ACCESSION, 'Gene1')])
        rows, rejected, _ = self._run()
        self.assertFalse(rejected)
        self.assertEqual([(row['window_start'], row['window_end']) for row in rows], [('1', '25'), ('16', '40')])
        self.assertTrue(all(len(row['mutant_peptide']) == 25 for row in rows))

    def testFirstSymbolValueRatherThanIsoformSum(self):
        with open(self.expression, 'a') as handle:
            handle.write('T2\tGene1\t100\n')
        rows, _, candidates = self._run()
        self.assertEqual(rows[0]['FPKM'], '2.5')
        self.assertEqual(candidates[0]['FPKM'], '2.5')

    def testSymbolJoinRatherThanTranscriptJoin(self):
        self._writeRows([('v1', 'different_transcript', '20', 'Y/A', ACCESSION, 'Gene1')])
        rows, _, _ = self._run()
        self.assertEqual(rows[0]['FPKM'], '2.5')

    def testUnmatchedSymbolHasZeroFPKM(self):
        self._writeRows([('v1', 'T1', '20', 'Y/A', ACCESSION, 'UnknownGene')])
        rows, rejected, candidates = self._run()
        self.assertFalse(rejected)
        self.assertEqual(rows[0]['FPKM'], '0.0')
        self.assertEqual(candidates[0]['FPKM'], '0.0')

    def testUniqueSequencesInCSVAndAllEffectsInTSV(self):
        self._writeRows([('v1', 'T1', '20', 'Y/A', ACCESSION, 'Gene1'),
                         ('v2', 'T2', '20', 'Y/A', ACCESSION, 'OtherGene')])
        rows, rejected, candidates = self._run()
        self.assertFalse(rejected)
        self.assertEqual(len(rows), 2)
        self.assertEqual(len(candidates), 1)
        self.assertEqual(candidates[0]['symbol'], 'Gene1')

    def testMissingAndAmbiguousUniParcOrSymbol(self):
        self._writeRows([('v1', 'T1', '20', 'Y/A', '', 'Gene1'),
                         ('v2', 'T1', '20', 'Y/A', ACCESSION + ',UPI0000000002', 'Gene1'),
                         ('v3', 'T1', '20', 'Y/A', ACCESSION, '')])
        rows, rejected, _ = self._run()
        self.assertFalse(rows)
        self.assertEqual([row['reason'] for row in rejected],
                         ['missing_or_ambiguous_uniparc', 'missing_or_ambiguous_uniparc', 'missing_symbol'])

    def testReferenceMismatchAndUnsupportedChanges(self):
        self._writeRows([('v1', 'T1', '20', 'A/Y', ACCESSION, 'Gene1'),
                         ('v2', 'T1', '20-21', 'Y/A', ACCESSION, 'Gene1'),
                         ('v3', 'T1', '20', 'YA/AA', ACCESSION, 'Gene1')])
        rows, rejected, _ = self._run()
        self.assertFalse(rows)
        self.assertEqual([row['reason'] for row in rejected],
                         ['reference_residue_mismatch', 'unsupported_protein_position', 'unsupported_amino_acid_change'])

    def testCuffdiffRequiresCondition(self):
        with open(self.expression, 'w') as handle:
            handle.write('tracking_id\tgene_short_name\ttumor_FPKM\tnormal_FPKM\nT1\tGene1\t4\t1\n')
        with self.assertRaisesRegex(RuntimeError, 'condition'):
            self._run()
        rows, _, _ = self._run(expressionColumn='normal_FPKM')
        self.assertEqual(rows[0]['FPKM'], '1.0')

    def testFPKMThreshold(self):
        rows, rejected, candidates = self._run(minimumFPKM=3)
        self.assertFalse(rows)
        self.assertFalse(candidates)
        self.assertEqual(rejected[0]['reason'], 'below_expression_threshold')

    def testNonMissenseAndHeaderOnly(self):
        self._writeRows([('v1', 'T1', '20', 'Y', ACCESSION, 'Gene1')], consequence='synonymous_variant')
        self.assertEqual(self._run()[1][0]['reason'], 'not_missense')
        self._writeRows([])
        self.assertEqual(self._run(), [[], [], []])

    def testInvalidLengthsAndMissingExpression(self):
        for length in (0, -1, 24):
            with self.subTest(length=length), self.assertRaises(ValueError):
                self._run(peptideLength=length)
        with self.assertRaises(ValueError):
            self._run(expression=None)

    def testProtocolOutputsAndSequenceSnapshot(self):
        protocol = ProtVEPToAminoacids(workingDir=os.path.join(self.work.name, 'protocol'))
        protocol.inputAnnotations.set(VariantAnnotation(filename=self.vep, assembly='GRCm39'))
        protocol.inputExpression.set(RNASeqExpression(filename=self.expression, assembly='GRCm39'))
        protocol.proteinFasta.set(self.fasta)
        self.assertEqual(protocol.peptideLength.get(), 25)
        self.assertEqual(protocol._validate(), [])
        protocol.convertStep()
        outputs = {}
        with patch.object(protocol, '_defineOutputs', side_effect=lambda **kw: outputs.update(kw)), \
                patch.object(protocol, '_defineSourceRelation'):
            protocol.createOutputStep()
        report = outputs['outputMutations']
        self.assertTrue(report.getFileName().endswith('candidates.csv'))
        self.assertEqual(report.getPeptideLength(), 25)
        self.assertEqual(report.getFlankSize(), 12)
        self.assertTrue(all(os.path.isfile(path) for path in report.getFiles()))
        with open(report.getProteinFasta()) as handle:
            self.assertIn(PROTEIN, handle.read())

    def testNetworkSequenceFetchedOnceAndCacheReused(self):
        cache = os.path.join(self.work.name, 'cache')
        session = Mock()
        response = Mock(status_code=200, text='>{}\n{}\n'.format(ACCESSION, PROTEIN))
        session.get.return_value = response
        for _ in range(2):
            self.assertEqual(fetchUniParcSequence(ACCESSION, cache, session), PROTEIN)
        session.get.assert_called_once()
        self.assertEqual(session.get.call_args.kwargs['timeout'], (10, 60))
        response.close.assert_called_once()

    def testRetryConnectionFailure(self):
        session = Mock()
        session.get.side_effect = [requests.ConnectionError('offline'),
            Mock(status_code=200, text='>{}\n{}\n'.format(ACCESSION, PROTEIN))]
        with patch('time.sleep'):
            result = fetchUniParcSequence(ACCESSION, os.path.join(self.work.name, 'cache'), session)
        self.assertEqual(result, PROTEIN)
        self.assertEqual(session.get.call_count, 2)

    def testDownloadFailureStopsWithoutStalePeptide(self):
        session = Mock()
        session.get.side_effect = requests.Timeout('timed out')
        cache = os.path.join(self.work.name, 'cache')
        with patch('time.sleep'), self.assertRaisesRegex(RuntimeError, 'download failed'):
            fetchUniParcSequence(ACCESSION, cache, session)
        self.assertEqual(session.get.call_count, 3)

    def testRejectWrongOrMultipleHTTPRecords(self):
        for fasta in ('>UPI0000000002\nAAAA\n',
                      '>{}\nAAAA\n>UPI0000000002\nAAAA\n'.format(ACCESSION)):
            session = Mock()
            session.get.return_value = Mock(status_code=200, text=fasta)
            cache = os.path.join(self.work.name, 'cache')
            with self.subTest(fasta=fasta), self.assertRaises(RuntimeError):
                fetchUniParcSequence(ACCESSION, cache, session)

    def testLocalSourceRemainsUnchanged(self):
        with open(self.fasta) as handle:
            original = handle.read()
        self._run()
        with open(self.fasta) as handle:
            self.assertEqual(handle.read(), original)

    def testRemoteConversionSnapshotAndCacheReuse(self):
        cache = os.path.join(self.work.name, 'remote_cache')
        directory = os.path.join(self.work.name, 'remote_report')
        self._writeRows([('v1', 'T1', '20', 'Y/A', ACCESSION, 'Gene1'),
                         ('v2', 'T2', '20', 'Y/A', ACCESSION, 'Gene1')])
        with patch('requests.Session') as factory:
            session = factory.return_value.__enter__.return_value
            session.get.return_value = Mock(status_code=200, text='>{}\n{}\n'.format(ACCESSION, PROTEIN))
            for _ in range(2):
                buildMutationReport(self.vep, None, directory, expression=self.expression, cacheDirectory=cache)
            session.get.assert_called_once()
        with open(os.path.join(directory, 'candidates.csv')) as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(len(rows), 1)
        self.assertEqual(len(rows[0]['seq']), 25)
        with open(os.path.join(directory, 'reference_sequences.fasta')) as handle:
            self.assertIn(PROTEIN, handle.read())

    def testExhaustedServerFailuresAreNotMissingExpression(self):
        session = Mock()
        response = Mock(status_code=503)
        session.get.return_value = response
        cache = os.path.join(self.work.name, 'cache')
        with patch('time.sleep'), self.assertRaisesRegex(RuntimeError, 'HTTP 503'):
            fetchUniParcSequence(ACCESSION, cache, session)
        self.assertEqual(session.get.call_count, 3)
        self.assertEqual(response.close.call_count, 3)
