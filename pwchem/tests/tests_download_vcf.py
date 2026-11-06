# **************************************************************************
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
# **************************************************************************

import csv
import os
import tempfile
import unittest
from unittest.mock import Mock, patch

from Bio import bgzf

from pyworkflow.tests import BaseTest, setupTestProject

from pwchem.protocols import ProtDownloadVCF
from pwchem.utils import assertHandle
from pwchem.utils.utilsRNA import assertOutputExists


class TestVCFRecordValidation(unittest.TestCase):
    """Small offline fixtures for malformed downloaded/reused records."""

    HEADER = ('##fileformat=VCFv4.2\n'
              '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')

    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.filename = os.path.join(self.directory.name, 'variants.vcf.gz')

    def writeVCF(self, records):
        with bgzf.BgzfWriter(self.filename, 'wb') as handle:
            handle.write((self.HEADER + records).encode('utf-8'))

    def testValidAlleleForms(self):
        records = ''.join('3\t16894810\trs1\tA\t{}\t.\t.\t.\n'.format(alt)
                          for alt in ('T', 'AC', 'T,G', '.', '*', '<DEL>', '<NON_REF>',
                                      'A[2:123[', 'A]2:123]', ']2:123]A', '[2:123[A',
                                      '.A', 'A.'))
        self.writeVCF(records)
        self.assertEqual(ProtDownloadVCF._validateVCFRecords(self.filename), set())

    def testRejectSymbolicReferenceWithLocation(self):
        self.writeVCF('3\t16894810\trs1553632783\t<W>\tTAA\t.\t.\t.\n')
        with self.assertRaisesRegex(RuntimeError,
                                    r'line 3 .*3:16894810; ID=rs1553632783.*invalid REF'):
            ProtDownloadVCF._validateVCFRecords(self.filename)

    def testRejectMalformedFields(self):
        for record in ('3\t1\trs1\tA\tT\t.\t.',
                       '3\tx\trs1\tA\tT\t.\t.\t.',
                       '3\t1\trs1\tR\tT\t.\t.\t.',
                       '3\t1\trs1\tA\t<W\t.\t.\t.',
                       '3\t1\trs1\tA\t\t.\t.\t.',
                       '3\t1\trs1\tA\t.,T\t.\t.\t.'):
            with self.subTest(record=record):
                self.writeVCF(record + '\n')
                with self.assertRaisesRegex(RuntimeError, 'Malformed VCF'):
                    ProtDownloadVCF._validateVCFRecords(self.filename)

    def testRejectTruncatedGzip(self):
        self.writeVCF('3\t1\trs1\tA\tT\t.\t.\t.\n')
        with open(self.filename, 'r+b') as handle:
            handle.truncate(os.path.getsize(self.filename) - 8)
        with self.assertRaisesRegex(RuntimeError, 'Invalid compressed VCF'):
            ProtDownloadVCF._validateVCFRecords(self.filename)

    def repairProtocol(self):
        protocol = Mock()
        protocol._validateVCFRecords = ProtDownloadVCF._validateVCFRecords
        protocol._validateVCFHeader = ProtDownloadVCF._validateVCFHeader
        protocol._repairVCF = lambda path: ProtDownloadVCF._repairVCF(protocol, path)
        protocol._getLocalVCFIndex.return_value = None
        return protocol

    def testRepairPreservesOriginalAndReportsExclusions(self):
        invalid = '3\t16894810\trs1553632783\t<W>\tTAA\t.\t.\t.\n'
        valid = '3\t16894810\trs2\tT\tA\t.\t.\t.\n'
        self.writeVCF(invalid + valid)
        with open(self.filename, 'rb') as handle:
            original = handle.read()
        protocol = self.repairProtocol()
        self.assertTrue(protocol._repairVCF(self.filename))
        with open(self.filename + '.original', 'rb') as handle:
            self.assertEqual(handle.read(), original)
        with bgzf.BgzfReader(self.filename, 'rb') as handle:
            self.assertEqual(b''.join(handle), (self.HEADER + valid).encode())
        with open(self.filename + '.excluded.tsv') as handle:
            report = list(csv.DictReader(handle, delimiter='\t'))
        self.assertIn('invalid REF', report[0]['reason'])
        self.assertEqual(report[0]['original_record'], invalid.rstrip('\n'))
        self.assertFalse(protocol._repairVCF(self.filename))
        self.assertEqual(ProtDownloadVCF._validateVCFRecords(self.filename), set())

    def testAllInvalidRecordsLeaveSourceUnchanged(self):
        self.writeVCF('3\t1\trs1\t<W>\tT\t.\t.\t.\n')
        with open(self.filename, 'rb') as handle:
            original = handle.read()
        with self.assertRaisesRegex(RuntimeError, 'no valid records'):
            self.repairProtocol()._repairVCF(self.filename)
        with open(self.filename, 'rb') as handle:
            self.assertEqual(handle.read(), original)
        self.assertFalse(os.path.exists(self.filename + '.original'))
        self.assertFalse(os.path.exists(self.filename + '.repair.tmp.vcf.gz'))

    def testRepairedVCFUsesLocalIndex(self):
        self.writeVCF('3\t1\trs1\t<W>\tT\t.\t.\t.\n'
                      '3\t2\trs2\tA\tT\t.\t.\t.\n')
        protocol = self.repairProtocol()
        ProtDownloadVCF._ensureVCFIndex(protocol, self.filename, 'https://example/vcf.gz')
        protocol._createLocalVCFIndex.assert_called_once_with(self.filename)
        protocol._downloadVCFIndex.assert_not_called()
        protocol._createLocalVCFIndex.reset_mock()
        # A resumed run must not reuse the provider's unfiltered index either.
        ProtDownloadVCF._ensureVCFIndex(protocol, self.filename, 'https://example/vcf.gz')
        protocol._createLocalVCFIndex.assert_called_once_with(self.filename)
        protocol._downloadVCFIndex.assert_not_called()

    def testRepairDoesNotPublishTruncatedSource(self):
        self.writeVCF('3\t1\trs1\t<W>\tT\t.\t.\t.\n'
                      '3\t2\trs2\tA\tT\t.\t.\t.\n')
        with open(self.filename, 'r+b') as handle:
            handle.truncate(os.path.getsize(self.filename) - 8)
        with open(self.filename, 'rb') as handle:
            original = handle.read()
        with self.assertRaisesRegex(RuntimeError, 'Invalid compressed VCF'):
            self.repairProtocol()._repairVCF(self.filename)
        with open(self.filename, 'rb') as handle:
            self.assertEqual(handle.read(), original)
        self.assertFalse(os.path.exists(self.filename + '.original'))
        self.assertFalse(os.path.exists(self.filename + '.excluded.tsv'))

    def testFailedReportPublicationRestoresSource(self):
        self.writeVCF('3\t1\trs1\t<W>\tT\t.\t.\t.\n'
                      '3\t2\trs2\tA\tT\t.\t.\t.\n')
        with open(self.filename, 'rb') as handle:
            original = handle.read()
        replace = os.replace

        def failReport(source, destination):
            if destination == self.filename + '.excluded.tsv':
                raise OSError('report publication failed')
            return replace(source, destination)

        with patch('os.replace', side_effect=failReport):
            with self.assertRaisesRegex(OSError, 'report publication failed'):
                self.repairProtocol()._repairVCF(self.filename)
        with open(self.filename, 'rb') as handle:
            self.assertEqual(handle.read(), original)
        self.assertFalse(os.path.exists(self.filename + '.excluded.tsv'))


class TestDownloadVCF(BaseTest):
    """Test known-variant VCF downloads from Ensembl, NCBI and EVA."""

    @classmethod
    def setUpClass(cls):
        setupTestProject(cls)

    def testEnsemblVCF(self):
        """Download a known-variant VCF from Ensembl."""
        print("\nKnown variants VCF: Ensembl")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_ENSEMBL,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Saccharomyces cerevisiae',
            assemblies='Latest',
            releases='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        self.assertVcfFile(
            prot,
            vcf
        )

    def testNcbiHumanVCF(self):
        """Download the human dbSNP VCF from NCBI."""
        print("\nKnown variants VCF: NCBI dbSNP")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Homo sapiens',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        self.assertVcfFile(
            prot,
            vcf
        )

    def testNcbiEvaFallback(self):
        """
        Check the automatic NCBI -> EVA fallback for a non-human species.

        When NCBI is selected for a non-human species, dbSNP does not
        provide the current VCF and the protocol delegates the download
        to EVA.
        """
        print("\nKnown variants VCF: NCBI -> EVA fallback")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Oryza sativa',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            1,
            cwd=prot.getWorkingDir()
        )

        vcf = output.getFirstItem()

        self.assertVcfFile(
            prot,
            vcf
        )

    def testMultipleEvaFallbackVCFs(self):
        """
        Download VCFs for multiple non-human species through EVA.
        """
        print("\nKnown variants VCF: multiple NCBI -> EVA fallbacks")

        prot = self.newProtocol(
            ProtDownloadVCF,
            source=ProtDownloadVCF.SOURCE_NCBI,
            speciesSelection=ProtDownloadVCF.SPECIES_CUSTOM,
            customSpecies='Bos taurus;Oryza sativa',
            ncbiAssemblies='Latest'
        )

        self.launchProtocol(prot)

        output = getattr(
            prot,
            'outputVCFs',
            None
        )

        assertOutputExists(
            self,
            prot,
            output
        )

        assertHandle(
            self.assertEqual,
            output.getSize(),
            2,
            cwd=prot.getWorkingDir()
        )

        for vcf in output:
            self.assertVcfFile(
                prot,
                vcf
            )


    def assertVcfFile(self, protocol, vcf):
        """Check that the downloaded VCF file exists and is not empty."""
        fileName = vcf.getFileName()

        assertHandle(
            self.assertTrue,
            bool(fileName),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            os.path.exists(fileName),
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertGreater,
            os.path.getsize(fileName),
            0,
            cwd=protocol.getWorkingDir()
        )

        assertHandle(
            self.assertTrue,
            fileName.endswith(('.vcf', '.vcf.gz')),
            cwd=protocol.getWorkingDir()
        )
