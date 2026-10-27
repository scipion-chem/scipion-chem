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

import os
import shutil
import subprocess
import tempfile
import unittest
from unittest.mock import patch

from pwchem.objects import VCFFile
from pwchem.protocols import ProtVEP
from pwchem.utils.utilsRNA import iterVepRows


VEP_HEADER = ('#Uploaded_variation\tLocation\tAllele\tGene\tFeature\tFeature_type\t'
              'Consequence\tProtein_position\tAmino_acids\tExtra\n')


class TestVEP(unittest.TestCase):
    """Check offline/cache configuration and VEP/filter_vep output contracts.

    Tool execution is mocked: a real Ensembl cache is not supplied with the
    genomics dataset. These tests do not establish annotation accuracy.
    Native table parsing and rejection of malformed records run directly.
    """

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(prefix='test_vep_')
        self.addCleanup(self.work.cleanup)
        self.protocol = ProtVEP(workingDir=self.work.name)
        os.makedirs(self.protocol._getExtraPath())
        cache = os.path.join(self.work.name, 'cache', 'homo_sapiens', '116_GRCh38')
        os.makedirs(cache)
        with open(os.path.join(cache, 'info.txt'), 'w') as handle:
            handle.write('test cache fixture; not a real annotation cache\n')
        self.protocol.cacheDirectory.set(os.path.join(self.work.name, 'cache'))
        path = os.path.join(self.work.name, 'input.vcf')
        with open(path, 'w') as handle:
            handle.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n'
                         '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')
        self.protocol.inputVCF.set(VCFFile(filename=path, scientificName='Homo sapiens', assembly='GRCh38'))
        self.commands = []

    def _runTool(self, program, args):
        self.commands.append((program, args))
        output = args[args.index('--output_file') + 1] if program == 'vep' else args[args.index('-o') + 1]
        with open(output, 'w') as handle:
            handle.write('## VEP test fixture\n' + VEP_HEADER)
            handle.write('v1\t1:10\tC\tG1\tT1\tTranscript\tmissense_variant\t2\tC/A\tENSP=P1;UNIPARC=UPI0000000001;SYMBOL=Gene1\n')
            if program == 'vep':
                handle.write('v2\t1:20\tA\tG1\tT1\tTranscript\tsynonymous_variant\t3\tD\tENSP=P1\n')

    def testOfflineAndMissenseCommands(self):
        self.assertEqual(self.protocol._validate(), [])
        with patch.object(self.protocol, '_runTool', side_effect=self._runTool):
            self.protocol.annotateStep()
        self.assertEqual([item[0] for item in self.commands], ['vep', 'filter_vep'])
        args = self.commands[0][1]
        for flag in ('--offline', '--cache', '--uniprot', '--protein', '--symbol', '--transcript_version', '--tsl'):
            self.assertIn(flag, args)
        self.assertIn('Consequence is missense_variant', self.commands[1][1])
        raw = list(iterVepRows(self.protocol._getExtraPath('annotation', 'all.tsv')))
        selected = list(iterVepRows(self.protocol._getExtraPath('annotation', 'missense.tsv')))
        self.assertEqual(len(raw), 2)
        self.assertEqual(len(selected), 1)
        self.assertEqual(selected[0]['ENSP'], 'P1')
        self.assertEqual(selected[0]['UNIPARC'], 'UPI0000000001')

    def testWithoutSelection(self):
        self.protocol.selectMissense.set(False)
        with patch.object(self.protocol, '_runTool', side_effect=self._runTool):
            self.protocol.annotateStep()
        self.assertEqual(len(self.commands), 1)
        self.assertEqual(len(list(iterVepRows(self.protocol._getExtraPath('annotation', 'all.tsv')))), 2)

    def testMissingCacheAndAssemblyMismatch(self):
        self.protocol.assembly.set('GRCh37')
        errors = self.protocol._validate()
        self.assertTrue(any('assemblies' in error for error in errors))
        self.assertTrue(any('info.txt' in error for error in errors))

    def testHeaderOnly(self):
        path = os.path.join(self.work.name, 'empty.tsv')
        with open(path, 'w') as handle:
            handle.write(VEP_HEADER)
        self.assertEqual(list(iterVepRows(path)), [])

    def testMalformedRows(self):
        path = os.path.join(self.work.name, 'invalid.tsv')
        for content in ('', 'v1\tbad\n', VEP_HEADER + 'v1\tbad\n'):
            with self.subTest(content=content), open(path, 'w') as handle:
                handle.write(content)
            with self.assertRaises(RuntimeError):
                list(iterVepRows(path))

    def testFailedRerunRetainsPreviousResults(self):
        with patch.object(self.protocol, '_runTool', side_effect=self._runTool):
            self.protocol.annotateStep()
        path = self.protocol._getExtraPath('annotation', 'all.tsv')
        with open(path) as handle:
            original = handle.read()
        with patch.object(self.protocol, '_runTool', side_effect=RuntimeError('Tool failed')):
            with self.assertRaises(RuntimeError):
                self.protocol.annotateStep()
        with open(path) as handle:
            self.assertEqual(handle.read(), original)

    def testRealAnnotationWhenResourcesAvailable(self):
        cache, vcf = os.environ.get('VEP_TEST_CACHE'), os.environ.get('VEP_TEST_VCF')
        searchPath = os.environ.get('RNASEQ_TEST_BIN', os.environ.get('PATH'))
        tools = {name: shutil.which(name, path=searchPath) for name in ('vep', 'filter_vep')}
        if not cache or not vcf or not all(tools.values()):
            self.skipTest('Set VEP_TEST_CACHE, VEP_TEST_VCF and provide vep/filter_vep for real annotation.')
        self.protocol.cacheDirectory.set(cache)
        self.protocol.inputVCF.set(VCFFile(filename=vcf))
        self.protocol.species.set(os.environ.get('VEP_TEST_SPECIES', 'homo_sapiens'))
        self.protocol.assembly.set(os.environ.get('VEP_TEST_ASSEMBLY', 'GRCh38'))
        self.protocol.cacheVersion.set(int(os.environ.get('VEP_TEST_RELEASE', '116')))
        self.protocol.numberOfThreads.set(1)
        def run(program, args):
            subprocess.run([tools[program]] + [str(arg) for arg in args], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with patch.object(self.protocol, '_runTool', side_effect=run):
            self.protocol.annotateStep()
        list(iterVepRows(self.protocol._getExtraPath('annotation', 'all.tsv')))
        for row in iterVepRows(self.protocol._getExtraPath('annotation', 'missense.tsv')):
            self.assertIn('missense_variant', row['Consequence'].replace('&', ',').split(','))
