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
from pwchem.protocols import ProtBCFToolsIsec


class TestBCFToolsIsec(unittest.TestCase):
    """Check VCF compatibility, partition provenance and real isec if available.

    Synthetic VCFs isolate allele matching and genotype provenance without
    downloading a genome. Set RNASEQ_TEST_BIN to the RNA-seq bin directory
    for the external-tool test; other cases exercise validation and output
    wiring without pretending to execute bcftools.
    """

    def setUp(self):
        self.work = tempfile.TemporaryDirectory(prefix='test_isec_')
        self.addCleanup(self.work.cleanup)
        self.protocol = ProtBCFToolsIsec(workingDir=self.work.name)
        os.makedirs(self.protocol._getExtraPath())
        self.inputs = []
        for index in (1, 2):
            path = os.path.join(self.work.name, '{}.vcf'.format(index))
            with open(path, 'w') as handle:
                handle.write('##fileformat=VCFv4.2\n##contig=<ID=1,length=100>\n'
                    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
                    '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS{}\n'.format(index))
                handle.write('1\t10\t.\tA\tC\t.\tPASS\t.\tGT\t{}\n'.format('0/1' if index == 1 else '1/1'))
                handle.write('1\t20\t.\tG\t{}\t.\tPASS\t.\tGT\t0/1\n'.format('A' if index == 1 else 'T'))
            obj = VCFFile(filename=path, scientificName='Mus musculus', assembly='GRCm39')
            self.inputs.append(obj)
        self.protocol.inputFirst.set(self.inputs[0])
        self.protocol.inputSecond.set(self.inputs[1])

    def testCompatibleInputs(self):
        self.assertEqual(self.protocol._validate(), [])

    def testRejectAssemblyMismatch(self):
        self.inputs[1].setAssembly('GRCh38')
        self.assertTrue(any('assembly' in message for message in self.protocol._validate()))

    def testRejectContigLengthMismatch(self):
        path = self.inputs[1].getFileName()
        with open(path) as handle:
            text = handle.read()
        with open(path, 'w') as handle:
            handle.write(text.replace('length=100', 'length=200'))
        self.assertTrue(any('length' in message for message in self.protocol._validate()))

    def testPartitionMetadata(self):
        for selection, sourceIndex in enumerate((1, 0, 0, 1)):
            with self.subTest(selection=selection):
                self.protocol.selection.set(selection)
                self.inputs[0].setTumorSample('sample_one')
                self.inputs[1].setTumorSample('sample_two')
                outputs = {}
                with patch.object(self.protocol, '_defineOutputs', side_effect=lambda **kw: outputs.update(kw)), \
                        patch.object(self.protocol, '_defineSourceRelation'):
                    self.protocol.createOutputStep()
                output = outputs['outputVCF']
                self.assertTrue(output.getFileName().endswith(self.protocol.FILES[selection]))
                self.assertEqual(output.getTumorSample(), self.inputs[sourceIndex].getTumorSample())
                self.assertFalse(output.isCompressed())
                self.assertFalse(output.getIndexFile())

    def testInvalidSelection(self):
        self.protocol.selection.set(99)
        self.assertTrue(self.protocol._validate())

    def testRealAlleleIntersection(self):
        binary = shutil.which('bcftools', path=os.environ.get('RNASEQ_TEST_BIN', os.environ.get('PATH')))
        if not binary:
            self.skipTest('bcftools not available; set RNASEQ_TEST_BIN.')
        def run(program, args):
            subprocess.run([binary] + [str(arg) for arg in args], check=True,
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        originals = []
        for obj in self.inputs:
            with open(obj.getFileName(), 'rb') as handle:
                originals.append(handle.read())
        with patch.object(self.protocol, 'runTool', side_effect=run):
            self.protocol.compareStep()
        for obj, original in zip(self.inputs, originals):
            with open(obj.getFileName(), 'rb') as handle:
                self.assertEqual(handle.read(), original)
        with open(self.protocol._getExtraPath('partitions', '0003.vcf')) as handle:
            text = handle.read()
        self.assertIn('S2', text)
        records = [line.split('\t') for line in text.splitlines() if not line.startswith('#')]
        self.assertEqual(len(records), 1)
        self.assertEqual(records[0][1], '10')
        self.assertEqual(records[0][-1], '1/1')
        for filename in ('0000.vcf', '0001.vcf'):
            with open(self.protocol._getExtraPath('partitions', filename)) as handle:
                records = [line for line in handle if not line.startswith('#')]
            self.assertEqual(len(records), 1)
            self.assertEqual(records[0].split('\t')[1], '20')
