# **************************************************************************
# *
# * Authors:    Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
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
import shutil
import tempfile

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import EnumParam, PointerParam

from pwchem.objects import VCFFile
from pwchem.utils.utilsRNA import GATKPreparationMixin


class ProtBCFToolsIsec(GATKPreparationMixin, EMProtocol):
    """Compare two VCFs by position and identical REF/ALT alleles.

    The default output reproduces the reference script's ``0003.vcf``:
    variants shared by both inputs, retaining the second input's sample,
    genotypes and annotations. Shared records from the first input, or
    records exclusive to either input, can also be selected. This operation
    does not combine genotypes or establish a tumor/normal relationship.

    Both inputs must belong to the same assembly. Available species and
    assembly metadata, and declared contig lengths, are checked before
    execution. Missing metadata cannot certify assembly compatibility.
    Each input is copied through bcftools sort into a local BGZF VCF and
    indexed with CSI. Original VCFs and indexes remain untouched. Plain VCF
    and BGZF VCF inputs are accepted, including files containing only headers.

    Matching uses ``--collapse none``: the complete REF/ALT alleles must
    match. Equivalent indels represented differently will not intersect;
    normalize both inputs against the same reference beforehand if needed.
    No PASS-only, population-frequency or known-variant filter is applied.
    Mutect2 filtering remains a separate, explicit upstream choice.

    outputVCF is the selected VCFFile; all four native isec partitions and
    its README remain in the protocol directory. Scientific metadata comes
    from the input whose records were selected. Requires bcftools in the
    plugin RNA-seq environment. Threads default to four.
    """

    _label = 'VCF comparison (bcftools isec)'
    PARTITIONS = ['Shared, records from second VCF', 'Shared, records from first VCF',
                  'Exclusive to first VCF', 'Exclusive to second VCF']
    FILES = ['0003.vcf', '0002.vcf', '0000.vcf', '0001.vcf']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputFirst', PointerParam, pointerClass='VCFFile', label='First VCF:')
        form.addParam('inputSecond', PointerParam, pointerClass='VCFFile', label='Second VCF:')
        form.addParam('selection', EnumParam, choices=self.PARTITIONS, default=0,
                      label='Records to retain:')
        form.addParallelSection(threads=4, mpi=1)

    def _getCommandFile(self):
        return self._getExtraPath('bcftools_commands.txt')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.compareStep)
        self._insertFunctionStep(self.createOutputStep)

    def _validate(self):
        errors = []
        inputs = [self.inputFirst.get(), self.inputSecond.get()]
        for index, obj in enumerate(inputs, 1):
            if obj is None or not obj.getFileName() or not os.path.isfile(obj.getFileName()):
                errors.append('VCF {} must be an existing file.'.format(index))
        if errors:
            return errors
        for getter, label in [('getScientificName', 'species'), ('getAssembly', 'assembly')]:
            a, b = [getattr(obj, getter)() for obj in inputs]
            if getter == 'getScientificName':
                a, b = [(value or '').strip().lower().replace('_', ' ') for value in (a, b)]
            if a and b and a != b:
                errors.append('VCF {} metadata differ.'.format(label))
        try:
            a, b = [self._readVCFContigs(obj.getFileName()) for obj in inputs]
            for name in a.keys() & b.keys():
                if a[name] and b[name] and a[name] != b[name]:
                    errors.append('Different lengths for contig {}.'.format(name))
            if a and b and not a.keys() & b.keys():
                errors.append('VCFs have no declared contigs in common.')
        except (OSError, ValueError, RuntimeError, EOFError) as error:
            errors.append(str(error))
        if self.selection.get() not in range(len(self.FILES)):
            errors.append('Invalid VCF partition.')
        if self.numberOfThreads.get() < 1:
            errors.append('At least one thread is required.')
        return errors

    def compareStep(self):
        os.makedirs(self._getExtraPath(), exist_ok=True)
        errors = self._validate()
        if errors:
            raise RuntimeError('; '.join(errors))
        with open(self._getCommandFile(), 'w'):
            pass
        with tempfile.TemporaryDirectory(prefix='isec_', dir=self._getExtraPath()) as work:
            paths = []
            for index, obj in enumerate([self.inputFirst.get(), self.inputSecond.get()]):
                target = os.path.join(work, 'input{}.vcf.gz'.format(index + 1))
                self._runTool('bcftools', ['sort', '-Oz', '-o', target, obj.getFileName()])
                self._runTool('bcftools', ['index', '--csi', '--threads',
                                          self.numberOfThreads.get(), target])
                paths.append(target)
            result = os.path.join(work, 'partitions')
            self._runTool('bcftools', ['isec', '--collapse', 'none', '-Ov', '-p', result] + paths)
            for filename in self.FILES:
                self._readVCFContigs(os.path.join(result, filename))
            destination = self._getExtraPath('partitions')
            if os.path.isdir(destination):
                shutil.rmtree(destination)
            shutil.move(result, destination)

    def createOutputStep(self):
        selection = self.selection.get()
        source = self.inputSecond.get() if selection in (0, 3) else self.inputFirst.get()
        output = source.clone()
        output.setFileName(self._getExtraPath('partitions', self.FILES[selection]))
        output.setIndexFile('')
        output.setIsCompressed(False)
        output.setSource('bcftools isec')
        self._defineOutputs(outputVCF=output)
        for obj in (self.inputFirst.get(), self.inputSecond.get()):
            self._defineSourceRelation(obj, output)

    def _summary(self):
        return [self.PARTITIONS[self.selection.get()]] if hasattr(self, 'outputVCF') else ['Not finished.']
