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
import math
import os
import shutil
import tempfile

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import FileParam, FloatParam, IntParam, PointerParam, StringParam

from pwchem.objects import MutationReport
from pwchem.utils.utilsRNA import buildMutationReport


class ProtVEPToAminoacids(EMProtocol):
    """Convert VEP missense effects to mutant peptides and Cufflinks FPKM.

    Reproduces the supplied VEPtoAminoacids.py data flow: retrieve the UniParc
    protein named by VEP Extra/UNIPARC, substitute the reported amino acid,
    match Extra/SYMBOL to Cufflinks gene_short_name, and write a CSV with
    seq,symbol,FPKM. The first Cufflinks row for each exact symbol is retained
    and unmatched symbols receive zero, as in the reference. This is an
    isoform FPKM associated with a symbol, not a sum of gene isoforms. Select
    an explicit condition FPKM column for Cuffdiff; plain FPKM is automatic
    for a single-sample Cufflinks table. No additional quantifier is needed.

    The requested default is exactly 25 amino acids: 12 reference residues,
    the mutant residue at position 13, and 12 reference residues. peptideLength
    may be changed to another positive odd length. The supplied Python script
    instead hardcodes a 12-residue asymmetric slice and does not define the
    -w or -n arguments used by the earlier shell script. Those incompatible
    slicing and argument conventions are not copied into this protocol.

    Proteins can be read from a local UniParc FASTA for offline execution.
    If proteinFasta is empty, each unique UPI accession is retrieved through
    the UniProt REST search endpoint, with timeout, bounded retries and a
    protocol-local cache. A snapshot of the retrieved proteins is published
    with the outputs so a later run can use it without the network. Ensembl
    ENSP identifiers alone cannot substitute for UniParc accessions. Inputs
    must use compatible species, assembly and expression annotations.

    The reference residue must match the UniParc sequence. Missing or
    ambiguous UPI accessions, absent symbols, multi-residue/uncertain changes,
    unsupported amino acids and mutations too near protein termini for a
    complete window are recorded in rejected.tsv. Windows are neither padded
    nor shortened. Download failures stop execution rather than reuse a
    previous variant's sequence, which could happen in the original script.
    Separate mutations are applied independently; haplotypes, frameshifts and
    indels are not reconstructed. An optional minimumFPKM explicitly filters
    candidates; the default zero retains absent-expression symbols.

    outputMutations is a MutationReport. Its primary candidates.csv contains
    one row per unique mutant peptide, retaining the first effect's symbol
    and FPKM, matching drop_duplicates(subset='seq'). Related mutations.tsv
    keeps every accepted transcript effect, coordinates, UniParc accession,
    expression and WT/MUT windows. peptides.fasta links paired sequences to
    effect identifiers; rejected.tsv preserves exclusion reasons. The protein
    snapshot reference_sequences.fasta and conversion settings are retained.
    Header-only annotations yield valid empty tables and FASTA files. Native
    VEP headers remain unchanged; no sed preprocessing or pandas is required.
    """

    _label = 'VEP mutations / amino-acid windows'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAnnotations', PointerParam, pointerClass='VariantAnnotation', label='VEP annotations:')
        form.addParam('inputExpression', PointerParam, pointerClass='RNASeqExpression', label='Cufflinks/Cuffdiff expression:')
        form.addParam('proteinFasta', FileParam, default='', allowsNull=True, label='UniParc FASTA (optional):',
                      help='Leave empty to download unique UPI proteins; a sequence snapshot is saved.')
        form.addParam('expressionColumn', StringParam, default='', label='FPKM column:',
                      help='Empty selects FPKM for Cufflinks. Cuffdiff requires an explicit condition column.')
        form.addSection(label='Peptide selection')
        form.addParam('peptideLength', IntParam, default=25, label='Peptide length (aa):',
                      help='Odd length; 25 gives 12 residues on each side of the mutation.')
        form.addParam('minimumFPKM', FloatParam, default=0, label='Minimum FPKM:')

    def _validate(self):
        errors = []
        obj, expression = self.inputAnnotations.get(), self.inputExpression.get()
        if obj is None or not obj.getFileName() or not os.path.isfile(obj.getFileName()):
            errors.append('An existing VEP table is required.')
        if self.proteinFasta.get() and not os.path.isfile(self.proteinFasta.get()):
            errors.append('UniParc FASTA not found.')
        length = self.peptideLength.get()
        if length is None or length < 1 or length % 2 == 0:
            errors.append('Peptide length must be a positive odd integer.')
        if self.minimumFPKM.get() is None or not math.isfinite(self.minimumFPKM.get()) or self.minimumFPKM.get() < 0:
            errors.append('Minimum FPKM must be finite and nonnegative.')
        if expression is None or not expression.getFileName() or not os.path.isfile(expression.getFileName()):
            errors.append('An existing Cufflinks/Cuffdiff expression table is required.')
        elif obj is not None:
            self._validateExpressionReference(obj, expression, errors)
        return errors

    @staticmethod
    def _validateExpressionReference(obj, expression, errors):
        for getter in ('getScientificName', 'getAssembly'):
            a, b = getattr(obj, getter)(), getattr(expression, getter)()
            if getter == 'getScientificName':
                a, b = [(value or '').strip().lower().replace('_', ' ') for value in (a, b)]
            if a and b and a != b:
                errors.append('Annotation and expression reference metadata differ.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertStep)
        self._insertFunctionStep(self.createOutputStep)

    def convertStep(self):
        os.makedirs(self._getExtraPath(), exist_ok=True)
        errors = self._validate()
        if errors:
            raise RuntimeError('; '.join(errors))
        with tempfile.TemporaryDirectory(prefix='peptides_', dir=self._getExtraPath()) as work:
            buildMutationReport(self.inputAnnotations.get().getFileName(), self.proteinFasta.get(), work,
                expression=self.inputExpression.get().getFileName(), expressionColumn=self.expressionColumn.get(),
                minimumFPKM=self.minimumFPKM.get(), peptideLength=self.peptideLength.get(),
                cacheDirectory=self._getExtraPath('uniparc_cache'))
            destination = self._getExtraPath('mutations')
            if os.path.isdir(destination):
                shutil.rmtree(destination)
            shutil.move(work, destination)

    def createOutputStep(self):
        obj, expression = self.inputAnnotations.get(), self.inputExpression.get()
        output = MutationReport(filename=self._getExtraPath('mutations', 'candidates.csv'),
            mutationsFile=self._getExtraPath('mutations', 'mutations.tsv'),
            peptidesFile=self._getExtraPath('mutations', 'peptides.fasta'),
            rejectedFile=self._getExtraPath('mutations', 'rejected.tsv'),
            scientificName=obj.getScientificName(), assembly=obj.getAssembly(),
            peptideLength=self.peptideLength.get(), flankSize=self.peptideLength.get() // 2,
            proteinFasta=self._getExtraPath('mutations', 'reference_sequences.fasta'),
            expressionColumn=self.expressionColumn.get() or 'FPKM', minimumFPKM=self.minimumFPKM.get())
        self._defineOutputs(outputMutations=output)
        self._defineSourceRelation(obj, output)
        self._defineSourceRelation(expression, output)

    def _summary(self):
        return ['UniParc mutant peptides: {} aa; Cufflinks FPKM matched by symbol.'.format(self.peptideLength.get())]
