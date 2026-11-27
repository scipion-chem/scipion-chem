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
import gzip
import math
import os
import re
import shutil
import tempfile

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import FileParam, FloatParam, IntParam, PointerParam, StringParam

from pwchem.objects import MutationReport
from pwchem.utils.utilsRNA import iterVepRows


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
            self.buildMutationReport(self.inputAnnotations.get().getFileName(), self.proteinFasta.get(), work,
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

    @classmethod
    def readSymbolExpression(cls, path, column=''):
        """Read first FPKM per exact gene symbol, matching the reference script.

        The source is isoforms.fpkm_tracking; retaining the first symbol is not
        a gene-level sum. Cuffdiff needs an explicit condition FPKM column.
        """
        import csv
        values = {}
        with open(path) as handle:
            reader = csv.DictReader(handle, delimiter='\t')
            column = cls.getExpressionColumn(reader.fieldnames or [], column)
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise RuntimeError('Malformed expression table row.')
                symbol = row['gene_short_name']
                if symbol in ('', '-') or symbol in values:
                    continue
                values[symbol] = cls.parseExpressionValue(row[column], symbol)
        return values


    @staticmethod
    def getExpressionColumn(fields, column):
        column = column or ('FPKM' if 'FPKM' in fields else '')
        if len(fields) != len(set(fields)) or 'gene_short_name' not in fields:
            raise RuntimeError('Expression table needs unique columns and gene_short_name.')
        if not column or column not in fields or not (column == 'FPKM' or column.endswith('_FPKM')):
            raise RuntimeError('Select an existing FPKM column (Cuffdiff requires a condition).')
        return column


    @staticmethod
    def parseExpressionValue(value, symbol):
        import math
        try:
            value = float(value)
            if not math.isfinite(value) or value < 0:
                raise ValueError
        except ValueError as error:
            raise RuntimeError('FPKM must be finite and nonnegative for {}.'.format(symbol)) from error
        return value


    @staticmethod
    def getUniParcId(value):
        """Return one unambiguous UniParc accession, never choose among isoforms."""
        accessions = set((value or '').split(','))
        if len(accessions) == 1:
            accession = next(iter(accessions))
            if re.fullmatch(r'UPI[0-9A-F]{10,}', accession):
                return accession
        return None


    @staticmethod
    def readUniParcFasta(path, requested):
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


    @classmethod
    def fetchUniParcSequence(cls, accession, cacheDirectory, session, attempts=3):
        """Fetch one exact UPI sequence, with timeout, bounded retries and cache.

        Queries the same UniProt search endpoint as the supplied script. Responses
        containing several FASTA entries or another accession are errors. A 404
        returns None; exhausted network/server failures stop the run rather than
        silently dropping mutations. Successfully validated responses are cached.
        """
        if cls.getUniParcId(accession) != accession:
            raise ValueError('Invalid UniParc accession.')
        os.makedirs(cacheDirectory, exist_ok=True)
        path = os.path.join(cacheDirectory, accession + '.fasta')
        if os.path.isfile(path):
            cached = cls.readUniParcFasta(path, {accession})
            if accession not in cached:
                raise RuntimeError('Cached UniParc file does not contain {}.'.format(accession))
            return cached[accession]
        if attempts < 1:
            raise ValueError('At least one download attempt is required.')
        sequence = cls.downloadUniParcSequence(accession, session, attempts)
        if sequence is not None:
            cls.cacheUniParcSequence(path, accession, sequence)
        return sequence


    @classmethod
    def downloadUniParcSequence(cls, accession, session, attempts):
        import requests
        import time
        for attempt in range(attempts):
            try:
                status, sequence = cls.getUniParcResponse(accession, session)
            except requests.RequestException as error:
                if attempt + 1 == attempts:
                    raise RuntimeError('UniParc download failed for {}.'.format(accession)) from error
            else:
                if status != 429 and status < 500:
                    return sequence
                if attempt + 1 == attempts:
                    raise RuntimeError('UniParc returned HTTP {} for {}.'.format(status, accession))
            time.sleep(attempt + 1)


    @classmethod
    def getUniParcResponse(cls, accession, session):
        import requests
        response = session.get('https://rest.uniprot.org/uniparc/search',
                               params={'query': accession, 'format': 'fasta'}, timeout=(10, 60))
        try:
            status = response.status_code
            if status == 404 or status == 429 or status >= 500:
                return status, None
            response.raise_for_status()
            return status, cls.parseUniParcResponse(response.text, accession)
        except requests.HTTPError as error:
            raise RuntimeError('UniParc rejected {} (HTTP {}).'.format(accession, response.status_code)) from error
        finally:
            response.close()


    @staticmethod
    def parseUniParcResponse(text, accession):
        import io
        from Bio import SeqIO
        records = list(SeqIO.parse(io.StringIO(text), 'fasta'))
        if not records and not text.strip():
            return None
        if len(records) != 1 or accession not in re.findall(r'\bUPI[0-9A-F]{10,}\b', records[0].id):
            raise RuntimeError('UniParc response is not a single matching protein for {}.'.format(accession))
        sequence = str(records[0].seq).upper()
        if not sequence:
            raise RuntimeError('Empty UniParc sequence for {}.'.format(accession))
        return sequence


    @staticmethod
    def cacheUniParcSequence(path, accession, sequence):
        with tempfile.NamedTemporaryFile(mode='w', dir=os.path.dirname(path), delete=False) as handle:
            temporary = handle.name
            handle.write('>{}\n{}\n'.format(accession, sequence))
        try:
            os.replace(temporary, path)
        finally:
            if os.path.exists(temporary):
                os.remove(temporary)


    @classmethod
    def buildMutationReport(cls, annotations, proteinFasta, directory, expression,
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
        import math
        if not isinstance(peptideLength, int) or isinstance(peptideLength, bool) or peptideLength < 1 or peptideLength % 2 == 0:
            raise ValueError('Peptide length must be a positive odd integer.')
        if not math.isfinite(minimumFPKM) or minimumFPKM < 0:
            raise ValueError('Minimum FPKM must be finite and nonnegative.')
        if not expression:
            raise ValueError('Cufflinks/Cuffdiff expression is required.')
        values = cls.readSymbolExpression(expression, expressionColumn)
        os.makedirs(directory, exist_ok=True)
        requested = cls.getRequestedUniParcIds(annotations)
        proteins = cls.getMutationProteins(proteinFasta, requested, cacheDirectory or os.path.join(directory, 'uniparc_cache'))
        with open(os.path.join(directory, 'reference_sequences.fasta'), 'w') as handle:
            for accession, sequence in sorted(proteins.items()):
                handle.write('>{}\n{}\n'.format(accession, sequence))
        return cls.writeMutationReport(annotations, proteins, values, directory, peptideLength, minimumFPKM)


    @classmethod
    def getRequestedUniParcIds(cls, annotations):
        requested = {cls.getUniParcId(row.get('UNIPARC')) for row in iterVepRows(annotations)
                     if row.get('SYMBOL') not in (None, '', '-')
                     and 'missense_variant' in row['Consequence'].replace('&', ',').split(',')}
        requested.discard(None)
        return requested


    @classmethod
    def getMutationProteins(cls, proteinFasta, requested, cacheDirectory):
        import requests
        if proteinFasta:
            return cls.readUniParcFasta(proteinFasta, requested)
        with requests.Session() as session:
            return {accession: sequence for accession in sorted(requested)
                    for sequence in [cls.fetchUniParcSequence(accession, cacheDirectory, session)]
                    if sequence is not None}


    @classmethod
    def getMutationChangeError(cls, row):
        aminoAcids = set('ACDEFGHIKLMNPQRSTVWY')
        change = row['Amino_acids'].split('/')
        position = row['Protein_position']
        if 'missense_variant' not in row['Consequence'].replace('&', ',').split(','):
            return 'not_missense'
        if row['Feature_type'] != 'Transcript':
            return 'not_transcript'
        if not cls.getUniParcId(row.get('UNIPARC')):
            return 'missing_or_ambiguous_uniparc'
        if not row.get('SYMBOL') or row['SYMBOL'] == '-':
            return 'missing_symbol'
        if len(change) != 2 or any(aa not in aminoAcids for aa in change) or change[0] == change[1]:
            return 'unsupported_amino_acid_change'
        if not position.isdigit() or int(position) < 1:
            return 'unsupported_protein_position'
        return ''


    @classmethod
    def getMutationPeptide(cls, row, proteins, values, flank, minimumFPKM):
        reason = cls.getMutationChangeError(row)
        if reason:
            return None, reason
        protein = cls.getUniParcId(row.get('UNIPARC'))
        sequence = proteins.get(protein, '')
        if not sequence:
            return None, 'missing_uniparc_protein'
        position = int(row['Protein_position']) - 1
        change = row['Amino_acids'].split('/')
        if position >= len(sequence) or sequence[position] != change[0]:
            return None, 'reference_residue_mismatch'
        start, end = position - flank, position + flank + 1
        if start < 0 or end > len(sequence):
            return None, 'incomplete_peptide_window'
        reference = sequence[start:end]
        if not set(reference).issubset(set('ACDEFGHIKLMNPQRSTVWY')):
            return None, 'unsupported_residue_in_window'
        fpkm = values.get(row['SYMBOL'], 0.0)
        if fpkm < minimumFPKM:
            return None, 'below_expression_threshold'
        mutant = reference[:flank] + change[1] + reference[flank + 1:]
        return (protein, position, change, fpkm, start, end, reference, mutant), ''


    @classmethod
    def writeMutationReport(cls, annotations, proteins, values, directory, peptideLength, minimumFPKM):
        import csv
        fields = ['effect_id', 'variant', 'location', 'allele', 'gene_id', 'symbol', 'transcript_id',
                  'protein_id', 'consequence', 'protein_position', 'reference_aa', 'mutant_aa',
                  'FPKM', 'window_start', 'window_end', 'mutation_offset', 'reference_peptide', 'mutant_peptide']
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
                peptide, reason = cls.getMutationPeptide(row, proteins, values, flank, minimumFPKM)
                if reason:
                    rejectWriter.writerow([row['Uploaded_variation'], row['Feature'], row.get('UNIPARC', ''), reason])
                    continue
                protein, position, change, fpkm, start, end, reference, mutant = peptide
                accepted += 1
                effect = 'effect_{:06d}'.format(accepted)
                writer.writerow([effect, row['Uploaded_variation'], row['Location'], row['Allele'],
                    row['Gene'], row['SYMBOL'], row['Feature'], protein, row['Consequence'], position + 1,
                    change[0], change[1], fpkm, start + 1, end, flank + 1, reference, mutant])
                fasta.write('>{}_WT\n{}\n>{}_MUT\n{}\n'.format(effect, reference, effect, mutant))
                if mutant not in seen:
                    candidateWriter.writerow([mutant, row['SYMBOL'], fpkm])
                    seen.add(mutant)
        return accepted
