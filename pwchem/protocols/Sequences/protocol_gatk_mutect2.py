
# **************************************************************************
# *
# * Authors:    Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia, CSIC
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
# *
# **************************************************************************

import gzip
import json
import os
import re

from pyworkflow.protocol.params import (
    BooleanParam, EnumParam, IntParam, PointerParam, StringParam
)
from pwem.protocols import EMProtocol

from pwchem.utils.utilsRNA import GATKPreparationMixin
from pwchem.objects import VCFFile


class ProtGATKMutect2(GATKPreparationMixin, EMProtocol):
    """
    Call candidate somatic SNVs and indels using GATK Mutect2.

    This protocol accepts a tumor ``AlignmentFile`` and optionally a
    matched normal alignment. The reference genome FASTA is obtained from
    the tumor alignment metadata. Results are stored as ``VCFFile`` objects.

    Calling can be supported by a germline resource and a panel of normals.
    ``FilterMutectCalls`` is disabled by default, matching the reference script. The reported variant types
    can be restricted to SNVs or insertions/deletions after calling and
    optional filtering.


    ---------------------------------------------------------------------
    Calling modes
    ---------------------------------------------------------------------

    Tumor-only
        Calls candidate variants using the tumor BAM without a matched
        normal. This is the default mode. Distinguishing somatic variants
        from inherited variants requires particular care in this mode.

    Tumor-normal
        Uses tumor and matched normal BAMs from the same individual.
        Enable ``useNormal`` and provide ``inputNormal``. Each BAM must contain
        exactly one sample, identified by its ``SM`` read-group tags, and tumor
        and normal sample names must differ. Multiple read groups for the same
        sample are supported; their IDs must be unique within each BAM.

        ``normalSample`` is optional: when empty, the single normal SM name is
        detected automatically. An explicit value must match the header exactly.

        Both alignments must refer to the same reference FASTA file, allowing
        symbolic links and hard links. The two BAMs must be different files;
        path aliases cannot be used to supply the tumor as its own normal.


    ---------------------------------------------------------------------
    RNA-seq calling options and reference-script reproduction
    ---------------------------------------------------------------------

    disableMateOnSameContigFilter : bool, default True
        Adds the following arguments, matching the reference RNA-seq script::

            --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter

        This disables that individual mate-placement filter, rather than all
        default GATK read filters. Disable the parameter to retain the tool's
        default read-filter configuration.

    dontUseSoftClippedBases : bool, default True
        Adds ``--dont-use-soft-clipped-bases`` as in the reference script.
        Disable the parameter to allow the default handling of soft-clipped bases.

    To reproduce the calling stage of the script, retain both options above,
    select ``SNVs and indels`` and leave ``filterCalls`` disabled (the default). Germline resource
    and panel of normals remain disabled unless additional resources are wanted.
    For tumor-normal calling, provide both prepared BAMs and their actual SM
    names; the normal is not required to be literally named ``normal``.

    ``FilterMutectCalls`` can be enabled explicitly as an optional extension
    to the script. It annotates filters without producing a PASS-only VCF.
    This protocol covers somatic calling; expression analysis, VCF intersections,
    VEP annotation and peptide generation belong to separate workflow stages.
    GATK versions differ from the original script, so identical settings do not
    guarantee byte-identical outputs.


    ---------------------------------------------------------------------
    Variant selection and filtering
    ---------------------------------------------------------------------

    variantType
        Selects which variant records are published in both output VCFs:

        - ``SNVs and indels``: retains all Mutect2 calls (default), including
          mixed records and multi-nucleotide substitutions if present.
        - ``SNVs only``: retains records classified as ``SNP`` by GATK.
        - ``Indels only``: retains records classified as ``INDEL`` by GATK.

        Exclusive modes use GATK ``SelectVariants``. Selection operates on
        complete VCF records rather than individual alternative alleles.
        Mixed records and multi-nucleotide substitutions are excluded from
        both exclusive modes.

        Mutect2 always performs variant calling before this selection.

    filterCalls
        Runs ``FilterMutectCalls`` on the complete Mutect2 call set before
        selecting variant types. Disabled by default.

        Variants that fail filtering remain in the VCF with ``FILTER``
        annotations. The protocol does not produce a PASS-only VCF.


    ---------------------------------------------------------------------
    Input
    ---------------------------------------------------------------------

    inputTumor : AlignmentFile
        Tumor BAM prepared for variant calling, with valid sample names in
        its read groups and an existing reference FASTA in its metadata.
        The BAM and reference must exist and be non-empty. Required indexes
        are copied or generated locally; input directories need not be writable.
        The BAM must already be coordinate-sorted and prepared for somatic calling;
        this protocol does not mark duplicates, split N CIGARs or perform BQSR.

    inputNormal : AlignmentFile, optional
        Matched normal BAM, required when ``useNormal`` is enabled.

    germlineVCFs : SetOfVCFFiles, optional
        Germline population resource used when ``useGermlineResource`` is
        enabled. Select one VCF using the zero-based ``germlineIndex``.
        The resource must be suitable for Mutect2, including the population
        allele-frequency annotations it expects, and match the assembly.
        Its header must declare INFO/AF with Number=A and Type=Float. This
        checks the declaration, not every annotation value or population quality.
        A dbSNP/known-sites VCF used for BQSR is not automatically suitable.

    ponVCFs : SetOfVCFFiles, optional
        Panel of normals used when ``usePanelOfNormals`` is enabled.
        Select one VCF using the zero-based ``ponIndex``. This must be a
        panel prepared for somatic calling, rather than an arbitrary VCF
        from a single normal sample, and must match the reference assembly.


    ---------------------------------------------------------------------
    Execution order
    ---------------------------------------------------------------------

    The protocol executes the following steps::

        Prepare reference and optional resource indexes
            -> Mutect2
            -> FilterMutectCalls (optional)
            -> SelectVariants (SNVs-only or indels-only modes)
            -> Create output VCFFile objects

    Reference, tumor, normal and resource inputs are linked into the protocol
    working directory. Preparation shares the same routines as GATK preprocessing,
    avoiding duplicated index and sequence-dictionary handling.

    Missing, stale, empty or structurally invalid reference ``.fai`` and ``.dict``
    files are regenerated with samtools ``faidx`` and GATK ``CreateSequenceDictionary``.
    Their sequence names, lengths and order must agree. BAM indexes accept both
    file.bam.bai and file.bai naming conventions; unusable BAI files are regenerated
    locally with samtools. Input BAMs undergo samtools quickcheck and header checks.

    Plain resource VCFs require a Tribble ``.idx``; BGZF resource VCFs require a
    tabix ``.tbi``. A lone ``.csi`` triggers GATK ``IndexFeatureFile``. Ordinary
    gzip is rejected because it is not suitable for tabix indexing. Source files
    and source indexes are not modified.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputUnfilteredVCF : VCFFile
        Mutect2 calls before ``FilterMutectCalls``, restricted to the
        selected variant types when an exclusive mode is used.

    outputVCF : VCFFile
        Calls annotated by ``FilterMutectCalls`` when filtering is enabled,
        restricted to the selected variant types. When filtering is
        disabled, this is a clone of ``outputUnfilteredVCF`` and refers to
        the same file.

    Both outputs store the prepared reference FASTA, tabix index, compression
    status, actual tumor/normal sample names, GATK/Mutect2 provenance, selected
    variant type and whether FilterMutectCalls was applied. The filtering flag
    means filters were annotated, not that only passing variants were retained.

    Scientific name and assembly are copied from structured tumor reference
    metadata when available. RNA-seq alignment now records these fields and
    Picard/GATK preserve them through cloning. Older alignments without structured
    reference metadata leave these fields unset rather than inferring them from
    filenames or resource databases.

    Output VCF headers and sample names are checked and a compatible index is
    ensured before publication. A VCF with a valid header and no variant records
    is a valid result. Source relations link both outputs to input alignments
    and enabled VCF resources.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - A tumor alignment is provided and its BAM exists and is not empty.
    - The tumor alignment contains an existing reference FASTA.
    - When matched-normal mode is enabled, the normal BAM exists and is
      not empty, its path differs from the tumor BAM path, its reference
      resolves to the same FASTA file. Sample-name checks occur during preparation.
    - Each enabled resource contains a non-empty ``SetOfVCFFiles``, its
      selection index is in range and the selected VCF exists and is non-empty.
    - A positive thread count and a valid variant-type choice are selected.

    During preparation, BAM read groups and sample identities are inspected.
    Sequence names, lengths and order must agree with the reference dictionary.
    Resource VCF contig declarations are compared to the reference when present;
    mismatched shared lengths or no common sequence names stop execution. Extra
    resource contigs generate a warning and remain subject to GATK validation.
    Missing VCF contig declarations generate a compatibility warning.

    These checks do not certify biological pairing, panel-of-normals construction,
    reference sequence identity or the quality of germline allele frequencies.
    Index signature, timestamp and header checks do not replace complete file
    validation. Mutect2/FilterMutectCalls retain their own validation checks.


    ---------------------------------------------------------------------
    Intermediate files
    ---------------------------------------------------------------------

    mutect2_unfiltered.vcf.gz
        Complete Mutect2 call set before filtering or type selection.

    mutect2_unfiltered.vcf.gz.stats
        Mutect2 statistics used by ``FilterMutectCalls``. Checked for existence
        and non-empty content and passed explicitly through ``--stats`` when
        filtering is enabled. Type selection never replaces this complete call set.

    mutect2_filtered.vcf.gz
        Complete call set annotated by ``FilterMutectCalls``, when enabled.

    mutect2_unfiltered_snvs.vcf.gz / mutect2_unfiltered_indels.vcf.gz
        Type-selected unfiltered calls, generated in the exclusive modes.

    mutect2_filtered_snvs.vcf.gz / mutect2_filtered_indels.vcf.gz
        Type-selected filtered calls, generated when filtering and an
        exclusive mode are enabled.

    gatk_mutect2_commands.txt
        Internal record of GATK and samtools commands executed by the protocol.

    Previously generated VCFs, indexes and statistics at the output paths
    are removed before the corresponding tool runs, preventing stale files
    from hiding a missing result. Input resources are not removed.

    mutect2_samples.json
        Actual SM names detected from tumor and optional normal BAM headers,
        persisted atomically for subsequent steps.

    reference/reference.fa, reference/reference.fa.fai, reference/reference.dict
        Local reference link, FASTA index and sequence dictionary.

    inputs/tumor.bam, inputs/normal.bam
        Local alignment links with their prepared indexes and extracted headers.

    inputs/germline.vcf[.gz], inputs/pon.vcf[.gz]
        Local links to enabled resources and their prepared GATK-compatible indexes.


    ---------------------------------------------------------------------
    Requirements
    ---------------------------------------------------------------------

    The following programs must be available through the RNA-seq
    environment configured by the Scipion-Chem plugin:

    - GATK
    - samtools

    The protocol working directory must be writable. Reference, alignment and
    resource directories only require read access. Reference FASTA must be
    uncompressed; resource VCFs can be plain or BGZF-compressed.


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Outputs contain candidate somatic variants, not confirmed mutations.
    - RNA-seq analysis requires appropriate alignment preparation and
      interpretation of expression, RNA editing and alignment artifacts.
    - Additional evaluation is especially relevant for tumor-only analyses.
    - This protocol does not estimate contamination or learn an orientation
      bias model for use by ``FilterMutectCalls``.
    - Structural variant, copy-number and gene-fusion calling and functional
      consequence annotation are outside the scope of this protocol.
    - Complete call sets remain available in the working directory when
      type-selected VCFs are published as outputs.
    """

    _label = 'GATK Mutect2'

    VARIANT_TYPES = ['SNVs and indels', 'SNVs only', 'Indels only']

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam(
            'inputTumor', PointerParam, pointerClass='AlignmentFile',
            label='Tumor alignment:',
            help='Tumor BAM (AlignmentFile) previously prepared for variant calling.'
        )
        form.addParam(
            'useNormal', BooleanParam, default=False,
            label='Use matched normal:',
            help='Enable tumor-normal variant calling.'
        )
        form.addParam(
            'inputNormal', PointerParam, pointerClass='AlignmentFile',
            condition='useNormal', label='Normal alignment:',
            help='Matched normal BAM with the same reference assembly.'
        )
        form.addParam(
            'normalSample', StringParam, default='', allowsNull=True,
            condition='useNormal', label='Normal sample name:',
            help='Optional SM name from the normal BAM header. Leave empty to detect '
                 'its single sample automatically. An explicit value must match exactly.'
        )

        form.addSection(label='Optional resources')
        form.addParam(
            'useGermlineResource', BooleanParam, default=False,
            label='Use germline resource:'
        )
        form.addParam(
            'germlineVCFs', PointerParam, pointerClass='SetOfVCFFiles',
            condition='useGermlineResource', label='Germline VCFs:'
        )
        form.addParam(
            'germlineIndex', IntParam, default=0,
            condition='useGermlineResource', label='Germline VCF index:',
            help='Zero-based index of the resource in the VCF set.'
        )
        form.addParam(
            'usePanelOfNormals', BooleanParam, default=False,
            label='Use panel of normals:'
        )
        form.addParam(
            'ponVCFs', PointerParam, pointerClass='SetOfVCFFiles',
            condition='usePanelOfNormals', label='Panel of normals VCFs:'
        )
        form.addParam(
            'ponIndex', IntParam, default=0,
            condition='usePanelOfNormals', label='Panel of normals VCF index:',
            help='Zero-based index of the PoN in the VCF set.'
        )

        form.addSection(label='RNA-seq calling options')
        form.addParam(
            'disableMateOnSameContigFilter', BooleanParam, default=True,
            label='Disable same-contig mate filter:',
            help='Pass --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter, '
                 'as in the reference RNA-seq script.'
        )
        form.addParam(
            'dontUseSoftClippedBases', BooleanParam, default=True,
            label='Exclude soft-clipped bases:',
            help='Pass --dont-use-soft-clipped-bases, as in the reference RNA-seq script.'
        )

        form.addSection(label='Variant filtering')
        form.addParam(
            'variantType', EnumParam, default=0,
            choices=self.VARIANT_TYPES, label='Variant types:',
            help='Select the variants reported in the output VCFs. '
                 'Selection is applied after calling and optional filtering. '
                 'SNVs only keeps GATK SNP records; indels only keeps INDEL '
                 'records. Mixed and multi-nucleotide substitution records '
                 'are excluded in either exclusive mode. SNVs and indels '
                 'retains all Mutect2 calls.'
        )
        form.addParam(
            'filterCalls', BooleanParam, default=False,
            label='Run FilterMutectCalls:',
            help='Optional filtering, disabled by default to match the reference script. '
                 'Apply GATK FilterMutectCalls to Mutect2 calls. '
                 'Variants failing filters remain in the VCF with FILTER tags.'
        )
        form.addParallelSection(threads=4, mpi=1)

    def _insertAllSteps(self):
        self._insertFunctionStep(self.prepareInputsStep)
        self._insertFunctionStep(self.mutect2Step)
        if self.filterCalls.get():
            self._insertFunctionStep(self.filterMutectCallsStep)
        if self.variantType.get() != 0:
            self._insertFunctionStep(self.selectVariantTypesStep)
        self._insertFunctionStep(self.createOutputStep)

    # -------------------------------------------------------
    # Input preparation
    # -------------------------------------------------------

    def prepareInputsStep(self):
        errors = self._validate()
        if errors:
            raise RuntimeError('Invalid Mutect2 inputs: ' + '; '.join(errors))
        os.makedirs(self._getExtraPath(), exist_ok=True)
        with open(self._getCommandFile(), 'w') as commandFile:
            commandFile.write('')
        self.removeIndex(self._getSamplesFile())
        tumor = self.inputTumor.get()
        contigs = self.prepareReference(self.getReferenceFasta(tumor))
        tumorSamples = self.prepareAlignment(tumor, self._getAlignmentPath('tumor'), contigs,
                                              'Tumor BAM', requireReadGroups=True)
        if len(tumorSamples) != 1:
            raise RuntimeError('Tumor BAM must contain exactly one SM sample.')
        samples = {'tumor': next(iter(tumorSamples)), 'normal': None}
        if self.useNormal.get():
            samples['normal'] = self._prepareNormalSample(contigs, samples['tumor'])

        for enabled, pointer, index, stem, option in self._resourceInputs():
            if enabled.get():
                vcf = self._getSelectedResource(pointer, index)
                path = self.getLocalVCFPath(vcf.getFileName(), stem)
                self.prepareVCF(vcf, path, contigs, stem,
                                 validateHeader=self._validateGermlineAF if stem == 'germline' else None)
        filename = self._getSamplesFile()
        with open(filename + '.part', 'w') as handle:
            json.dump(samples, handle, indent=2, sort_keys=True)
        os.replace(filename + '.part', filename)

    def _prepareNormalSample(self, contigs, tumorSample):
        normalSamples = self.prepareAlignment(
            self.inputNormal.get(), self._getAlignmentPath('normal'), contigs,
            'Normal BAM', requireReadGroups=True
        )
        if len(normalSamples) != 1:
            raise RuntimeError('Normal BAM must contain exactly one SM sample.')
        normalSample = next(iter(normalSamples))
        requested = (self.normalSample.get() or '').strip()
        if requested and requested != normalSample:
            raise RuntimeError('Normal sample {} does not match BAM SM {}.'.format(requested, normalSample))
        if normalSample == tumorSample:
            raise RuntimeError('Tumor and normal BAMs must have distinct SM sample names.')
        return normalSample

    def _resourceInputs(self):
        return [
            (self.useGermlineResource, self.germlineVCFs, self.germlineIndex,
             'germline', '--germline-resource'),
            (self.usePanelOfNormals, self.ponVCFs, self.ponIndex,
             'pon', '--panel-of-normals')
        ]

    @staticmethod
    def _getSelectedResource(pointer, indexParam):
        vcfSet = pointer.get()
        if vcfSet is None or vcfSet.getSize() == 0:
            raise RuntimeError('A non-empty SetOfVCFFiles is required.')
        index = indexParam.get()
        if index is None or index < 0 or index >= vcfSet.getSize():
            raise RuntimeError('VCF index {} is out of range (size {}).'.format(index, vcfSet.getSize()))
        for position, vcf in enumerate(vcfSet):
            if position == index:
                GATKPreparationMixin.checkNonempty(vcf.getFileName(), 'Selected resource VCF')
                return vcf
        raise RuntimeError('Unable to retrieve the selected resource VCF.')

    @staticmethod
    def _validateGermlineAF(path):
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            for line in handle:
                if line.startswith('##INFO=<') and re.search(r'[<,]ID=AF[,>]', line):
                    if not re.search(r',Number=A[,>]', line) or not re.search(r',Type=Float[,>]', line):
                        raise RuntimeError('Germline INFO/AF must declare Number=A and Type=Float.')
                    return
                if line.startswith('#CHROM'):
                    break
        raise RuntimeError('Germline resource must declare INFO/AF population allele frequencies; '
                           'an arbitrary known-sites VCF is insufficient.')

    # -------------------------------------------------------
    # Calling, filtering and type selection
    # -------------------------------------------------------

    def mutect2Step(self):
        self.clearVCFOutput(self._getUnfilteredVcf())
        self.removeIndex(self._getStatsFile())
        args = [
            'Mutect2', '-R', self.getPreparedReference(),
            '-I', self._getAlignmentPath('tumor'), '-O', self._getUnfilteredVcf(),
            '--native-pair-hmm-threads', self.numberOfThreads.get(),
            '--create-output-variant-index', 'true'
        ]
        if self.useNormal.get():
            args.extend(['-I', self._getAlignmentPath('normal'),
                         '-normal', self._readSamples()['normal']])
        if self.disableMateOnSameContigFilter.get():
            args.extend(['--disable-read-filter', 'MateOnSameContigOrNoMappedMateReadFilter'])
        if self.dontUseSoftClippedBases.get():
            args.append('--dont-use-soft-clipped-bases')
        for enabled, pointer, index, stem, option in self._resourceInputs():
            if enabled.get():
                vcf = self._getSelectedResource(pointer, index)
                args.extend([option, self.getLocalVCFPath(vcf.getFileName(), stem)])
        self._runGatk(args)
        self._checkOutputVCF(self._getUnfilteredVcf())
        if self.filterCalls.get():
            self.checkNonempty(self._getStatsFile(), 'Mutect2 statistics')

    def filterMutectCallsStep(self):
        self.checkNonempty(self._getStatsFile(), 'Mutect2 statistics')
        self.clearVCFOutput(self._getFilteredVcf())
        self._runGatk([
            'FilterMutectCalls', '-R', self.getPreparedReference(),
            '-V', self._getUnfilteredVcf(), '--stats', self._getStatsFile(),
            '-O', self._getFilteredVcf(), '--create-output-variant-index', 'true'
        ])
        self._checkOutputVCF(self._getFilteredVcf())

    def selectVariantTypesStep(self):
        """Select full VCF record types after optional filtering."""
        variantType = {1: 'SNP', 2: 'INDEL'}[self.variantType.get()]
        sources = [self._getUnfilteredVcf()]
        if self.filterCalls.get():
            sources.append(self._getFilteredVcf())
        for source in sources:
            target = self._getSelectedVcf(source)
            self.clearVCFOutput(target)
            self._runGatk([
                'SelectVariants', '-R', self.getPreparedReference(),
                '-V', source, '--select-type-to-include', variantType,
                '-O', target, '--create-output-variant-index', 'true'
            ])
            self._checkOutputVCF(target)

    def _getSelectedVcf(self, source):
        if self.variantType.get() == 0:
            return source
        suffix = {1: 'snvs', 2: 'indels'}[self.variantType.get()]
        return source[:-len('.vcf.gz')] + '_' + suffix + '.vcf.gz'

    def _checkOutputVCF(self, path):
        self.checkNonempty(path, 'Mutect2 VCF')
        self.readVCFContigs(path)
        with gzip.open(path, 'rt') as handle:
            for line in handle:
                if line.startswith('#CHROM\t'):
                    samples = line.rstrip().split('\t')[9:]
                    break
            else:
                raise RuntimeError('Mutect2 VCF has no sample header: {}'.format(path))
        expected = {sample for sample in self._readSamples().values() if sample is not None}
        if len(samples) != len(expected) or set(samples) != expected:
            raise RuntimeError('VCF samples {} do not match input samples {}.'.format(samples, sorted(expected)))
        self.ensureFeatureIndex(path)

    # -------------------------------------------------------
    # Outputs and metadata
    # -------------------------------------------------------

    def createOutputStep(self):
        unfiltered = self._makeVcf(self._getSelectedVcf(self._getUnfilteredVcf()), filteringApplied=False)
        filtered = (self._makeVcf(self._getSelectedVcf(self._getFilteredVcf()), filteringApplied=True)
                    if self.filterCalls.get() else unfiltered.clone())
        outputs = {'outputUnfilteredVCF': unfiltered, 'outputVCF': filtered}
        self._defineOutputs(**outputs)
        inputs = [self.inputTumor]
        if self.useNormal.get():
            inputs.append(self.inputNormal)
        inputs.extend(pointer for enabled, pointer, index, stem, option in self._resourceInputs() if enabled.get())
        for output in outputs.values():
            for pointer in inputs:
                self._defineSourceRelation(pointer, output)

    def _makeVcf(self, path, filteringApplied):
        self._checkOutputVCF(path)
        tumor = self.inputTumor.get()
        samples = self._readSamples()
        vcf = VCFFile(filename=path)
        vcf.setIndexFile(path + '.tbi')
        vcf.setReferenceFasta(self.getPreparedReference())
        vcf.setIsCompressed(True)
        vcf.setScientificName(tumor.getReferenceScientificName())
        vcf.setAssembly(tumor.getReferenceAssembly())
        vcf.setSource('GATK')
        vcf.setDatabase('Mutect2')
        vcf.setVariantType(['somatic candidates', 'SNP', 'INDEL'][self.variantType.get()])
        vcf.setTumorSample(samples['tumor'])
        vcf.setNormalSample(samples['normal'])
        vcf.setFilteringApplied(filteringApplied)
        return vcf

    # -------------------------------------------------------
    # Paths and command history
    # -------------------------------------------------------

    def _getAlignmentPath(self, label):
        return self._getExtraPath('inputs', label + '.bam')

    def _getSamplesFile(self):
        return self._getExtraPath('mutect2_samples.json')

    def _readSamples(self):
        with open(self._getSamplesFile()) as handle:
            return json.load(handle)

    def _getUnfilteredVcf(self):
        return self._getExtraPath('mutect2_unfiltered.vcf.gz')

    def _getFilteredVcf(self):
        return self._getExtraPath('mutect2_filtered.vcf.gz')

    def _getStatsFile(self):
        return self._getUnfilteredVcf() + '.stats'

    def _getCommandFile(self):
        return self._getExtraPath('gatk_mutect2_commands.txt')

    def _runGatk(self, arguments):
        self.runTool('gatk', arguments)

    # -------------------------------------------------------
    # Validation, summary and methods
    # -------------------------------------------------------

    def _validate(self):
        errors = []
        tumor = self.inputTumor.get()
        if tumor is None:
            return ['A tumor AlignmentFile is required.']
        alignments = [('Tumor', tumor)]
        if self.useNormal.get():
            normal = self.inputNormal.get()
            if normal is None:
                errors.append('A normal AlignmentFile is required.')
            else:
                alignments.append(('Normal', normal))
        self._validateMutectAlignments(alignments, errors)
        if len(alignments) == 2:
            self._validateMatchedPair(tumor, alignments[1][1], errors)
        self._validateResources(errors)
        if self.numberOfThreads.get() < 1:
            errors.append('At least one thread is required.')
        if self.variantType.get() not in range(len(self.VARIANT_TYPES)):
            errors.append('Select a valid variant type.')
        return errors

    def _validateMutectAlignments(self, alignments, errors):
        for label, alignment in alignments:
            try:
                self.checkNonempty(alignment.getFileName(), label + ' BAM')
                self.checkNonempty(self.getReferenceFasta(alignment), label + ' reference FASTA')
            except RuntimeError as error:
                errors.append(str(error))
            if alignment.getFormat() != 'BAM':
                errors.append(label + ' alignment must be BAM.')

    def _validateMatchedPair(self, tumor, normal, errors):
        tumorBam, normalBam = tumor.getFileName(), normal.getFileName()
        tumorRef, normalRef = self.getReferenceFasta(tumor), self.getReferenceFasta(normal)
        if all(path and os.path.isfile(path) for path in (tumorBam, normalBam)) and os.path.samefile(tumorBam, normalBam):
            errors.append('Tumor and normal must be different BAM files, including aliases.')
        if all(path and os.path.isfile(path) for path in (tumorRef, normalRef)) and not os.path.samefile(tumorRef, normalRef):
            errors.append('Tumor and normal must refer to the same reference FASTA file.')

    def _validateResources(self, errors):
        for enabled, pointer, index, stem, option in self._resourceInputs():
            if enabled.get():
                try:
                    self._getSelectedResource(pointer, index)
                except RuntimeError as error:
                    errors.append('{} resource: {}'.format(stem, error))

    def _summary(self):
        summary = [
            'Variant types: {}'.format(self.VARIANT_TYPES[self.variantType.get()]),
            'Calling mode: {}'.format('Tumor-normal' if self.useNormal.get() else 'Tumor-only'),
            'FilterMutectCalls: {}'.format('enabled' if self.filterCalls.get() else 'disabled'),
            'Same-contig mate filter disabled: {}'.format(self.disableMateOnSameContigFilter.get()),
            'Soft-clipped bases excluded: {}'.format(self.dontUseSoftClippedBases.get())
        ]
        if hasattr(self, 'outputVCF'):
            summary.append('Tumor sample: {}'.format(self.outputVCF.getTumorSample()))
            if self.outputVCF.getNormalSample():
                summary.append('Normal sample: {}'.format(self.outputVCF.getNormalSample()))
            summary.append('Output VCF: {}'.format(self.outputVCF.getFileName()))
        return summary

    def _methods(self):
        if not hasattr(self, 'outputVCF'):
            return []
        return [
            'Candidate somatic variants were called with GATK Mutect2 using {} analysis{}. '
            'Variant selection: {}. The result was stored in an indexed VCF file.'
            .format('tumor-normal' if self.useNormal.get() else 'tumor-only',
                    ' and FilterMutectCalls' if self.filterCalls.get() else '',
                    self.VARIANT_TYPES[self.variantType.get()])
        ]
