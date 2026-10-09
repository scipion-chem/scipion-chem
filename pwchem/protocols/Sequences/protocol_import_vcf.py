# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import json
import os
import shlex

from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import VCFFile, SetOfVCFFiles


class ProtImportVCF(EMProtocol):
    """
    Import one or more known-variant VCF files.

    This protocol imports local Variant Call Format (VCF) files into
    Scipion and creates a ``SetOfVCFFiles`` that can be used as input for
    downstream genomic and variant-processing protocols.

    Multiple VCF files can be imported in a single protocol execution.
    Each imported file is represented by an independent ``VCFFile`` object
    containing the VCF path and the available variant metadata.

    Existing indexes (.tbi, .csi, .idx) can be imported. If no index is
    specified, the protocol looks for a neighboring index and, if needed,
    creates one for supported VCF formats.


    ---------------------------------------------------------------------
    Input parameters
    ---------------------------------------------------------------------

    vcfsData : str
        JSON configuration describing the VCF files to import.

        Each entry can contain the following fields:

        scientificName
            Scientific name of the organism associated with the variants.

            For example::

                Homo sapiens
                Mus musculus

        assembly
            Genome assembly associated with the VCF.

            For example::

                GRCh38
                GRCm39

        vcfFile
            Path to the local VCF file to import.

            This field is required.

        indexFile
            Optional path to an existing VCF index.

            Tabix (``.tbi``) and CSI (``.csi``) indexes are explicitly
            supported. Other index extensions are preserved when imported.

        The ``vcfsData`` parameter is normally configured through the
        corresponding Scipion wizard rather than edited manually.


    ---------------------------------------------------------------------
    Imported VCF
    ---------------------------------------------------------------------

    Each configured VCF file is copied into the protocol output directory.

    Imported files are renamed using the scientific name and, when
    available, the genome assembly.

    For example, a VCF configured as::

        scientificName = Mus musculus
        assembly = GRCm39

    is stored using a name based on::

        Mus_musculus_GRCm39_variants

    The original VCF extension is preserved. Compressed ``.vcf.gz`` files
    therefore remain compressed after import.

    Each imported VCF is represented by a ``VCFFile`` object with the
    following metadata:

    scientificName
        Scientific name provided during import, when available.

    assembly
        Genome assembly provided during import, when available.

    source
        Set to ``Imported``.

    database
        Set to ``Local``.

    variantType
        Set to ``known``.

    isCompressed
        Determined from whether the imported VCF filename ends in ``.gz``.


    ---------------------------------------------------------------------
    VCF index
    ---------------------------------------------------------------------

    An existing VCF index can optionally be imported together with the
    corresponding VCF.

    For Tabix indexes (``.tbi``), the imported index is stored as::

        <imported_vcf>.tbi

    For CSI indexes (``.csi``), the imported index is stored as::

        <imported_vcf>.csi

    Other index extensions are preserved and stored using the generated
    VCF base name.

    The resulting index path is stored in the corresponding ``VCFFile``
    object.

    If no index is supplied, the protocol checks for a neighboring index.
    If none exists, it generates a .tbi index for BGZF-compressed .vcf.gz
    using bcftools, or an .idx index for plain .vcf using GATK.
    The input VCF is never modified.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputVCFs : SetOfVCFFiles
        Set containing all imported VCF files.

        A ``SetOfVCFFiles`` is created independently of the number of VCF
        files configured for import.

        Each ``VCFFile`` contains:

        - Path to the copied VCF file.
        - Scientific name, when provided.
        - Genome assembly, when provided.
        - Source set to ``Imported``.
        - Database set to ``Local``.
        - Variant type set to ``known``.
        - Compression status.
        - VCF index path (imported or generated).

        Each object is labelled using the scientific name and genome
        assembly when available.


    ---------------------------------------------------------------------
    Workflow
    ---------------------------------------------------------------------

    The protocol performs the following steps for each configured VCF:

    1. Read the VCF configuration from ``vcfsData``.

    2. Generate an output filename using the scientific name and genome
       assembly when available.

    3. Copy the VCF file into the protocol output directory.

    4. Create a ``VCFFile`` object and assign the available metadata.

    5. Copy an existing index (explicit or adjacent) or generate one,
       then associate it with the ``VCFFile``.

    6. Add the resulting ``VCFFile`` to the output ``SetOfVCFFiles``.


    ---------------------------------------------------------------------
    Validation
    ---------------------------------------------------------------------

    Before execution, the protocol verifies that:

    - The VCF configuration can be parsed.
    - At least one VCF has been configured.
    - Every configured entry contains a VCF file.
    - Every VCF file exists.
    - Every specified index file exists.

    Scientific name, genome assembly and VCF index are optional.


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Multiple VCF files can be imported in a single protocol execution.

    - Input VCF files are copied into the protocol output directory.

    - The original VCF files are not modified.

    - Both compressed and uncompressed VCF files can be imported.

    - ``.vcf.gz`` extensions are preserved when generating the imported
      filename.

    - Imported VCF files are classified as known variants.

    - The source of imported VCF objects is stored as ``Imported`` and
      their database as ``Local``.

    - Scientific name and genome assembly metadata are optional.

    - VCF index files are optional inputs; missing indexes are generated.

    - Existing ``.tbi`` and ``.csi`` indexes are associated directly with
      the copied VCF using the standard VCF index naming convention.

    - Index generation requires bcftools (.vcf.gz) or GATK (.vcf) in
      the RNASEQ_DIC environment. BGZF VCFs must be sorted.

    - The protocol always returns a ``SetOfVCFFiles`` to provide a
      consistent output type for downstream Scipion protocols.
    """

    VCF_GZ_EXTENSION = '.vcf.gz'

    _label = 'import vcf'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'vcfsData',
            params.StringParam,
            default='[]',
            label='VCFs',
            help='Configure the VCF files to import.'
        )

    def _insertAllSteps(self):
        self._insertFunctionStep(
            self.createOutputStep
        )

    def _createVCF(self, vcfData):
        """Create a VCFFile object from imported VCF data."""
        scientificName = vcfData.get(
            'scientificName',
            ''
        )
        assembly = vcfData.get(
            'assembly',
            ''
        )
        vcfFile = vcfData.get(
            'vcfFile',
            ''
        )
        indexFile = vcfData.get(
            'indexFile',
            ''
        )

        safeName = (
            scientificName.strip().replace(' ', '_')
            if scientificName
            else 'variants'
        )

        if assembly:
            safeName += '_{}'.format(
                assembly.strip().replace(' ', '_')
            )

        if vcfFile.endswith(self.VCF_GZ_EXTENSION):
            extension = self.VCF_GZ_EXTENSION
        else:
            extension = os.path.splitext(
                vcfFile
            )[1]

        importedVCF = self._getExtraPath(
            '{}_variants{}'.format(
                safeName,
                extension
            )
        )

        copyFile(
            vcfFile,
            importedVCF
        )

        vcf = VCFFile(
            filename=importedVCF
        )

        if scientificName:
            vcf.setScientificName(
                scientificName
            )

        if assembly:
            vcf.setAssembly(
                assembly
            )

        vcf.setSource(
            'Imported'
        )

        vcf.setDatabase(
            'Local'
        )

        vcf.setVariantType(
            'known'
        )

        vcf.setIsCompressed(
            importedVCF.endswith('.gz')
        )

        importedIndex = self._ensureVCFIndex(
            vcfFile,
            importedVCF,
            indexFile
        )
        vcf.setIndexFile(importedIndex)

        label = scientificName or 'VCF'

        if assembly:
            label += ' ({})'.format(
                assembly
            )

        vcf.setObjLabel(
            label
        )

        return vcf

    def _findAdjacentVCFIndex(self, vcfFile):
        """Return an adjacent VCF index if available."""
        for suffix in ('.tbi', '.csi', '.idx'):
            candidate = vcfFile + suffix
            if os.path.isfile(candidate):
                return candidate

        return None

    def _ensureVCFIndex(self, sourceVCF, importedVCF, indexFile):
        """Copy a supplied/adjacent index, or create one for the copied VCF."""
        if not indexFile:
            indexFile = self._findAdjacentVCFIndex(sourceVCF)

        if indexFile:
            suffix = os.path.splitext(indexFile)[1].lower()

            if suffix in ('.tbi', '.csi', '.idx'):
                importedIndex = importedVCF + suffix
            else:
                importedIndex = os.path.join(
                    os.path.dirname(importedVCF),
                    os.path.basename(importedVCF) + suffix
                )

            copyFile(indexFile, importedIndex)
            return importedIndex

        # Reuse an existing index beside the imported VCF if present.
        existingIndex = self._findAdjacentVCFIndex(importedVCF)
        if existingIndex:
            return existingIndex

        if importedVCF.endswith(self.VCF_GZ_EXTENSION):
            # bcftools requires a sorted, BGZF-compressed VCF.
            args = 'index -f -t {}'.format(
                shlex.quote(importedVCF)
            )
            expectedIndex = importedVCF + '.tbi'
            program = 'bcftools'

        elif importedVCF.endswith('.vcf'):
            args = 'IndexFeatureFile -I {}'.format(
                shlex.quote(importedVCF)
            )
            expectedIndex = importedVCF + '.idx'
            program = 'gatk'

        else:
            raise RuntimeError(
                'Cannot index unsupported VCF format: {}'.format(
                    importedVCF
                )
            )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            program
        )

        if not os.path.isfile(expectedIndex):
            raise RuntimeError(
                'Index was not generated for {} (expected {}).'.format(
                    importedVCF,
                    expectedIndex
                )
            )

        return expectedIndex

    def createOutputStep(self):
        vcfsData = json.loads(
            self.vcfsData.get()
        )

        outputVCFs = SetOfVCFFiles().create(
            outputPath=self._getPath()
        )

        outputVCFs.setObjLabel(
            'Known variants VCFs'
        )

        for vcfData in vcfsData:
            outputVCFs.append(
                self._createVCF(
                    vcfData
                )
            )

        self._defineOutputs(
            outputVCFs=outputVCFs
        )

    def _validate(self):
        errors = []

        try:
            vcfsData = json.loads(
                self.vcfsData.get()
            )
        except (TypeError, ValueError):
            return [
                'Invalid VCF configuration.'
            ]

        if not vcfsData:
            errors.append(
                'At least one VCF must be configured.'
            )
            return errors

        for index, vcfData in enumerate(
            vcfsData,
            start=1
        ):
            vcfFile = vcfData.get(
                'vcfFile',
                ''
            )

            if not vcfFile:
                errors.append(
                    'VCF {}: VCF file is required.'.format(
                        index
                    )
                )

            elif not os.path.isfile(vcfFile):
                errors.append(
                    'VCF {}: VCF file does not exist: {}'.format(
                        index,
                        vcfFile
                    )
                )

            indexFile = vcfData.get(
                'indexFile',
                ''
            )

            if (
                indexFile
                and not os.path.isfile(indexFile)
            ):
                errors.append(
                    'VCF {}: index file does not exist: {}'.format(
                        index,
                        indexFile
                    )
                )

        return errors