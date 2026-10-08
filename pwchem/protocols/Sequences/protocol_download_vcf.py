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
import json
import os
import re
import shlex
import shutil
import urllib.error
import urllib.request
from urllib.parse import urlparse

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import (
    BooleanParam,
    EnumParam,
    StringParam
)

from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import VCFFile, SetOfVCFFiles
from pwchem.utils.utilsRNA import COMMON_SPECIES,getProviderSpeciesName,parseCustomSpecies


class ProtDownloadVCF(EMProtocol):
    """
    Download known-variant VCF files from Ensembl or NCBI/EVA.

    This protocol downloads one or more reference VCF files containing
    known genomic variants and creates a ``SetOfVCFFiles`` that can be
    used as input for downstream variant-processing protocols.

    Two reference sources can be selected:

    - Ensembl
    - NCBI

    When NCBI is selected, the protocol uses NCBI dbSNP when a suitable
    VCF is available. For species or assemblies not available through
    dbSNP, the European Variation Archive (EVA) is used as an alternative
    source.

    Species can be selected from a predefined list of commonly used
    organisms or entered manually using their scientific names. Multiple
    species can be processed in a single protocol execution.

    The protocol always produces a ``SetOfVCFFiles``, even when only one
    VCF file is successfully downloaded.


    ---------------------------------------------------------------------
    VCF sources
    ---------------------------------------------------------------------

    Ensembl
        Known-variant VCF files are downloaded from Ensembl Variation.

        The protocol first searches the main Ensembl release corresponding
        to the requested species and release.

        If no suitable VCF is available in the main Ensembl collection,
        the protocol automatically searches the Ensembl Genomes divisions:

        - plants
        - metazoa
        - fungi
        - protists
        - bacteria

        Ensembl Genomes uses its ``current`` release tree because its
        release numbering does not necessarily correspond to the release
        numbering of the main Ensembl site.

        When several VCF files are available, files containing consequence,
        phenotype or somatic-specific data are excluded when possible in
        order to select a general known-variation VCF.

        If a remote ``.tbi`` or ``.csi`` index is available, it is downloaded.
        Otherwise a local ``.csi`` index is generated.

    NCBI
        NCBI dbSNP is used as the primary source when a VCF is available
        for the selected assembly.

        Current NCBI dbSNP VCF releases are primarily available for human
        reference assemblies. Consequently, the protocol automatically
        uses the European Variation Archive (EVA) when an appropriate
        dbSNP VCF is not available.

        For non-human species requested with ``Latest``, the protocol
        directly searches EVA and selects the newest assembly for which
        EVA actually provides a current-ID VCF.

        When an explicit NCBI assembly accession is provided, the protocol
        first resolves its assembly metadata through NCBI. If dbSNP does
        not provide a VCF for that accession, EVA is searched for a
        compatible assembly.

    European Variation Archive (EVA)
        EVA is used as a fallback source for NCBI selections when dbSNP
        does not provide an appropriate VCF.

        The protocol downloads the ``*_current_ids.vcf.gz`` file associated
        with the selected species and assembly from the configured EVA
        release.

        If a remote index is available, it is downloaded. Otherwise
        a local ``.csi`` index is generated.

        For ``Latest`` non-human requests, the protocol inspects the
        assemblies actually available in EVA. Assembly accessions embedded
        in the VCF filenames are used to query their NCBI release dates,
        and the most recently released assembly is selected.

        If release dates cannot be resolved, the highest versioned assembly
        accession is used as a deterministic fallback.


    ---------------------------------------------------------------------
    Species selection
    ---------------------------------------------------------------------

    Species can be selected using one of two modes:

    Common species
        Select one or more species from the predefined common species
        list using the corresponding wizard.

        The selected values correspond to scientific species names.
        Provider-specific identifiers and default assemblies are resolved
        internally when required.

    Custom species
        Enter one or more species manually using their scientific names.

        Multiple species can be separated by semicolons or commas.

        Examples::

            Bos taurus

            Bos taurus;Oryza sativa

        Duplicate species entries are handled by the common species parsing
        utilities used by the genomic protocols.


    ---------------------------------------------------------------------
    Ensembl assembly selection
    ---------------------------------------------------------------------

    For Ensembl, the ``assemblies`` parameter determines the genome
    assembly associated with each VCF.

    ``Latest``
        Automatically resolves the assembly independently for each selected
        species.

        For species included in the common species configuration, the
        predefined reference assembly is used when available.

        Otherwise, the current assembly is obtained through the Ensembl
        REST API.

    Explicit assembly
        An assembly can be specified manually, for example::

            GRCh38

        A single assembly value is applied to every selected species.

        Alternatively, one assembly can be provided per species::

            GRCh38;GRCm39


    ---------------------------------------------------------------------
    Ensembl release selection
    ---------------------------------------------------------------------

    The ``releases`` parameter determines the Ensembl release used to
    locate the variation VCF.

    ``Latest``
        The most recent available Ensembl release is obtained through the
        Ensembl REST API.

    Explicit release
        A positive integer can be provided, for example::

            116

        A single release value is applied to every selected species.

        Alternatively, one release can be provided per species::

            116;116

    If the VCF is obtained from Ensembl Genomes rather than the main
    Ensembl collection, the stored release is ``current``.


    ---------------------------------------------------------------------
    NCBI assembly selection
    ---------------------------------------------------------------------

    For NCBI, the ``ncbiAssemblies`` parameter determines the requested
    reference assembly.

    ``Latest``
        Assembly resolution depends on the selected species.

        For Homo sapiens, the protocol resolves the current NCBI assembly,
        preferring RefSeq and using GenBank as a fallback when necessary.
        The corresponding dbSNP VCF is used when available.

        For non-human species, the protocol searches EVA directly and
        selects the newest assembly for which the configured EVA release
        provides a current-ID VCF.

    Explicit accession
        A versioned NCBI Assembly accession can be provided manually.

        Supported formats are::

            GCF_<accession>.<version>
            GCA_<accession>.<version>

        For example::

            GCF_000001405.40

        The assembly name is resolved using NCBI Datasets.

        If dbSNP does not provide a VCF corresponding to the requested
        accession, the protocol attempts to identify a compatible assembly
        in EVA.

        EVA matching first considers assembly names and exact accessions.
        When necessary, RefSeq and GenBank accessions sharing the same
        stable numeric assembly identifier can also be matched.

        A single accession can be applied to every selected species or one
        accession can be provided per species.


    ---------------------------------------------------------------------
    Downloaded files
    ---------------------------------------------------------------------

    VCF
        Downloaded variant files are stored as compressed VCF files::

            variants.vcf.gz

        Each successfully downloaded VCF is represented by a ``VCFFile``
        object.

    Index
        When the remote source provides an index, it is downloaded together
        with the VCF.

        Ensembl and NCBI dbSNP may provide::

            .tbi

        EVA may provide::

            .csi

        The index path is stored in the corresponding ``VCFFile`` object.

        When no remote index is available, bcftools generates a local
        ``.csi`` index (the VCF must be BGZF-compressed and sorted).


    ---------------------------------------------------------------------
    Existing files and overwrite behaviour
    ---------------------------------------------------------------------

    By default, previously downloaded files are reused when possible.

    If ``overwrite`` is enabled, existing destination files are removed
    before being downloaded again.

    This avoids unnecessary downloads while allowing users to explicitly
    refresh previously retrieved VCF files.


    ---------------------------------------------------------------------
    Input parameters
    ---------------------------------------------------------------------

    source : Enum
        Source requested for known-variant VCF retrieval.

        Available values:

        - Ensembl
        - NCBI

        Selecting NCBI may result in an EVA download when dbSNP does not
        provide an appropriate VCF.

    speciesSelection : Enum
        Method used to select species.

        Available values:

        - Common species
        - Custom species

    commonSpecies : str
        One or more scientific species names selected from the common
        species list using the wizard.

        Only available when ``speciesSelection`` is set to
        ``Common species``.

    customSpecies : str
        One or more manually specified scientific species names.

        Multiple species can be separated by semicolons or commas.

        Examples::

            Bos taurus

            Bos taurus;Oryza sativa

        Only available when ``speciesSelection`` is set to
        ``Custom species``.

    assemblies : str
        Ensembl assembly selection.

        Use ``Latest`` for automatic resolution or provide explicit
        assembly names.

        A single value is propagated to all selected species. Multiple
        values can be provided when one assembly is required per species.

        Only available when the selected source is Ensembl.

    releases : str
        Ensembl release selection.

        Use ``Latest`` for automatic resolution or provide positive
        integer release numbers.

        A single value is propagated to all selected species. Multiple
        values can be provided when one release is required per species.

        Only available when the selected source is Ensembl.

    ncbiAssemblies : str
        NCBI assembly selection.

        Use ``Latest`` for automatic source-specific resolution or provide
        a versioned ``GCF_`` or ``GCA_`` accession.

        For non-human ``Latest`` requests, the newest assembly actually
        available in EVA is selected.

        A single value is propagated to all selected species. Multiple
        values can be provided when one accession is required per species.

        Only available when the selected source is NCBI.

    overwrite : bool
        If True, existing downloaded files are replaced.


    ---------------------------------------------------------------------
    Output
    ---------------------------------------------------------------------

    outputVCFs : SetOfVCFFiles
        Set containing all successfully downloaded known-variant VCF files.

        A ``SetOfVCFFiles`` is always created independently of whether one
        or multiple species were selected.

        Each ``VCFFile`` contains, when available:

        - Scientific species name.
        - Genome assembly.
        - Effective data source.
        - Variant database.
        - Release identifier.
        - Variant type.
        - Compressed VCF file.
        - VCF index file.

        The effective source and database reflect where the VCF was
        actually obtained.

        Examples include:

        - ``Ensembl`` / ``Ensembl Variation``
        - ``NCBI`` / ``dbSNP``
        - ``EVA`` / ``European Variation Archive``

        The variant type is stored as ``known`` and downloaded VCF files
        are marked as compressed.


    ---------------------------------------------------------------------
    Workflow
    ---------------------------------------------------------------------

    The protocol is executed in three main steps:

    1. Resolve VCF information

       The selected species, assemblies and releases are interpreted and
       converted into source-specific metadata.

       For Ensembl, ``Latest`` releases and assemblies are resolved using
       Ensembl metadata.

       For NCBI, required assembly metadata are resolved using NCBI
       Datasets. Non-human ``Latest`` requests are deferred to EVA so that
       the newest assembly actually containing an EVA VCF can be selected.

    2. Download VCF files

       For Ensembl, the protocol first searches the main Ensembl Variation
       collection and then Ensembl Genomes when necessary.

       For NCBI, the protocol searches the current dbSNP VCF release.

       When dbSNP does not provide an appropriate VCF, or for non-human
       ``Latest`` requests, EVA is used.

       Available remote VCF indexes are downloaded together with their
       corresponding VCF files.

    3. Create output

       Successfully downloaded files and their metadata are converted into
       individual ``VCFFile`` objects and stored in the final
       ``SetOfVCFFiles``.

       Species for which no suitable VCF could be resolved or downloaded
       are skipped and are not included in the output set.


    ---------------------------------------------------------------------
    Requirements
    ---------------------------------------------------------------------

    Ensembl downloads require:

    - Internet access to the Ensembl REST API.
    - Internet access to the Ensembl FTP server.
    - Internet access to the Ensembl Genomes FTP server when fallback is
      required.

    NCBI downloads require:

    - Internet access to NCBI.
    - NCBI Datasets CLI (``datasets``).

    EVA downloads require:

    - Internet access to the European Variation Archive FTP server.

    Required command-line programs must be available through the RNA-seq
    environment configured by the Scipion-Chem plugin.


    ---------------------------------------------------------------------
    Notes
    ---------------------------------------------------------------------

    - Multiple species can be processed in a single protocol execution.

    - Species are specified using scientific names.

    - Species values can be separated by semicolons or commas.

    - A single assembly or release parameter is automatically propagated
      to all selected species.

    - When multiple assembly or release values are provided, the number
      of values must match the number of selected species.

    - Ensembl searches the main Variation collection first and Ensembl
      Genomes as a fallback.

    - Selecting NCBI does not guarantee that the resulting VCF originates
      from NCBI. EVA is automatically used when dbSNP does not provide a
      suitable VCF.

    - For non-human NCBI ``Latest`` requests, EVA determines the newest
      assembly that actually provides a current-ID VCF.

    - Explicit NCBI assembly accessions must be versioned ``GCF_`` or
      ``GCA_`` accessions.

    - EVA assembly matching can account for equivalent GCF/GCA accessions
      and version differences sharing the same stable assembly identifier.

    - Remote VCF indexes are downloaded when available. Otherwise a local
      ``.csi`` index is generated with bcftools.

    - VCF information is persisted internally in ``vcfs.json`` between
      protocol steps.

    - Errors affecting one species do not necessarily prevent other
      selected species from being processed.

    - Species for which no VCF is successfully downloaded are omitted from
      the final output set.

    - The protocol always returns a ``SetOfVCFFiles`` to provide a
      consistent output type for downstream Scipion protocols.
    """

    # ---------------------------------------------------------------------
    # Sources
    # ---------------------------------------------------------------------

    SOURCE_ENSEMBL = 0
    SOURCE_NCBI = 1

    SOURCE_CONDITION = 'source == %d'

    # ---------------------------------------------------------------------
    # Species selection
    # ---------------------------------------------------------------------

    SPECIES_COMMON = 0
    SPECIES_CUSTOM = 1

    # ---------------------------------------------------------------------
    # Ensembl REST
    # ---------------------------------------------------------------------

    ENSEMBL_REST_DATA_URL = (
        'https://rest.ensembl.org/info/data'
        '?content-type=application/json'
    )

    ENSEMBL_REST_ASSEMBLY_URL = (
        'https://rest.ensembl.org/info/assembly/{}'
        '?content-type=application/json'
    )

    ENSEMBL_FTP_BASE_URL = (
        'https://ftp.ensembl.org/pub/'
    )

    # Ensembl Genomes hosts non-vertebrate divisions such as plants.
    # The 'current' alias is used because Ensembl Genomes releases do not
    # necessarily have the same release number as the main Ensembl site.
    ENSEMBL_GENOMES_FTP_BASE_URL = (
        'https://ftp.ebi.ac.uk/ensemblgenomes/pub/'
    )

    ENSEMBL_GENOMES_DIVISIONS = (
        'plants',
        'metazoa',
        'fungi',
        'protists',
        'bacteria'
    )

    # ---------------------------------------------------------------------
    # NCBI dbSNP
    # ---------------------------------------------------------------------

    NCBI_DBSNP_VCF_URL = (
        'https://ftp.ncbi.nih.gov/snp/latest_release/VCF/'
    )

    # ---------------------------------------------------------------------
    # European Variation Archive (EVA)
    # ---------------------------------------------------------------------

    EVA_RELEASE = '9'

    EVA_BASE_URL = (
        'https://ftp.ebi.ac.uk/pub/databases/eva/rs_releases/'
        'release_{}/by_species/'.format(EVA_RELEASE)
    )

    # ---------------------------------------------------------------------
    # Allowed remote hosts
    # ---------------------------------------------------------------------

    ALLOWED_HOSTS = {
        'rest.ensembl.org',
        'ftp.ensembl.org',
        'ftp.ebi.ac.uk',
        'ftp.ncbi.nih.gov'
    }
    EMPTY_PARAMETER_ERROR = '{} cannot be empty.'
    ENSEMBL_SPECIES = 'Ensembl species'

    VCF_FILENAME = 'variants.vcf.gz'
    HREF_PATTERN = r'href="([^"]+)"'

    _label = 'download vcf'

    # =====================================================================
    # Parameters
    # =====================================================================

    def _defineParams(self, form):

        # ---------------------------------------------------------
        # Source
        # ---------------------------------------------------------

        form.addSection(label='Source')

        form.addParam(
            'source',
            EnumParam,
            choices=[
                'Ensembl',
                'NCBI'
            ],
            default=self.SOURCE_ENSEMBL,
            label='Source: ',
            help=(
                'Select the database used to download known variant VCF '
                'files.\n\n'
                'When NCBI is selected, dbSNP is used when a VCF is '
                'available. NCBI currently maintains dbSNP VCF releases '
                'for human assemblies; for non-human species, variant '
                'data are automatically retrieved from the European '
                'Variation Archive (EVA).'
            )
        )

        # ---------------------------------------------------------
        # Species selection
        # ---------------------------------------------------------

        form.addSection(label='Species selection')

        form.addParam(
            'speciesSelection',
            EnumParam,
            choices=[
                'Common species',
                'Custom species'
            ],
            default=self.SPECIES_COMMON,
            label='Species selection: ',
            help=(
                'Select species from the common species list or '
                'enter scientific species names manually.'
            )
        )

        form.addParam(
            'commonSpecies',
            StringParam,
            default='',
            condition=(
                'speciesSelection == %d'
                % self.SPECIES_COMMON
            ),
            label='Common species: ',
            help=(
                'Select one or more common species using the wizard.'
            )
        )

        form.addParam(
            'customSpecies',
            StringParam,
            default='',
            condition=(
                'speciesSelection == %d'
                % self.SPECIES_CUSTOM
            ),
            label='Species: ',
            help=(
                'Enter one or more scientific species names separated '
                'by semicolons.\n\n'
                'Examples:\n'
                'Bos taurus\n'
                'Bos taurus;Oryza sativa'
            )
        )

        # ---------------------------------------------------------
        # Ensembl
        # ---------------------------------------------------------

        form.addParam(
            'assemblies',
            StringParam,
            default='Latest',
            condition=(
                self.SOURCE_CONDITION
                % self.SOURCE_ENSEMBL
            ),
            label='Assemblies: ',
            help=(
                'Use "Latest" to resolve the current assembly for '
                'every selected species.\n\n'
                'Alternatively, provide one assembly per species.\n\n'
                'Examples:\n'
                'Latest\n'
                'GRCh38\n'
                'GRCh38;GRCm39'
            )
        )

        form.addParam(
            'releases',
            StringParam,
            default='Latest',
            condition=(
                self.SOURCE_CONDITION
                % self.SOURCE_ENSEMBL
            ),
            label='Ensembl releases: ',
            help=(
                'Use "Latest" to query the current Ensembl release.\n\n'
                'Alternatively, provide one release per species.\n\n'
                'Examples:\n'
                'Latest\n'
                '116\n'
                '116;116'
            )
        )

        # ---------------------------------------------------------
        # NCBI
        # ---------------------------------------------------------

        form.addParam(
            'ncbiAssemblies',
            StringParam,
            default='Latest',
            condition=(
                self.SOURCE_CONDITION
                % self.SOURCE_NCBI
            ),
            label='NCBI assemblies: ',
            help=(
                'Use "Latest" to resolve the current NCBI assembly for '
                'every selected species.\n\n'
                'Alternatively, provide a versioned NCBI assembly '
                'accession (GCF_ or GCA_).\n\n'
                'If dbSNP does not provide a VCF for the selected species '
                'and assembly, the protocol automatically searches EVA.\n\n'
                'Example:\n'
                'GCF_000001405.40'
            )
        )

        # ---------------------------------------------------------
        # Download
        # ---------------------------------------------------------

        form.addSection(label='Download options')

        form.addParam(
            'overwrite',
            BooleanParam,
            default=False,
            label='Overwrite existing files: '
        )

        form.addParallelSection(
            threads=1,
            mpi=0
        )

    # =====================================================================
    # Steps
    # =====================================================================

    def _insertAllSteps(self):

        self._insertFunctionStep(
            self.resolveVCFsStep
        )

        self._insertFunctionStep(
            self.downloadVCFsStep
        )

        self._insertFunctionStep(
            self.createOutputStep
        )

    # =====================================================================
    # Resolve
    # =====================================================================

    def resolveVCFsStep(self):

        vcfInfo = self._getSelectedVCFs()

        if self.source.get() == self.SOURCE_ENSEMBL:
            self._resolveEnsemblVCFs(vcfInfo)
        else:
            self._resolveNcbiVCFs(vcfInfo)

        self._writeVCFInfo(vcfInfo)

    # =====================================================================
    # Species
    # =====================================================================

    @staticmethod
    def _splitParameterValues(value):

        return [
            item.strip()
            for item in re.split(r'[;,]', value)
            if item.strip()
        ]

    def _getSelectedSpecies(self):

        if (
            self.speciesSelection.get()
            == self.SPECIES_COMMON
        ):
            selected = self._splitParameterValues(
                self.commonSpecies.get()
            )
        else:
            selected = parseCustomSpecies(
                self.customSpecies.get()
            )

        if not selected:
            raise ValueError(
                'At least one species must be selected.'
            )

        return selected

    # =====================================================================
    # VCF information
    # =====================================================================

    def _getSelectedVCFs(self):

        species = self._getSelectedSpecies()

        if self.source.get() == self.SOURCE_ENSEMBL:
            return self._getSelectedEnsemblVCFs(
                species
            )

        return self._getSelectedNcbiVCFs(
            species
        )

    # ---------------------------------------------------------------------
    # Ensembl information
    # ---------------------------------------------------------------------

    def _getSelectedEnsemblVCFs(
        self,
        species
    ):

        assemblies = self._expandParameterValues(
            self.assemblies.get(),
            len(species),
            'assemblies'
        )

        releases = self._expandParameterValues(
            self.releases.get(),
            len(species),
            'releases'
        )

        vcfInfo = []

        for speciesName, assembly, release in zip(
            species,
            assemblies,
            releases
        ):

            scientificName = speciesName.strip()

            ensemblName = getProviderSpeciesName(
                scientificName,
                'ensembl'
            )

            speciesInfo = COMMON_SPECIES.get(
                scientificName
            )

            defaultAssembly = (
                speciesInfo.get('assembly')
                if speciesInfo
                else None
            )

            vcfInfo.append({
                'scientificName': scientificName,
                'ensemblName': ensemblName,
                'defaultAssembly': defaultAssembly,
                'assembly': assembly,
                'release': release,
                'source': 'Ensembl',
                'database': 'Ensembl Variation',
                'variantType': 'known',
                'vcfFile': None,
                'indexFile': None,
                'error': None
            })

        return vcfInfo

    # ---------------------------------------------------------------------
    # NCBI information
    # ---------------------------------------------------------------------

    def _getSelectedNcbiVCFs(
        self,
        species
    ):

        assemblies = self._expandParameterValues(
            self.ncbiAssemblies.get(),
            len(species),
            'NCBI assemblies'
        )

        vcfInfo = []

        for speciesName, assemblyRequest in zip(
            species,
            assemblies
        ):

            scientificName = speciesName.strip()

            ncbiTaxon = getProviderSpeciesName(
                scientificName,
                'ncbi'
            )

            accession = None

            if assemblyRequest.lower() != 'latest':
                accession = assemblyRequest.upper()

            vcfInfo.append({
                'scientificName': scientificName,
                'ncbiTaxon': ncbiTaxon,
                'assemblyRequest': assemblyRequest,
                'accession': accession,
                'assembly': 'Latest',
                'release': 'latest',
                'requestedSource': 'NCBI',
                'source': None,
                'database': None,
                'variantType': 'known',
                'vcfFile': None,
                'indexFile': None,
                'error': None
            })

        return vcfInfo

    # ---------------------------------------------------------------------
    # Parameter utilities
    # ---------------------------------------------------------------------

    @classmethod
    def _expandParameterValues(
        cls,
        value,
        numberOfSpecies,
        parameterName
    ):

        values = cls._splitParameterValues(
            value
        )

        if not values:
            raise ValueError(
                cls.EMPTY_PARAMETER_ERROR.format(
                    parameterName.capitalize()
                )
            )

        if len(values) == 1:
            return values * numberOfSpecies

        if len(values) != numberOfSpecies:
            raise ValueError(
                'The number of {} must be either 1 or equal '
                'to the number of species.'.format(
                    parameterName
                )
            )

        return values

    # =====================================================================
    # Ensembl resolution
    # =====================================================================

    def _resolveEnsemblVCFs(
            self,
            vcfInfo
    ):

        latestRelease = None

        for info in vcfInfo:

            try:
                if info['release'].lower() == 'latest':

                    if latestRelease is None:
                        latestRelease = (
                            self._getLatestEnsemblRelease()
                        )

                    info['release'] = str(
                        latestRelease
                    )

                info['assembly'] = (
                    self._resolveEnsemblAssembly(
                        info['assembly'],
                        info
                    )
                )

            except RuntimeError as error:
                info['error'] = str(error)

                self.warning(
                    'Could not resolve Ensembl data for {}: {}'
                    .format(
                        info['scientificName'],
                        error
                    )
                )

    def _getLatestEnsemblRelease(self):

        data = self._requestJson(
            self.ENSEMBL_REST_DATA_URL,
            'latest Ensembl release'
        )

        releases = data.get(
            'releases',
            []
        )

        if not releases:
            raise RuntimeError(
                'The Ensembl REST response does not '
                'contain releases.'
            )

        return max(
            int(release)
            for release in releases
        )

    def _resolveEnsemblAssembly(
        self,
        assembly,
        info
    ):

        assembly = assembly.strip()

        if assembly.lower() != 'latest':
            return assembly

        defaultAssembly = info.get(
            'defaultAssembly'
        )

        if defaultAssembly:
            return defaultAssembly

        return self._getLatestEnsemblAssembly(
            info['ensemblName']
        )

    def _getLatestEnsemblAssembly(
        self,
        species
    ):

        species = self._validateUrlComponent(
            species,
            self.ENSEMBL_SPECIES
        )

        url = (
            self.ENSEMBL_REST_ASSEMBLY_URL
            .format(species)
        )

        data = self._requestJson(
            url,
            'latest assembly for {}'.format(
                species
            )
        )

        assembly = data.get(
            'assembly_name'
        )

        if not assembly:
            raise RuntimeError(
                'The Ensembl REST response does not contain '
                'an assembly name for {}.'.format(
                    species
                )
            )

        return assembly

    # =====================================================================
    # NCBI resolution
    # =====================================================================

    def _resolveNcbiVCFs(
            self,
            vcfInfo
    ):
        """Resolve assemblies needed by the NCBI/EVA source.

        For an explicit GCF_/GCA_ accession, resolve its assembly name through
        NCBI as before.

        For ``Latest``:
          * Homo sapiens is resolved through NCBI because dbSNP provides the
            human VCF there.
          * Non-human species are left unresolved here. They are handled
            directly by EVA during the download step, where the newest
            assembly actually available in the selected EVA release is chosen.
        """

        for info in vcfInfo:

            try:
                assemblyRequest = str(
                    info.get('assemblyRequest', 'Latest')
                ).strip()

                isLatest = (
                    assemblyRequest.lower() == 'latest'
                )

                isHuman = (
                    self._normaliseEvaName(
                        info['scientificName']
                    )
                    == 'homosapiens'
                )

                if isLatest and not isHuman:
                    info['accession'] = None
                    info['assembly'] = 'Latest'
                    continue

                if info.get('accession') is None:
                    info['accession'] = (
                        self._getNcbiReferenceAccession(
                            info['ncbiTaxon']
                        )
                    )

                info['assembly'] = (
                    self._getNcbiAssemblyName(
                        info['accession']
                    )
                )

            except RuntimeError as error:
                info['error'] = str(error)

                self.warning(
                    'Could not resolve NCBI data for {}: {}'
                    .format(
                        info['scientificName'],
                        error
                    )
                )

    def _getNcbiReferenceAccession(
        self,
        taxon
    ):
        """Return the newest NCBI assembly, preferring RefSeq.

        dbSNP VCF availability is checked later. Resolving a GenBank
        assembly here therefore does not imply that dbSNP provides a VCF
        for that assembly.
        """

        reports = self._getNcbiAssemblyReports(
            taxon,
            'RefSeq'
        )

        if not reports:
            self.info(
                'No RefSeq assembly found for {}. Trying GenBank.'
                .format(taxon)
            )

            reports = self._getNcbiAssemblyReports(
                taxon,
                'GenBank'
            )

        if not reports:
            raise RuntimeError(
                'NCBI does not provide a current RefSeq or GenBank '
                'assembly for {}.'.format(taxon)
            )

        def getReleaseDate(report):

            assemblyInfo = (
                report.get('assemblyInfo')
                or report.get('assembly_info')
                or {}
            )

            return (
                assemblyInfo.get('releaseDate')
                or assemblyInfo.get('release_date')
                or assemblyInfo.get('submissionDate')
                or assemblyInfo.get('submission_date')
                or ''
            )

        latestReport = max(
            reports,
            key=getReleaseDate
        )

        accession = (
            latestReport.get('accession')
            or latestReport.get('currentAccession')
            or latestReport.get('current_accession')
        )

        if not accession:
            raise RuntimeError(
                'The NCBI assembly report for {} does not '
                'contain an accession.'.format(taxon)
            )

        return accession

    def _getNcbiAssemblyReports(
        self,
        taxon,
        assemblySource
    ):
        """Retrieve NCBI assembly reports for one taxon and source."""

        outputFile = self._getTmpPath(
            '{}_{}_ncbi_summary.jsonl'.format(
                self._safeName(taxon),
                assemblySource.lower()
            )
        )

        arguments = (
            'summary genome taxon "{}" '
            '--assembly-source {} '
            '--tax-exact-match '
            '--as-json-lines '
            '> "{}"'
        ).format(
            taxon,
            assemblySource,
            outputFile
        )

        Plugin.runCondaCommand(
            self,
            arguments,
            RNASEQ_DIC,
            'datasets'
        )

        return self._readJsonLines(
            outputFile
        )

    def _getNcbiAssemblyName(
        self,
        accession
    ):

        outputFile = self._getTmpPath(
            '{}_assembly.jsonl'.format(
                self._safeName(accession)
            )
        )

        arguments = (
            'summary genome accession "{}" '
            '--as-json-lines '
            '> "{}"'
        ).format(
            accession,
            outputFile
        )

        Plugin.runCondaCommand(
            self,
            arguments,
            RNASEQ_DIC,
            'datasets'
        )

        reports = self._readJsonLines(
            outputFile
        )

        if not reports:
            raise RuntimeError(
                'Could not retrieve NCBI assembly '
                'information for {}.'.format(
                    accession
                )
            )

        report = reports[0]

        assemblyInfo = (
            report.get('assemblyInfo')
            or report.get('assembly_info')
            or {}
        )

        return (
            assemblyInfo.get('assemblyName')
            or assemblyInfo.get('assembly_name')
            or report.get('assemblyName')
            or report.get('assembly_name')
            or accession
        )

    # =====================================================================
    # Download
    # =====================================================================

    def downloadVCFsStep(self):

        vcfInfo = self._readVCFInfo()

        if self.source.get() == self.SOURCE_ENSEMBL:
            self._downloadEnsemblVCFs(
                vcfInfo
            )
        else:
            self._downloadNcbiVCFs(
                vcfInfo
            )

        self._writeVCFInfo(
            vcfInfo
        )

    # =====================================================================
    # Ensembl download
    # =====================================================================

    def _downloadEnsemblVCFs(
            self,
            vcfInfo
    ):

        for info in vcfInfo:

            if info.get('error'):
                continue

            try:
                filename, baseUrl, provider = self._findEnsemblVCF(
                    info
                )

                filename = self._validateUrlComponent(
                    filename,
                    'Ensembl VCF filename'
                )

                info['database'] = provider

                if provider == 'Ensembl Genomes':
                    info['release'] = 'current'

                ensemblName = self._validateUrlComponent(
                    info['ensemblName'],
                    self.ENSEMBL_SPECIES
                )

                assembly = self._validateUrlComponent(
                    info['assembly'],
                    'assembly'
                )

                outputDir = self._getExtraPath(
                    '{}_{}_release-{}'.format(
                        ensemblName,
                        assembly,
                        info['release']
                    )
                )

                os.makedirs(
                    outputDir,
                    exist_ok=True
                )

                vcfUrl = baseUrl + filename

                outputVCF = os.path.join(
                    outputDir,
                    self.VCF_FILENAME
                )

                self._downloadFile(
                    vcfUrl,
                    outputVCF
                )

                indexFile = self._ensureVCFIndex(
                    outputVCF,
                    vcfUrl
                )

                info['vcfFile'] = outputVCF
                info['indexFile'] = indexFile

            except (
                    RuntimeError,
                    ValueError
            ) as error:

                info['error'] = str(error)

                self.warning(
                    'Could not download Ensembl VCF for {}: {}'
                    .format(
                        info['scientificName'],
                        error
                    )
                )

    def _findEnsemblVCF(
        self,
        info
    ):
        """Find a VCF in Ensembl main, then in Ensembl Genomes.

        Returns
        -------
        tuple
            ``(filename, base_url, provider)``.
        """

        ensemblName = self._validateUrlComponent(
            info['ensemblName'],
            self.ENSEMBL_SPECIES
        )

        release = self._validateUrlComponent(
            info['release'],
            'Ensembl release'
        )

        mainUrl = (
            '{}release-{}/variation/vcf/{}/'
        ).format(
            self.ENSEMBL_FTP_BASE_URL,
            release,
            ensemblName
        )

        try:
            filename = self._findVCFInRemoteDirectory(
                mainUrl,
                info['scientificName']
            )
            return filename, mainUrl, 'Ensembl Variation'

        except RuntimeError as mainError:
            self.info(
                'No usable VCF found for {} in Ensembl main. '
                'Trying Ensembl Genomes.'.format(
                    info['scientificName']
                )
            )

        genomesResult = self._findEnsemblGenomesVCF(
            info
        )

        if genomesResult:
            filename, baseUrl, division = genomesResult
            self.info(
                'Using Ensembl Genomes {} VCF for {}: {}'
                .format(
                    division,
                    info['scientificName'],
                    filename
                )
            )
            return filename, baseUrl, 'Ensembl Genomes'

        raise RuntimeError(
            'No VCF file was found for {} in Ensembl main or '
            'Ensembl Genomes.'.format(
                info['scientificName']
            )
        ) from mainError

    def _findEnsemblGenomesVCF(
        self,
        info
    ):
        """Search the Ensembl Genomes divisions for a species VCF."""

        ensemblName = self._validateUrlComponent(
            info['ensemblName'],
            self.ENSEMBL_SPECIES
        )

        for division in self.ENSEMBL_GENOMES_DIVISIONS:

            division = self._validateUrlComponent(
                division,
                'Ensembl Genomes division'
            )

            baseUrl = (
                '{}{}/current/variation/vcf/{}/'
            ).format(
                self.ENSEMBL_GENOMES_FTP_BASE_URL,
                division,
                ensemblName
            )

            try:
                filename = self._findVCFInRemoteDirectory(
                    baseUrl,
                    info['scientificName']
                )
                return filename, baseUrl, division

            except RuntimeError:
                continue

        return None

    def _findVCFInRemoteDirectory(
        self,
        url,
        scientificName
    ):
        """Select a suitable VCF from a remote directory listing."""

        html = self._readRemoteDirectory(
            url
        )

        matches = re.findall(
            self.HREF_PATTERN,
            html
        )

        candidates = [
            filename
            for filename in matches
            if filename.endswith('.vcf.gz')
        ]

        if not candidates:
            raise RuntimeError(
                'No VCF files found in remote directory: {}'
                .format(url)
            )

        preferred = [
            filename
            for filename in candidates
            if 'incl_consequences' not in filename.lower()
            and 'phenotype' not in filename.lower()
            and 'somatic' not in filename.lower()
        ]

        if preferred:
            candidates = preferred

        candidates = sorted(
            candidates
        )

        if len(candidates) > 1:
            self.info(
                'Multiple VCF files found for {}. Using {}.'
                .format(
                    scientificName,
                    candidates[0]
                )
            )

        return self._validateUrlComponent(
            candidates[0],
            'Ensembl VCF filename'
        )

    # =====================================================================
    # NCBI dbSNP / EVA download
    # =====================================================================

    def _downloadNcbiVCFs(
            self,
            vcfInfo
    ):
        """Download from NCBI dbSNP, using EVA for non-human Latest requests."""

        try:
            directoryHtml = self._readRemoteDirectory(
                self.NCBI_DBSNP_VCF_URL
            )

            availableFiles = set(
                re.findall(
                    self.HREF_PATTERN,
                    directoryHtml
                )
            )

        except RuntimeError as error:
            self.warning(
                'Could not read the NCBI dbSNP VCF directory: {}. '
                'Trying EVA for the selected species.'.format(error)
            )
            availableFiles = set()

        for info in vcfInfo:

            if info.get('error'):
                continue

            try:
                assemblyRequest = str(
                    info.get('assemblyRequest', 'Latest')
                ).strip()

                isLatest = (
                    assemblyRequest.lower() == 'latest'
                )

                isHuman = (
                    self._normaliseEvaName(
                        info['scientificName']
                    )
                    == 'homosapiens'
                )

                # NCBI no longer maintains non-human dbSNP VCF releases.
                # If the user requested Latest, let EVA choose the newest
                # assembly that EVA itself actually provides.
                if isLatest and not isHuman:
                    self.info(
                        'NCBI does not maintain a current dbSNP VCF for '
                        '{}. Selecting the latest assembly available in EVA.'
                        .format(
                            info['scientificName']
                        )
                    )

                    self._downloadEvaVCF(
                        info
                    )
                    continue

                accession = self._validateUrlComponent(
                    info['accession'],
                    'NCBI accession'
                )

                filename = '{}.gz'.format(
                    accession
                )

                indexName = filename + '.tbi'

                if filename in availableFiles:
                    self._downloadNcbiVCF(
                        info,
                        accession,
                        filename,
                        indexName,
                        availableFiles
                    )
                    continue

                self.info(
                    'NCBI dbSNP does not provide a VCF for {} '
                    '({}, {}). Trying EVA.'.format(
                        info['scientificName'],
                        info['assembly'],
                        accession
                    )
                )

                self._downloadEvaVCF(
                    info
                )

            except (
                    RuntimeError,
                    ValueError
            ) as error:

                info['error'] = str(error)

                self.warning(
                    'Could not download a reference VCF for {}: {}'
                    .format(
                        info['scientificName'],
                        error
                    )
                )

    def _downloadNcbiVCF(
            self,
            info,
            accession,
            filename,
            indexName,
            availableFiles
    ):
        """Download one VCF from the NCBI dbSNP latest release."""

        outputDir = self._getExtraPath(
            'ncbi_{}_{}'.format(
                self._safeName(
                    info['scientificName']
                ),
                accession
            )
        )

        os.makedirs(
            outputDir,
            exist_ok=True
        )

        outputVCF = os.path.join(
            outputDir,
             self.VCF_FILENAME
        )

        self._downloadFile(
            self.NCBI_DBSNP_VCF_URL
            + filename,
            outputVCF
        )

        indexFile = self._ensureVCFIndex(
            outputVCF,
            self.NCBI_DBSNP_VCF_URL + filename,
            availableFiles=availableFiles,
            remoteFilename=filename
        )

        info['requestedSource'] = 'NCBI'
        info['source'] = 'NCBI'
        info['database'] = 'dbSNP'
        info['release'] = 'latest'
        info['vcfFile'] = outputVCF
        info['indexFile'] = indexFile
        info['error'] = None

    # =====================================================================
    # EVA fallback
    # =====================================================================

    def _downloadEvaVCF(
            self,
            info
    ):
        """Download the matching current-ID VCF from EVA."""

        filename, indexName, baseUrl, evaAssembly = (
            self._findEvaVCF(
                info
            )
        )

        outputDir = self._getExtraPath(
            'eva_{}_{}_release-{}'.format(
                self._safeName(
                    info['scientificName']
                ),
                self._safeName(
                    evaAssembly
                ),
                self.EVA_RELEASE
            )
        )

        os.makedirs(
            outputDir,
            exist_ok=True
        )

        outputVCF = os.path.join(
            outputDir,
             self.VCF_FILENAME
        )

        self._downloadFile(
            baseUrl + filename,
            outputVCF
        )

        indexFile = self._ensureVCFIndex(
            outputVCF,
            baseUrl + filename,
            availableFiles={indexName} if indexName else set(),
            remoteFilename=filename
        )

        info['requestedSource'] = 'NCBI'
        info['source'] = 'EVA'
        info['database'] = 'European Variation Archive'
        info['release'] = self.EVA_RELEASE
        info['evaAssembly'] = evaAssembly
        info['vcfFile'] = outputVCF
        info['indexFile'] = indexFile
        info['error'] = None

        self.info(
            'Using EVA release {} for {} ({}, EVA assembly directory {}).'
            .format(
                self.EVA_RELEASE,
                info['scientificName'],
                info['assembly'],
                evaAssembly
            )
        )

    def _findEvaVCF(
            self,
            info
    ):
        """Find a suitable current-ID VCF in EVA.

        If the user requested ``Latest``, select the newest assembly among the
        assemblies that are actually present in EVA. For an explicit assembly
        request, preserve the previous exact/matching behaviour.
        """

        speciesDirectory = self._getEvaSpeciesDirectory(
            info['scientificName']
        )

        speciesUrl = (
            self.EVA_BASE_URL
            + speciesDirectory
            + '/'
        )

        html = self._readRemoteDirectory(
            speciesUrl
        )

        entries = re.findall(
            self.HREF_PATTERN,
            html
        )

        assemblyDirectories = [
            entry.rstrip('/')
            for entry in entries
            if entry.endswith('/')
            and not entry.startswith('/')
            and not entry.startswith('?')
        ]

        if not assemblyDirectories:
            raise RuntimeError(
                'EVA does not provide assembly directories for {} '
                'in release {}.'.format(
                    info['scientificName'],
                    self.EVA_RELEASE
                )
            )

        assemblyRequest = str(
            info.get('assemblyRequest', 'Latest')
        ).strip()

        if assemblyRequest.lower() == 'latest':
            evaAssembly = self._getLatestEvaAssembly(
                speciesUrl,
                assemblyDirectories,
                info
            )
        else:
            evaAssembly = self._matchEvaAssemblyDirectory(
                assemblyDirectories,
                info
            )

        assemblyUrl = (
            speciesUrl
            + self._validateUrlComponent(
                evaAssembly,
                'EVA assembly'
            )
            + '/'
        )

        assemblyHtml = self._readRemoteDirectory(
            assemblyUrl
        )

        files = set(
            re.findall(
                self.HREF_PATTERN,
                assemblyHtml
            )
        )

        currentVCFs = sorted(
            filename
            for filename in files
            if filename.endswith('_current_ids.vcf.gz')
        )

        if not currentVCFs:
            raise RuntimeError(
                'EVA does not provide a current-ID VCF for {} '
                'assembly {}.'.format(
                    info['scientificName'],
                    evaAssembly
                )
            )

        filename = currentVCFs[0]

        csiName = filename + '.csi'
        indexName = (
            csiName
            if csiName in files
            else None
        )

        # The EVA VCF filename contains the assembly accession, e.g.
        # 9913_GCA_002263795.2_current_ids.vcf.gz.
        evaAccession = self._getEvaAccessionFromFilename(
            filename
        )

        info['assembly'] = evaAssembly

        if evaAccession:
            info['accession'] = evaAccession

        return (
            filename,
            indexName,
            assemblyUrl,
            evaAssembly
        )

    def _getEvaAssemblyCandidate(self, speciesUrl, directory):
        """Return the EVA VCF candidate for an assembly directory."""
        validatedDirectory = self._validateUrlComponent(
            directory,
            'EVA assembly'
        )

        assemblyUrl = speciesUrl + validatedDirectory + '/'

        assemblyHtml = self._readRemoteDirectory(
            assemblyUrl
        )

        files = set(
            re.findall(
                self.HREF_PATTERN,
                assemblyHtml
            )
        )

        currentVCFs = sorted(
            filename
            for filename in files
            if filename.endswith('_current_ids.vcf.gz')
        )

        if not currentVCFs:
            return None

        filename = currentVCFs[0]

        accession = self._getEvaAccessionFromFilename(
            filename
        )

        return {
            'directory': directory,
            'filename': filename,
            'accession': accession,
            'releaseDate': ''
        }

    def _setEvaCandidateReleaseDates(self, candidates):
        """Add NCBI release dates to EVA assembly candidates."""

        for candidate in candidates.values():
            accession = candidate.get('accession')

            if not accession:
                continue

            try:
                candidate['releaseDate'] = (
                    self._getNcbiAssemblyReleaseDate(
                        accession
                    )
                )
            except RuntimeError:
                candidate['releaseDate'] = ''

    def _getLatestEvaAssembly(
            self,
            speciesUrl,
            directories,
            info
    ):
        """Return the newest assembly that actually has a VCF in EVA.

        EVA exposes several aliases for the same assembly (assembly names and
        GCA accessions). Each candidate directory is inspected and the
        accession embedded in its ``*_current_ids.vcf.gz`` filename is used to
        obtain the NCBI assembly release date. The newest dated assembly wins.

        If release dates cannot be resolved, the candidate with the highest
        versioned accession is used as a deterministic fallback.
        """

        candidates = {}

        for directory in directories:
            try:
                candidate = self._getEvaAssemblyCandidate(
                    speciesUrl,
                    directory
                )
            except (RuntimeError, ValueError):
                continue

            if not candidate:
                continue

            key = candidate['accession'] or directory

            # EVA can expose both an assembly-name directory and a GCA
            # directory pointing to the same VCF. Keep only one.
            if key not in candidates:
                candidates[key] = candidate

        if not candidates:
            raise RuntimeError(
                'EVA does not provide a current-ID VCF for {} '
                'in release {}.'.format(
                    info['scientificName'],
                    self.EVA_RELEASE
                )
            )

        self._setEvaCandidateReleaseDates(candidates)

        datedCandidates = [
            candidate
            for candidate in candidates.values()
            if candidate.get('releaseDate')
        ]

        if datedCandidates:
            selected = max(
                datedCandidates,
                key=lambda candidate: (
                    candidate['releaseDate'],
                    self._evaAccessionSortKey(
                        candidate.get('accession')
                    )
                )
            )
        else:
            selected = max(
                candidates.values(),
                key=lambda candidate: (
                    self._evaAccessionSortKey(
                        candidate.get('accession')
                    ),
                    self._naturalSortKey(
                        candidate['directory']
                    )
                )
            )

        accession = selected.get('accession')
        releaseDate = selected.get('releaseDate')

        if accession and releaseDate:
            assemblyDetails = ' ({}, released {})'.format(
                accession,
                releaseDate
            )
        elif accession:
            assemblyDetails = ' ({})'.format(
                accession
            )
        else:
            assemblyDetails = ''

        self.info(
            'Latest EVA assembly selected for {}: {}{}.'
            .format(
                info['scientificName'],
                selected['directory'],
                assemblyDetails
            )
        )

        return selected['directory']

    @staticmethod
    def _getEvaAccessionFromFilename(
            filename
    ):
        """Extract a versioned GCA/GCF accession from an EVA VCF filename."""

        match = re.search(
            r'((?:GCA|GCF)_\d+\.\d+)',
            str(filename),
            re.IGNORECASE
        )

        return (
            match.group(1).upper()
            if match
            else None
        )

    def _getNcbiAssemblyReleaseDate(
            self,
            accession
    ):
        """Return the NCBI release/submission date for one assembly."""

        outputFile = self._getTmpPath(
            '{}_eva_assembly.jsonl'.format(
                self._safeName(accession)
            )
        )

        arguments = (
            'summary genome accession "{}" '
            '--as-json-lines '
            '> "{}"'
        ).format(
            accession,
            outputFile
        )

        Plugin.runCondaCommand(
            self,
            arguments,
            RNASEQ_DIC,
            'datasets'
        )

        reports = self._readJsonLines(
            outputFile
        )

        if not reports:
            raise RuntimeError(
                'Could not retrieve NCBI assembly information for {}.'
                .format(
                    accession
                )
            )

        report = reports[0]

        assemblyInfo = (
            report.get('assemblyInfo')
            or report.get('assembly_info')
            or {}
        )

        releaseDate = (
            assemblyInfo.get('releaseDate')
            or assemblyInfo.get('release_date')
            or assemblyInfo.get('submissionDate')
            or assemblyInfo.get('submission_date')
            or ''
        )

        return str(releaseDate)

    @staticmethod
    def _evaAccessionSortKey(
            accession
    ):
        """Create a deterministic sortable key for a GCA/GCF accession."""

        match = re.fullmatch(
            r'(?:GCA|GCF)_(\d+)\.(\d+)',
            str(accession or '').strip(),
            re.IGNORECASE
        )

        if not match:
            return (0, 0)

        return (
            int(match.group(1)),
            int(match.group(2))
        )

    @staticmethod
    def _naturalSortKey(
            value
    ):
        """Natural-sort strings containing numeric assembly versions."""

        return tuple(
            int(part) if part.isdigit() else part.lower()
            for part in re.split(
                r'(\d+)',
                str(value)
            )
        )

    def _getEvaSpeciesDirectory(
            self,
            scientificName
    ):
        """Return the EVA by_species directory for a scientific name."""

        expected = self._normaliseEvaName(
            scientificName
        )

        expectedDirectory = (
            scientificName.strip()
            .lower()
            .replace(' ', '_')
        )

        expectedDirectory = self._validateUrlComponent(
            expectedDirectory,
            'EVA species'
        )

        directUrl = (
            self.EVA_BASE_URL
            + expectedDirectory
            + '/'
        )

        try:
            self._readRemoteDirectory(
                directUrl
            )
            return expectedDirectory

        except RuntimeError:
            pass

        html = self._readRemoteDirectory(
            self.EVA_BASE_URL
        )

        directories = [
            entry.rstrip('/')
            for entry in re.findall(
                self.HREF_PATTERN,
                html
            )
            if entry.endswith('/')
            and not entry.startswith('/')
            and not entry.startswith('?')
        ]

        for directory in directories:
            if self._normaliseEvaName(directory) == expected:
                return self._validateUrlComponent(
                    directory,
                    'EVA species'
                )

        raise RuntimeError(
            'Species {} was not found in EVA release {}.'
            .format(
                scientificName,
                self.EVA_RELEASE
            )
        )

    def _matchEvaAssemblyName(self, directories, assemblyName):
        """Match an EVA directory using the assembly name."""
        if (
                not assemblyName
                or str(assemblyName).lower() == 'latest'
                or self._isNcbiAssemblyAccession(assemblyName)
        ):
            return None

        normalisedName = self._normaliseEvaName(
            assemblyName
        )

        for directory in directories:
            if (
                    self._normaliseEvaName(directory)
                    == normalisedName
            ):
                return directory

        return None

    def _matchEvaExactAccession(self, directories, candidates):
        """Match an EVA directory using an exact NCBI accession."""
        for candidate in candidates:
            candidate = str(candidate).strip()

            for directory in directories:
                if directory.lower() == candidate.lower():
                    return directory

        return None

    def _matchEvaAssemblyDirectory(
            self,
            directories,
            info
    ):
        """Match NCBI assembly metadata to an EVA assembly directory.

        EVA and NCBI do not always expose exactly the same assembly
        identifier. Prefer assembly-name matches, then exact accessions,
        and finally accessions sharing the same numeric assembly identifier
        while ignoring GCF/GCA and version differences.
        """

        assemblyName = info.get('assembly')
        accession = info.get('accession')
        assemblyRequest = info.get('assemblyRequest')

        # -------------------------------------------------------------
        # 1. Match by assembly name.
        # -------------------------------------------------------------

        directory = self._matchEvaAssemblyName(
            directories,
            assemblyName
        )

        if directory:
            return directory

        # -------------------------------------------------------------
        # 2. Try exact NCBI accession matches.
        # -------------------------------------------------------------

        accessionCandidates = [
            candidate
            for candidate in (
                accession,
                assemblyRequest,
                assemblyName
            )
            if candidate
               and str(candidate).lower() != 'latest'
               and self._isNcbiAssemblyAccession(
                candidate
            )
        ]

        directory = self._matchEvaExactAccession(
            directories,
            accessionCandidates
        )

        if directory:
            return directory

        # -------------------------------------------------------------
        # 3. Try normalised accession matches.
        # -------------------------------------------------------------

        normalisedCandidates = {
            self._normaliseEvaName(
                candidate
            )
            for candidate in accessionCandidates
        }

        for directory in directories:
            if (
                    self._normaliseEvaName(
                        directory
                    )
                    in normalisedCandidates
            ):
                return directory

        # -------------------------------------------------------------
        # 4. Match the stable numeric accession.
        #
        # NCBI RefSeq and GenBank accessions can describe the same
        # assembly using GCF_ and GCA_ respectively. EVA may also retain
        # another version of that assembly.
        #
        # Example:
        #
        # NCBI:
        #   GCF_002263795.3
        #
        # EVA:
        #   GCA_002263795.2
        #
        # Both have stable identifier:
        #   002263795
        # -------------------------------------------------------------

        accessionKeys = {
            self._getNcbiAssemblyAccessionKey(
                candidate
            )
            for candidate in accessionCandidates
        }

        accessionKeys.discard(
            None
        )

        if accessionKeys:
            for directory in directories:
                directoryKey = (
                    self._getNcbiAssemblyAccessionKey(
                        directory
                    )
                )

                if directoryKey in accessionKeys:
                    return directory

        raise RuntimeError(
            'EVA contains {}, but no assembly matching '
            '{} ({}) was found in release {}.'.format(
                info['scientificName'],
                info.get('assembly'),
                info.get('accession'),
                self.EVA_RELEASE
            )
        )

    @staticmethod
    def _isNcbiAssemblyAccession(
        value
    ):
        """Return whether a value is a versioned GCF/GCA accession."""

        return bool(
            re.fullmatch(
                r'(?:GCF|GCA)_\d+\.\d+',
                str(value).strip(),
                re.IGNORECASE
            )
        )

    @staticmethod
    def _getNcbiAssemblyAccessionKey(
        value
    ):
        """Return the stable numeric part of a GCF/GCA accession."""

        match = re.fullmatch(
            r'(?:GCF|GCA)_(\d+)\.\d+',
            str(value).strip(),
            re.IGNORECASE
        )

        return (
            match.group(1)
            if match
            else None
        )

    @staticmethod
    def _normaliseEvaName(
            value
    ):
        """Normalise species/assembly names for EVA directory matching."""

        return re.sub(
            r'[^a-z0-9]+',
            '',
            str(value).lower()
        )


    # =====================================================================
    # HTTP utilities
    # =====================================================================

    HTTP_TIMEOUT = 120
    DOWNLOAD_TIMEOUT = 300

    HTTP_HEADERS = {
        'User-Agent': 'Scipion-Chem'
    }

    JSON_HEADERS = {
        'Accept': 'application/json',
        'User-Agent': 'Scipion-Chem'
    }

    @classmethod
    def _validateUrlComponent(cls,value, componentName):
        """Validate a value before using it in a remote URL."""

        value = str(value).strip()

        if not value:
            raise ValueError(
                cls.EMPTY_PARAMETER_ERROR.format(
                    componentName
                )
            )

        if not re.fullmatch(
                r'[A-Za-z0-9_.-]+',
                value
        ):
            raise ValueError(
                'Invalid {}: {}'.format(
                    componentName,
                    value
                )
            )

        return value

    @classmethod
    def _validateRemoteUrl(cls, url):
        """Validate an HTTPS URL against the remote-host allowlist."""

        parsedUrl = urlparse(url)

        if parsedUrl.scheme != 'https':
            raise ValueError(
                'Only HTTPS URLs are allowed.'
            )

        if parsedUrl.hostname not in cls.ALLOWED_HOSTS:
            raise ValueError(
                'Remote host is not allowed: {}'
                .format(
                    parsedUrl.hostname
                )
            )

        if parsedUrl.username or parsedUrl.password:
            raise ValueError(
                'Credentials are not allowed in remote URLs.'
            )

        if parsedUrl.port not in (None, 443):
            raise ValueError(
                'Only the default HTTPS port is allowed.'
            )

        pathParts = [
            part
            for part in parsedUrl.path.split('/')
            if part
        ]

        if any(
                part in ('.', '..')
                for part in pathParts
        ):
            raise ValueError(
                'Path traversal is not allowed in remote URLs.'
            )

        return url

    @classmethod
    def _openRemoteRequest(
            cls,
            url,
            *,
            headers=None,
            method=None,
            timeout=None
    ):
        """Open a validated request to an allowed remote server."""

        validatedUrl = cls._validateRemoteUrl(url)

        request = urllib.request.Request(
            validatedUrl,
            headers=headers or cls.HTTP_HEADERS,
            method=method
        )

        return urllib.request.urlopen(  # NOSONAR
            request,
            timeout=timeout or cls.HTTP_TIMEOUT
        )

    def _downloadFile(
            self,
            url,
            outputFile
    ):
        """Download a file from an allowed remote server."""

        if (
                os.path.exists(outputFile)
                and not self.overwrite.get()
        ):
            return

        if (
                self.overwrite.get()
                and os.path.exists(outputFile)
        ):
            os.remove(
                outputFile
            )

        try:
            with self._openRemoteRequest(
                    url,
                    timeout=self.DOWNLOAD_TIMEOUT
            ) as response, open(
                outputFile,
                'wb'
            ) as output:
                shutil.copyfileobj(
                    response,
                    output
                )

        except (
                urllib.error.URLError,
                TimeoutError,
                ValueError
        ) as error:

            if os.path.exists(
                    outputFile
            ):
                os.remove(
                    outputFile
                )

            raise RuntimeError(
                'Could not download VCF file: {}'
                .format(
                    error
                )
            ) from error

        if not os.path.exists(
                outputFile
        ):
            raise RuntimeError(
                'Downloaded VCF file was not created: {}'
                .format(
                    outputFile
                )
            )

    @classmethod
    def _readRemoteDirectory(
            cls,
            url
    ):
        """Read an allowed remote directory."""

        try:
            with cls._openRemoteRequest(
                    url
            ) as response:

                return response.read().decode(
                    'utf-8'
                )

        except (
                urllib.error.URLError,
                TimeoutError,
                ValueError
        ) as error:

            raise RuntimeError(
                'Cannot access remote directory: {}'
                .format(
                    error
                )
            ) from error

    @classmethod
    def _remoteFileExists(
            cls,
            url
    ):
        """Return whether a remote file exists."""

        try:
            with cls._openRemoteRequest(
                    url,
                    method='HEAD'
            ):
                return True

        except (
                urllib.error.URLError,
                TimeoutError,
                ValueError
        ):
            return False

    def _ensureVCFIndex(
            self,
            vcfFile,
            vcfUrl,
            availableFiles=None,
            remoteFilename=None
    ):
        """Reuse/download a VCF index or create a CSI index locally.

        A remote TBI is preferred for compatibility with downstream tools.
        When the source has no index, bcftools creates a CSI index.
        """
        if remoteFilename is None:
            remoteFilename = vcfUrl.rsplit('/', 1)[-1]

        # An overwritten VCF must not retain an index from an older file.
        if self.overwrite.get():
            for extension in ('.tbi', '.csi'):
                indexPath = vcfFile + extension
                if os.path.isfile(indexPath):
                    os.remove(indexPath)
        else:
            for extension in ('.tbi', '.csi'):
                indexPath = vcfFile + extension
                if os.path.isfile(indexPath) and os.path.getsize(indexPath) > 0:
                    return indexPath

        for extension in ('.tbi', '.csi'):
            indexName = remoteFilename + extension
            indexUrl = vcfUrl + extension
            if availableFiles is not None:
                exists = indexName in availableFiles
            else:
                exists = self._remoteFileExists(indexUrl)

            if not exists:
                continue

            indexPath = vcfFile + extension
            try:
                self._downloadFile(indexUrl, indexPath)
                if os.path.isfile(indexPath) and os.path.getsize(indexPath) > 0:
                    return indexPath
            except RuntimeError as error:
                self.warning(
                    'Could not download VCF index {}: {}. Trying another '
                    'index or generating one locally.'.format(indexUrl, error)
                )

        self.info('No remote VCF index available. Generating CSI index.')
        arguments = 'index -f -c {}'.format(shlex.quote(vcfFile))
        try:
            Plugin.runCondaCommand(
                self,
                arguments,
                RNASEQ_DIC,
                'bcftools'
            )
        except Exception as error:
            raise RuntimeError(
                'Could not generate a VCF index for {}. Ensure bcftools '
                'is installed and the VCF is BGZF-compressed and sorted: {}'
                .format(vcfFile, error)
            ) from error

        indexPath = vcfFile + '.csi'
        if not os.path.isfile(indexPath) or os.path.getsize(indexPath) == 0:
            raise RuntimeError(
                'bcftools did not create the expected VCF index: {}'
                .format(indexPath)
            )
        return indexPath

    # =====================================================================
    # JSON utilities
    # =====================================================================

    @classmethod
    def _requestJson(
        cls,
        url,
        description
    ):

        try:
            with cls._openRemoteRequest(
                url,
                headers=cls.JSON_HEADERS
            ) as response:

                return json.loads(
                    response.read().decode(
                        'utf-8'
                    )
                )

        except (
            urllib.error.URLError,
            TimeoutError,
            ValueError
        ) as error:
            raise RuntimeError(
                'Could not retrieve {}: {}'
                .format(
                    description,
                    error
                )
            ) from error

    @staticmethod
    def _readJsonLines(
        filename
    ):

        if not os.path.exists(
            filename
        ):
            return []

        reports = []

        with open(
            filename
        ) as inputFile:

            for line in inputFile:

                line = line.strip()

                if line:
                    reports.append(
                        json.loads(
                            line
                        )
                    )

        return reports

    @staticmethod
    def _safeName(
        value
    ):

        return re.sub(
            r'[^A-Za-z0-9_.-]+',
            '_',
            value
        )

    # =====================================================================
    # Persistent VCF information
    # =====================================================================

    def _getVCFInfoFile(self):

        return self._getExtraPath(
            'vcfs.json'
        )

    def _writeVCFInfo(
        self,
        vcfInfo
    ):

        with open(
            self._getVCFInfoFile(),
            'w'
        ) as outputFile:

            json.dump(
                vcfInfo,
                outputFile,
                indent=2,
                sort_keys=True
            )

    def _readVCFInfo(self):

        with open(
            self._getVCFInfoFile()
        ) as inputFile:

            return json.load(
                inputFile
            )

    # =====================================================================
    # Output
    # =====================================================================

    def createOutputStep(self):

        vcfsInfo = self._readVCFInfo()

        outputVCFs = SetOfVCFFiles().create(
            outputPath=self._getPath()
        )

        outputVCFs.setObjLabel(
            'Known variants VCFs'
        )

        for info in vcfsInfo:

            if not info.get('vcfFile'):
                continue

            vcf = VCFFile(
                filename=info['vcfFile']
            )

            vcf.setScientificName(
                info['scientificName']
            )

            vcf.setAssembly(
                info['assembly']
            )

            vcf.setSource(
                info['source']
            )

            vcf.setDatabase(
                info['database']
            )

            vcf.setRelease(
                info['release']
            )

            vcf.setVariantType(
                info['variantType']
            )

            vcf.setIsCompressed(
                True
            )

            if info.get(
                    'indexFile'
            ):
                vcf.setIndexFile(
                    info['indexFile']
                )

            vcf.setObjLabel(
                '{} ({})'.format(
                    info['scientificName'],
                    info['assembly']
                )
            )

            outputVCFs.append(
                vcf
            )

        self._defineOutputs(
            outputVCFs=outputVCFs
        )

    # =====================================================================
    # Validation
    # =====================================================================

    def _validateEnsemblParams(self, species, errors):
        assemblies = self._splitParameterValues(
            self.assemblies.get()
        )

        releases = self._splitParameterValues(
            self.releases.get()
        )

        self._validateParameterCount(
            assemblies,
            species,
            'Assemblies',
            errors
        )

        self._validateParameterCount(
            releases,
            species,
            'Ensembl releases',
            errors
        )

        for release in releases:
            if release.lower() == 'latest':
                continue

            try:
                releaseValue = int(release)

                if releaseValue <= 0:
                    raise ValueError

            except ValueError:
                errors.append(
                    'Ensembl release must be "Latest" '
                    'or a positive integer: {}'.format(
                        release
                    )
                )

    def _validateNcbiParams(self, species, errors):
        assemblies = self._splitParameterValues(
            self.ncbiAssemblies.get()
        )

        self._validateParameterCount(
            assemblies,
            species,
            'NCBI assemblies',
            errors
        )

        accessionPattern = re.compile(
            r'^(GCF|GCA)_\d+\.\d+$',
            re.IGNORECASE
        )

        for assembly in assemblies:
            if assembly.lower() == 'latest':
                continue

            if not accessionPattern.match(assembly):
                errors.append(
                    'NCBI assembly must be "Latest" or '
                    'a valid versioned GCF_/GCA_ accession: {}'
                    .format(
                        assembly
                    )
                )

    @classmethod
    def _validateParameterCount(
            cls,
            values,
            species,
            parameterName,
            errors
    ):
        if not values:
            errors.append(
                cls.EMPTY_PARAMETER_ERROR.format(
                    parameterName
                )
            )
            return

        if (
                species
                and len(values) not in (
                1,
                len(species)
        )
        ):
            errors.append(
                '{} must contain one value or '
                'one value per species.'.format(
                    parameterName
                )
            )

    def _validate(self):
        errors = []

        try:
            species = self._getSelectedSpecies()
        except ValueError as error:
            errors.append(str(error))
            species = []

        if self.source.get() == self.SOURCE_ENSEMBL:
            self._validateEnsemblParams(species, errors)
        else:
            self._validateNcbiParams(species, errors)

        return errors


    # =====================================================================
    # Summary
    # =====================================================================

    def _summary(self):

        if not hasattr(
                self,
                'outputVCFs'
        ):
            return [
                'No VCF files have been downloaded.'
            ]

        summary = [
            'Downloaded VCF files: {}'.format(
                self.outputVCFs.getSize()
            )
        ]

        for vcf in self.outputVCFs:
            summary.append(
                '{} ({}) - source: {} ({})'.format(
                    vcf.getScientificName(),
                    vcf.getAssembly(),
                    vcf.getSource(),
                    vcf.getDatabase()
                )
            )

        vcfsInfo = self._readVCFInfo()

        skipped = [
            info
            for info in vcfsInfo
            if info.get('error')
        ]

        if skipped:

            summary.append(
                'Skipped species: {}'.format(
                    len(skipped)
                )
            )

            for info in skipped:
                summary.append(
                    '{}: {}'.format(
                        info['scientificName'],
                        info['error']
                    )
                )

        return summary
    # =====================================================================
    # Methods
    # =====================================================================

    def _methods(self):

        if not hasattr(
            self,
            'outputVCFs'
        ):
            return []

        sources = sorted({
            '{} ({})'.format(
                vcf.getSource(),
                vcf.getDatabase()
            )
            for vcf in self.outputVCFs
        })

        return [
            '{} known-variant VCF file(s) were downloaded '
            'from {}.'.format(
                self.outputVCFs.getSize(),
                ', '.join(sources)
            )
        ]