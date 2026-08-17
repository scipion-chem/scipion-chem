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
import gzip
import os
import re
import shlex
import shutil
import urllib.error
import urllib.request
import zipfile
from urllib.parse import urlparse, unquote

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import BooleanParam,EnumParam,StringParam
from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC
from pwchem.objects import Genome, SetOfGenomes


class ProtReferenceGenomes(EMProtocol):
    """
       Download reference genomes and genome annotations from Ensembl or NCBI.

       This protocol downloads one or more reference genomes together with
       their associated genome annotations and creates a ``SetOfGenomes``
       that can be used as input for downstream genomic and RNA-seq
       protocols.

       Two reference genome sources are supported:

       - Ensembl
       - NCBI

       The source is selected before defining the genomes to download. The
       available parameters for assembly and release selection depend on the
       selected source.

       Species can be selected from a predefined list of commonly used
       organisms or entered manually. Multiple genomes can be downloaded in
       a single protocol execution.

       The protocol always produces a ``SetOfGenomes``, even when only one
       reference genome is downloaded.


       ---------------------------------------------------------------------
       Genome source
       ---------------------------------------------------------------------

       Ensembl
           Reference genomes are downloaded directly from the Ensembl FTP
           server.

           The protocol resolves the requested Ensembl release and genome
           assembly before constructing the corresponding FASTA and GTF
           download URLs.

           Both the assembly and Ensembl release can be specified explicitly
           or resolved automatically using ``Latest``.

       NCBI
           Reference genomes are downloaded using the NCBI Datasets command
           line interface.

           Genome assemblies are identified by their versioned NCBI Assembly
           accession, for example::

               GCF_000001405.40

           or::

               GCA_000001405.29

           When ``Latest`` is requested, the protocol queries the available
           RefSeq assemblies for the selected taxon and selects the most
           recently released assembly according to the NCBI assembly
           metadata.

           The resolved accession is then used to download the genome data
           package.


       ---------------------------------------------------------------------
       Genome selection
       ---------------------------------------------------------------------

       Genomes can be selected using one of two modes:

       Common genomes
           Select one or more species from the predefined list of commonly
           used reference organisms.

           Currently supported common genomes include:

           - Homo sapiens
           - Mus musculus
           - Rattus norvegicus
           - Danio rerio
           - Drosophila melanogaster
           - Caenorhabditis elegans
           - Saccharomyces cerevisiae
           - Arabidopsis thaliana

       Custom genomes
           Enter one or more species manually.

           Multiple species can be separated by semicolons or commas.

           For Ensembl, species must be provided using Ensembl identifiers,
           for example::

               bos_taurus;oryza_sativa

           For NCBI, scientific names can be provided, for example::

               Bos taurus;Oryza sativa

           Duplicate species entries are automatically removed while
           preserving their original order.


       ---------------------------------------------------------------------
       Ensembl assembly selection
       ---------------------------------------------------------------------

       For Ensembl genomes, the ``assemblies`` parameter determines which
       genome assembly is downloaded.

       ``Latest``
           Automatically resolves the assembly.

           For genomes included in the common genome list, the predefined
           reference assembly is used.

           For custom genomes, the current assembly is obtained through the
           Ensembl REST API.

       Explicit assembly
           An assembly can be specified manually, for example::

               GRCh38

           When multiple genomes are selected, a single assembly value is
           applied to every genome.

           Alternatively, one assembly per genome can be provided::

               GRCh38;GRCm39


       ---------------------------------------------------------------------
       Ensembl release selection
       ---------------------------------------------------------------------

       The ``releases`` parameter determines the Ensembl release used for
       downloading the reference files.

       ``Latest``
           The most recent available Ensembl release is obtained through the
           Ensembl REST API.

       Explicit release
           A positive integer can be provided, for example::

               116

           A single release value can be applied to all selected genomes.

           Alternatively, one release can be specified per genome::

               116;116


       ---------------------------------------------------------------------
       NCBI assembly selection
       ---------------------------------------------------------------------

       For NCBI genomes, the ``ncbiAssemblies`` parameter determines the
       assembly to download.

       ``Latest``
           The protocol queries NCBI Datasets for RefSeq assemblies matching
           each selected taxon exactly.

           If several assemblies are returned for a taxon, the protocol
           compares their release dates and selects the most recently
           released assembly.

           The submission date is used as a fallback when a release date is
           not available.

           This resolution is performed independently for every selected
           species.

           The corresponding versioned NCBI Assembly accession is stored and
           subsequently used for downloading the genome.

       Explicit accession
           A specific versioned NCBI Assembly accession can be provided.

           Supported accession formats are::

               GCF_<accession>.<version>
               GCA_<accession>.<version>

           For example::

               GCF_000001405.40

           When multiple genomes are selected, a single value can be applied
           to every genome or one accession can be provided per genome::

               GCF_000001405.40;GCF_000001635.27


       ---------------------------------------------------------------------
       Genome FASTA
       ---------------------------------------------------------------------

       The genomic nucleotide sequence is always downloaded.

       For Ensembl, the protocol searches the corresponding release
       directory for a file ending in::

           .dna.toplevel.fa.gz

       The downloaded file is automatically decompressed before being stored
       in the output ``Genome`` object.

       For NCBI, the genome FASTA is retrieved from the NCBI Datasets package
       corresponding to the resolved assembly accession.


       ---------------------------------------------------------------------
       Genome annotation
       ---------------------------------------------------------------------

       Genome annotations can optionally be downloaded using the
       ``downloadAnnotation`` parameter.

       When enabled:

       Ensembl
           The protocol downloads the GTF annotation associated with the
           selected Ensembl release.

       NCBI
           The GTF annotation included in the NCBI genome data package is
           downloaded.

       The resulting GTF path is stored in the corresponding ``Genome``
       object.

       When annotation download is disabled, only the reference genome FASTA
       is downloaded.


       ---------------------------------------------------------------------
       Existing files and overwrite behaviour
       ---------------------------------------------------------------------

       By default, existing downloaded files are reused whenever possible.

       If ``overwrite`` is enabled, existing files and intermediate download
       directories are removed and the requested genome files are downloaded
       again.

       This behaviour allows interrupted or previously completed protocol
       executions to avoid unnecessary downloads.


       ---------------------------------------------------------------------
       Input parameters
       ---------------------------------------------------------------------

       source : Enum
           Database used to obtain the reference genomes.

           Available values:

           - Ensembl
           - NCBI

       genomeSelection : Enum
           Method used to select species.

           Available values:

           - Common genomes
           - Custom genomes

       commonGenomes : str
           One or more genomes selected from the predefined common genome
           list.

           Only available when ``genomeSelection`` is set to
           ``Common genomes``.

       customGenomes : str
           One or more manually specified species.

           Multiple species can be separated by semicolons or commas.

           Only available when ``genomeSelection`` is set to
           ``Custom genomes``.

       assemblies : str
           Ensembl assembly or assemblies.

           Use ``Latest`` for automatic resolution or provide explicit
           assembly names.

           A single value is propagated to all selected genomes. Multiple
           values can be provided when one assembly is required per genome.

           Only available when the selected source is Ensembl.

       releases : str
           Ensembl release or releases.

           Use ``Latest`` for automatic resolution or provide positive
           integer release numbers.

           A single value is propagated to all selected genomes. Multiple
           values can be provided when one release is required per genome.

           Only available when the selected source is Ensembl.

       ncbiAssemblies : str
           NCBI assembly selection.

           Use ``Latest`` to select the most recently released RefSeq
           assembly independently for each selected species, or provide a
           versioned ``GCF_`` or ``GCA_`` accession.

           A single value is propagated to all selected genomes. Multiple
           values can be provided when one accession is required per genome.

           Only available when the selected source is NCBI.

       downloadAnnotation : bool
           If True, download the GTF genome annotation in addition to the
           genome FASTA.

       overwrite : bool
           If True, existing genome files are removed and downloaded again.


       ---------------------------------------------------------------------
       Output
       ---------------------------------------------------------------------

       referenceGenomes : SetOfGenomes
           Set containing all downloaded reference genomes.

           A ``SetOfGenomes`` is always generated, independently of whether
           one or multiple genomes were selected.

           Each ``Genome`` object contains, when available:

           - Scientific name.
           - Genome assembly.
           - Data source.
           - Release or version identifier.
           - Genome FASTA file.
           - GTF annotation file.
           - Ensembl species identifier for Ensembl genomes.

           For Ensembl genomes, ``release`` corresponds to the Ensembl
           release number.

           For NCBI genomes, ``release`` stores the resolved versioned NCBI
           Assembly accession.


       ---------------------------------------------------------------------
       Workflow
       ---------------------------------------------------------------------

       The protocol is executed in three main steps:

       1. Resolve genomes

          The selected species, assemblies and releases are interpreted and
          converted into source-specific genome metadata.

          Each selected species is represented independently.

          For ``Latest`` selections, the corresponding remote database is
          queried to resolve the current release or assembly.

       2. Download genomes

          Genome FASTA files and, optionally, GTF annotations are downloaded
          independently for each selected species from Ensembl or NCBI.

          Compressed Ensembl files are automatically decompressed.

          NCBI data packages are extracted and the relevant files are copied
          to the protocol output directory.

       3. Create output

          Downloaded genome metadata and file paths are converted into
          individual ``Genome`` objects and stored in the final
          ``SetOfGenomes``.


       ---------------------------------------------------------------------
       Requirements
       ---------------------------------------------------------------------

       Ensembl downloads require:

       - Internet access to the Ensembl REST API.
       - Internet access to the Ensembl FTP server.
       - curl.
       - gunzip.

       NCBI downloads require:

       - Internet access to NCBI.
       - NCBI Datasets CLI (``datasets``).

       Required command-line programs must be available through the RNA-seq
       environment configured by the Scipion-Chem plugin.


       ---------------------------------------------------------------------
       Notes
       ---------------------------------------------------------------------

       - Multiple genomes can be downloaded in a single protocol execution.

       - Species values can be separated by semicolons or commas.

       - Duplicate species are removed while preserving their original order.

       - A single assembly or release parameter is automatically propagated
         to all selected genomes.

       - When multiple assembly or release values are provided, the number
         of values must match the number of selected genomes.

       - For NCBI ``Latest``, the most recently released RefSeq assembly is
         resolved independently for each selected taxon.

       - NCBI assembly accessions must include their version number.

       - Ensembl species identifiers use underscores, whereas NCBI custom
         species are specified using scientific names.

       - Genome information is persisted internally in ``genomes.json``
         between protocol steps.

       - Each selected species generates an independent ``Genome`` object.

       - The protocol always returns a ``SetOfGenomes`` to provide a
         consistent output type for downstream Scipion protocols.
       """

    # ---------------------------------------------------------------------
    # Sources
    # ---------------------------------------------------------------------

    SOURCE_ENSEMBL = 0
    SOURCE_NCBI = 1
    SOURCE_CONDITION = 'source == %d'

    # ---------------------------------------------------------------------
    # Genome selection
    # ---------------------------------------------------------------------

    GENOME_COMMON = 0
    GENOME_CUSTOM = 1

    # ---------------------------------------------------------------------
    # Ensembl
    # ---------------------------------------------------------------------

    ENSEMBL_REST_DATA_URL = (
        'https://rest.ensembl.org/info/data'
        '?content-type=application/json'
    )

    ENSEMBL_REST_ASSEMBLY_URL = (
        'https://rest.ensembl.org/info/assembly/{}'
        '?content-type=application/json'
    )

    # ---------------------------------------------------------------------
    # Common genomes
    # ---------------------------------------------------------------------

    COMMON_GENOMES = {
        'homo_sapiens': {
            'scientificName': 'Homo sapiens',
            'assembly': 'GRCh38'
        },
        'mus_musculus': {
            'scientificName': 'Mus musculus',
            'assembly': 'GRCm39'
        },
        'rattus_norvegicus': {
            'scientificName': 'Rattus norvegicus',
            'assembly': 'mRatBN7.2'
        },
        'danio_rerio': {
            'scientificName': 'Danio rerio',
            'assembly': 'GRCz11'
        },
        'drosophila_melanogaster': {
            'scientificName': 'Drosophila melanogaster',
            'assembly': 'BDGP6.46'
        },
        'caenorhabditis_elegans': {
            'scientificName': 'Caenorhabditis elegans',
            'assembly': 'WBcel235'
        },
        'saccharomyces_cerevisiae': {
            'scientificName': 'Saccharomyces cerevisiae',
            'assembly': 'R64-1-1'
        },
        'arabidopsis_thaliana': {
            'scientificName': 'Arabidopsis thaliana',
            'assembly': 'TAIR10'
        }
    }

    _label = 'reference genomes'

    # ---------------------------------------------------------------------
    # Parameters
    # ---------------------------------------------------------------------

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
                'Select the database used to download the reference '
                'genomes.'
            )
        )

        # ---------------------------------------------------------
        # Genome selection
        # ---------------------------------------------------------

        form.addSection(label='Genome selection')

        form.addParam(
            'genomeSelection',
            EnumParam,
            choices=[
                'Common genomes',
                'Custom genomes'
            ],
            default=self.GENOME_COMMON,
            label='Genome selection: ',
            help=(
                'Select genomes from the common species list or enter '
                'species manually.'
            )
        )

        form.addParam(
            'commonGenomes',
            StringParam,
            default='',
            condition='genomeSelection == %d' % self.GENOME_COMMON,
            label='Common genomes: ',
            help=(
                'Select one or more common genomes using the wizard.'
            )
        )

        form.addParam(
            'customGenomes',
            StringParam,
            default='',
            condition='genomeSelection == %d' % self.GENOME_CUSTOM,
            label='Species: ',
            help=(
                'Enter one or more species separated by semicolons.\n\n'
                'For Ensembl use species identifiers:\n'
                'bos_taurus;oryza_sativa\n\n'
                'For NCBI use scientific names:\n'
                'Bos taurus;Oryza sativa'
            )
        )

        # ---------------------------------------------------------
        # Ensembl-specific genome options
        # ---------------------------------------------------------

        form.addParam(
            'assemblies',
            StringParam,
            default='Latest',
            condition=self.SOURCE_CONDITION % self.SOURCE_ENSEMBL,
            label='Assemblies: ',
            help=(
                'Use "Latest" to resolve the current assembly for every '
                'selected species.\n\n'
                'A single value is applied to all genomes. Alternatively, '
                'provide one value per genome separated by semicolons.\n\n'
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
            condition=self.SOURCE_CONDITION % self.SOURCE_ENSEMBL,
            label='Ensembl releases: ',
            help=(
                'Use "Latest" to query the current Ensembl release.\n\n'
                'A single value is applied to all selected genomes. '
                'Alternatively, provide one release per genome separated '
                'by semicolons.\n\n'
                'Examples:\n'
                'Latest\n'
                '116\n'
                '116;116'
            )
        )

        # ---------------------------------------------------------
        # NCBI-specific genome options
        # ---------------------------------------------------------

        form.addParam(
            'ncbiAssemblies',
            StringParam,
            default='Latest',
            condition=self.SOURCE_CONDITION % self.SOURCE_NCBI,
            label='NCBI assemblies: ',
            help=(
                'Use "Latest" to download the most recently released '
                'RefSeq assembly for each selected species.\n\n'
                'Alternatively, provide a specific NCBI Assembly '
                'accession (GCF_ or GCA_).\n\n'
                'A single value is applied to every selected species or '
                'one value can be provided per genome.\n\n'
                'Examples:\n'
                'Latest\n'
                'GCF_000001405.40\n'
                'GCF_000001405.40;GCF_000001635.27'
            )
        )

        # ---------------------------------------------------------
        # Download options
        # ---------------------------------------------------------

        form.addSection(label='Download options')

        form.addParam(
            'downloadAnnotation',
            BooleanParam,
            default=True,
            label='Download annotation GTF: '
        )

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

    # ---------------------------------------------------------------------
    # Steps
    # ---------------------------------------------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep(self.resolveGenomesStep)
        self._insertFunctionStep(self.downloadGenomesStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------------------------------------------
    # Resolve genomes
    # ---------------------------------------------------------------------

    def resolveGenomesStep(self):
        genomes = self._getSelectedGenomes()

        if self.source.get() == self.SOURCE_ENSEMBL:
            self._resolveEnsemblGenomes(genomes)
        else:
            self._resolveNcbiGenomes(genomes)

        self._writeGenomeInfo(genomes)

    # ---------------------------------------------------------------------
    # Species selection
    # ---------------------------------------------------------------------

    @staticmethod
    def _splitParameterValues(value):
        """Split comma- or semicolon-separated parameter values."""
        return [
            item.strip()
            for item in re.split(r'[;,]', value)
            if item.strip()
        ]

    def _getSelectedSpecies(self):
        """Return selected species according to the selection mode."""

        value = (
            self.commonGenomes.get()
            if self.genomeSelection.get() == self.GENOME_COMMON
            else self.customGenomes.get()
        )

        selected = self._splitParameterValues(value)

        if (
                self.genomeSelection.get() == self.GENOME_COMMON
                or self.source.get() == self.SOURCE_ENSEMBL
        ):
            selected = [
                item.lower()
                for item in selected
            ]

        selected = list(dict.fromkeys(selected))

        if not selected:
            raise ValueError(
                'At least one genome must be selected.'
            )

        return selected

    # ---------------------------------------------------------------------
    # Genome information
    # ---------------------------------------------------------------------

    def _getSelectedGenomes(self):
        species = self._getSelectedSpecies()

        if self.source.get() == self.SOURCE_ENSEMBL:
            return self._getSelectedEnsemblGenomes(species)

        return self._getSelectedNcbiGenomes(species)

    # ---------------------------------------------------------------------
    # Ensembl genome information
    # ---------------------------------------------------------------------

    def _getSelectedEnsemblGenomes(self, species):

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

        genomes = []

        for speciesName, assembly, release in zip(
                species,
                assemblies,
                releases):

            if speciesName in self.COMMON_GENOMES:
                info = self.COMMON_GENOMES[speciesName]

                scientificName = info['scientificName']
                defaultAssembly = info['assembly']

            else:
                scientificName = (
                    self._speciesIdentifierToScientificName(
                        speciesName
                    )
                )
                defaultAssembly = None

            genomes.append({
                'scientificName': scientificName,
                'ensemblName': speciesName,
                'defaultAssembly': defaultAssembly,
                'assembly': assembly,
                'release': release,
                'source': 'Ensembl'
            })

        return genomes

    # ---------------------------------------------------------------------
    # NCBI genome information
    # ---------------------------------------------------------------------

    def _getSelectedNcbiGenomes(self, species):
        """Build NCBI genome information for the selected species."""

        assemblies = self._expandParameterValues(
            self.ncbiAssemblies.get(),
            len(species),
            'NCBI assemblies'
        )

        genomes = []

        for speciesName, assemblyRequest in zip(
                species,
                assemblies):

            commonKey = (
                speciesName
                .strip()
                .lower()
                .replace(' ', '_')
            )

            if commonKey in self.COMMON_GENOMES:
                scientificName = (
                    self.COMMON_GENOMES[
                        commonKey
                    ]['scientificName']
                )
            else:
                scientificName = speciesName.strip()

            if assemblyRequest.lower() == 'latest':
                accession = None
            else:
                accession = assemblyRequest.upper()

            genomes.append({
                'scientificName': scientificName,
                'ncbiTaxon': scientificName,
                'assemblyRequest': assemblyRequest,
                'accession': accession,
                'assembly': 'Latest',
                'release': None,
                'source': 'NCBI'
            })

        return genomes

    # ---------------------------------------------------------------------
    # Parameter utilities
    # ---------------------------------------------------------------------

    @classmethod
    def _expandParameterValues(
            cls,
            value,
            numberOfGenomes,
            parameterName):

        values = cls._splitParameterValues(value)

        if not values:
            raise ValueError(
                '{} cannot be empty.'.format(
                    parameterName.capitalize()
                )
            )

        if len(values) == 1:
            return values * numberOfGenomes

        if len(values) != numberOfGenomes:
            raise ValueError(
                'The number of {} must be either 1 or equal to the '
                'number of genomes.'.format(parameterName)
            )

        return values

    # =====================================================================
    # Ensembl
    # =====================================================================

    def _resolveEnsemblGenomes(self, genomes):

        latestRelease = None

        for genomeInfo in genomes:

            releaseValue = genomeInfo['release']

            if releaseValue.lower() == 'latest':

                if latestRelease is None:
                    latestRelease = (
                        self._getLatestEnsemblRelease()
                    )

                genomeInfo['release'] = str(latestRelease)

            genomeInfo['assembly'] = self._resolveAssembly(
                genomeInfo['assembly'],
                genomeInfo
            )

    def _getLatestEnsemblRelease(self):

        data = self._requestJson(
            self.ENSEMBL_REST_DATA_URL,
            'latest Ensembl release'
        )

        releases = data.get('releases', [])

        if not releases:
            raise RuntimeError(
                'The Ensembl REST response does not contain releases.'
            )

        return max(
            int(release)
            for release in releases
        )

    def _resolveAssembly(
            self,
            assembly,
            genomeInfo):

        assembly = assembly.strip()

        if assembly.lower() != 'latest':
            return assembly

        defaultAssembly = genomeInfo.get(
            'defaultAssembly'
        )

        if defaultAssembly:
            return defaultAssembly

        return self._getLatestEnsemblAssembly(
            genomeInfo['ensemblName']
        )

    def _getLatestEnsemblAssembly(self, species):

        url = self.ENSEMBL_REST_ASSEMBLY_URL.format(
            species
        )

        data = self._requestJson(
            url,
            'latest assembly for {}'.format(species)
        )

        assembly = data.get('assembly_name')

        if not assembly:
            raise RuntimeError(
                'The Ensembl REST response does not contain '
                'an assembly name for {}.'.format(species)
            )

        return assembly

    def _findEnsemblFile(
            self,
            release,
            species,
            folder,
            suffix):

        if folder == 'fasta':
            url = (
                'https://ftp.ensembl.org/pub/release-{}/'
                'fasta/{}/dna/'
            ).format(
                release,
                species
            )

        else:
            url = (
                'https://ftp.ensembl.org/pub/release-{}/'
                'gtf/{}/'
            ).format(
                release,
                species
            )

        request = urllib.request.Request(
            url,
            headers={
                'User-Agent': 'Scipion-Chem'
            }
        )

        try:
            with urllib.request.urlopen(
                    request,
                    timeout=30) as response:

                html = response.read().decode()

        except (
            urllib.error.URLError,
            TimeoutError
        ) as error:

            raise RuntimeError(
                'Cannot access Ensembl directory:\n{}\n\n{}'
                .format(
                    url,
                    error
                )
            )

        matches = re.findall(
            r'href="([^"]+)"',
            html
        )

        candidates = [
            filename
            for filename in matches
            if filename.endswith(suffix)
        ]

        if not candidates:
            raise RuntimeError(
                "No file ending with '{}' found in\n{}"
                .format(
                    suffix,
                    url
                )
            )

        candidates = sorted(candidates)

        if len(candidates) > 1:
            self.info(
                "Multiple Ensembl files found for '{}'. "
                "Using '{}'.".format(
                    suffix,
                    candidates[0]
                )
            )

        return candidates[0]

    def _buildEnsemblDownloadUrls(self, genomeInfo):

        release = self._validateEnsemblUrlComponent(
            genomeInfo['release'],
            'release'
        )

        species = self._validateEnsemblUrlComponent(
            genomeInfo['ensemblName'],
            'species'
        )

        fastaName = self._findEnsemblFile(
            release,
            species,
            'fasta',
            '.dna.toplevel.fa.gz'
        )

        fastaName = self._validateEnsemblUrlComponent(
            fastaName,
            'FASTA filename'
        )

        base = (
            'https://ftp.ensembl.org/pub/release-{}'
            .format(release)
        )

        fastaUrl = (
            '{}/fasta/{}/dna/{}'
            .format(
                base,
                species,
                fastaName
            )
        )

        gtfUrl = None

        if self.downloadAnnotation.get():
            gtfName = self._findEnsemblFile(
                release,
                species,
                'gtf',
                '.gtf.gz'
            )

            gtfName = self._validateEnsemblUrlComponent(
                gtfName,
                'GTF filename'
            )

            gtfUrl = (
                '{}/gtf/{}/{}'
                .format(
                    base,
                    species,
                    gtfName
                )
            )

        return fastaUrl, gtfUrl

    # =====================================================================
    # NCBI
    # =====================================================================

    def _resolveNcbiGenomes(self, genomes):
        """Resolve NCBI assembly accessions before downloading."""

        for genomeInfo in genomes:
            if genomeInfo.get('accession') is None:
                genomeInfo['accession'] = (
                    self._getNcbiReferenceAccession(
                        genomeInfo['ncbiTaxon']
                    )
                )

    def _getNcbiReferenceAccession(self, taxon):
        """Return the latest current NCBI RefSeq assembly accession for a taxon."""

        outputFile = self._getTmpPath(
            '{}_ncbi_summary.jsonl'.format(
                self._safeName(taxon)
            )
        )

        arguments = (
            'summary genome taxon "{}" '
            '--assembly-source RefSeq '
            '--tax-exact-match '
            '--as-json-lines '
            '> "{}"'
        ).format(
            taxon,
            outputFile
        )

        Plugin.runCondaCommand(
            self,
            arguments,
            RNASEQ_DIC,
            'datasets'
        )

        if not os.path.exists(outputFile):
            raise RuntimeError(
                'NCBI Datasets did not generate a genome summary for {}.'
                .format(taxon)
            )

        reports = []

        with open(outputFile) as inputFile:
            for line in inputFile:
                line = line.strip()

                if line:
                    reports.append(
                        json.loads(line)
                    )

        if not reports:
            raise RuntimeError(
                'NCBI does not provide a current RefSeq assembly for {}.'
                .format(taxon)
            )

        def getReleaseDate(report):
            assemblyInfo = report.get('assemblyInfo', {})

            return (
                    assemblyInfo.get('releaseDate')
                    or assemblyInfo.get('submissionDate')
                    or ''
            )

        latestReport = max(
            reports,
            key=getReleaseDate
        )

        accession = (
                latestReport.get('accession')
                or latestReport.get('currentAccession')
        )

        if not accession:
            raise RuntimeError(
                'The latest NCBI genome assembly for {} does not contain '
                'an assembly accession.'
                .format(taxon)
            )

        return accession

    # ---------------------------------------------------------------------
    # Download
    # ---------------------------------------------------------------------

    def downloadGenomesStep(self):

        genomes = self._readGenomeInfo()

        if self.source.get() == self.SOURCE_ENSEMBL:
            self._downloadEnsemblGenomes(genomes)
        else:
            self._downloadNcbiGenomes(genomes)

        self._writeGenomeInfo(genomes)

    # ---------------------------------------------------------------------
    # Ensembl download
    # ---------------------------------------------------------------------

    def _downloadEnsemblGenomes(self, genomes):

        for genomeInfo in genomes:

            outputDir = self._getExtraPath(
                '{}_{}_release-{}'.format(
                    genomeInfo['ensemblName'],
                    genomeInfo['assembly'],
                    genomeInfo['release']
                )
            )

            os.makedirs(
                outputDir,
                exist_ok=True
            )

            fastaUrl, gtfUrl = self._buildEnsemblDownloadUrls(
                genomeInfo
            )

            #FASTA
            fastaGz = os.path.join(
                outputDir,
                'genome.fa.gz'
            )

            fastaFile = os.path.join(
                outputDir,
                'genome.fa'
            )

            self._downloadAndUncompress(
                url=fastaUrl,
                compressedFile=fastaGz,
                outputFile=fastaFile
            )

            genomeInfo['fastaFile'] = fastaFile
            genomeInfo['gtfFile'] = None

            # GTF
            if self.downloadAnnotation.get():

                gtfGz = os.path.join(
                    outputDir,
                    'annotation.gtf.gz'
                )

                gtfFile = os.path.join(
                    outputDir,
                    'annotation.gtf'
                )

                self._downloadAndUncompress(
                    url=gtfUrl,
                    compressedFile=gtfGz,
                    outputFile=gtfFile
                )

                genomeInfo['gtfFile'] = gtfFile

    # ---------------------------------------------------------------------
    # NCBI download
    # ---------------------------------------------------------------------
    def _canReuseNcbiDownload(self, genomeInfo):
        existingFasta = genomeInfo.get('fastaFile')
        existingGtf = genomeInfo.get('gtfFile')

        fastaExists = (
                existingFasta
                and os.path.exists(existingFasta)
        )

        gtfExists = (
                existingGtf
                and os.path.exists(existingGtf)
        )

        return (
                not self.overwrite.get()
                and fastaExists
                and (
                        not self.downloadAnnotation.get()
                        or gtfExists
                )
        )

    def _prepareNcbiDownloadPaths(
            self,
            outputDir,
            extractDir,
            zipFile):

        if (
                self.overwrite.get()
                and os.path.exists(outputDir)
        ):
            shutil.rmtree(outputDir)

        os.makedirs(
            outputDir,
            exist_ok=True
        )

        if os.path.exists(extractDir):
            shutil.rmtree(extractDir)

        if os.path.exists(zipFile):
            os.remove(zipFile)

    @staticmethod
    def _getNcbiQueryName(genomeInfo):
        return (
                genomeInfo.get('accession')
                or genomeInfo['ncbiTaxon']
        )

    def _getNcbiGenomesToDownload(self, genomes):
        return [
            genomeInfo
            for genomeInfo in genomes
            if not self._canReuseNcbiDownload(genomeInfo)
        ]


    def _downloadNcbiGenomes(self, genomes):

        for genomeInfo in self._getNcbiGenomesToDownload(genomes):

            queryName = self._getNcbiQueryName(
                genomeInfo
            )

            outputDir = self._getExtraPath(
                'ncbi_{}_{}'.format(
                    self._safeName(
                        genomeInfo['scientificName']
                    ),
                    self._safeName(queryName)
                )
            )

            zipFile = os.path.join(
                outputDir,
                'ncbi_dataset.zip'
            )

            extractDir = os.path.join(
                outputDir,
                'package'
            )

            self._prepareNcbiDownloadPaths(
                outputDir,
                extractDir,
                zipFile
            )

            includeFiles = (
                'genome,gtf'
                if self.downloadAnnotation.get()
                else 'genome'
            )

            # -------------------------------------------------------------
            # Build NCBI Datasets command
            # -------------------------------------------------------------

            arguments = (
                'download genome accession "{}" '
                '--include {} '
                '--filename "{}" '
                '--no-progressbar'
            ).format(
                genomeInfo['accession'],
                includeFiles,
                zipFile
            )

            self._runNcbiDownload(
                arguments,
                zipFile,
                retries=3
            )

            # -------------------------------------------------------------
            # Extract NCBI package
            # -------------------------------------------------------------

            try:
                with zipfile.ZipFile(
                        zipFile,
                        'r') as archive:

                    archive.extractall(extractDir)

            except zipfile.BadZipFile as error:
                raise RuntimeError(
                    'Invalid NCBI data package: {}'
                    .format(error)
                )

            # -------------------------------------------------------------
            # Metadata
            # -------------------------------------------------------------

            reportFile = os.path.join(
                extractDir,
                'ncbi_dataset',
                'data',
                'assembly_data_report.jsonl'
            )

            report = self._readNcbiAssemblyReport(
                reportFile
            )

            accession = (
                    report.get('accession')
                    or report.get('currentAccession')
            )

            if not accession:
                raise RuntimeError(
                    'NCBI assembly report does not contain '
                    'an assembly accession.'
                )

            assemblyInfo = report.get(
                'assemblyInfo',
                {}
            )

            assemblyName = (
                    assemblyInfo.get('assemblyName')
                    or accession
            )

            organismInfo = report.get(
                'organism',
                {}
            )

            scientificName = (
                    organismInfo.get('organismName')
                    or genomeInfo['scientificName']
            )

            annotationInfo = report.get(
                'annotationInfo',
                {}
            )

            # -------------------------------------------------------------
            # Locate downloaded files
            # -------------------------------------------------------------

            fastaSource = self._findNcbiPackageFile(
                extractDir,
                ('.fna',),
                'genome FASTA'
            )

            fastaFile = os.path.join(
                outputDir,
                '{}_genomic.fna'.format(
                    accession
                )
            )

            shutil.copy2(
                fastaSource,
                fastaFile
            )

            gtfFile = None

            if self.downloadAnnotation.get():
                gtfSource = self._findNcbiPackageFile(
                    extractDir,
                    ('.gtf',),
                    'GTF annotation'
                )

                gtfFile = os.path.join(
                    outputDir,
                    '{}_genomic.gtf'.format(
                        accession
                    )
                )

                shutil.copy2(
                    gtfSource,
                    gtfFile
                )

            # -------------------------------------------------------------
            # Store resolved metadata
            # -------------------------------------------------------------

            genomeInfo['scientificName'] = scientificName
            genomeInfo['accession'] = accession
            genomeInfo['assembly'] = assemblyName

            # NCBI accession acts as the versioned genome identifier.
            genomeInfo['release'] = accession

            genomeInfo['annotationRelease'] = (
                annotationInfo.get('name')
            )

            genomeInfo['fastaFile'] = fastaFile
            genomeInfo['gtfFile'] = gtfFile

            # Package is no longer required once files and metadata have
            # been copied.
            shutil.rmtree(
                extractDir,
                ignore_errors=True
            )

            if os.path.exists(zipFile):
                os.remove(zipFile)

    # ---------------------------------------------------------------------
    # NCBI package utilities
    # ---------------------------------------------------------------------

    @staticmethod
    def _readNcbiAssemblyReport(reportFile):

        if not os.path.exists(reportFile):
            raise RuntimeError(
                'NCBI assembly report was not found: {}'
                .format(reportFile)
            )

        reports = []

        with open(reportFile) as inputFile:

            for line in inputFile:
                line = line.strip()

                if line:
                    reports.append(
                        json.loads(line)
                    )

        if not reports:
            raise RuntimeError(
                'NCBI assembly report is empty.'
            )

        if len(reports) > 1:
            raise RuntimeError(
                'The NCBI query returned more than one reference '
                'assembly. Please specify an explicit GCF_ or GCA_ '
                'assembly accession.'
            )

        return reports[0]

    @staticmethod
    def _findNcbiPackageFile(
            rootDir,
            extensions,
            description):

        candidates = []

        for root, _, files in os.walk(rootDir):

            for filename in files:

                if filename.endswith(extensions):
                    candidates.append(
                        os.path.join(
                            root,
                            filename
                        )
                    )

        if not candidates:
            raise RuntimeError(
                'NCBI package does not contain {}.'
                .format(description)
            )

        if len(candidates) > 1:
            raise RuntimeError(
                'NCBI package contains multiple files matching {}: {}'
                .format(
                    description,
                    ', '.join(candidates)
                )
            )

        return candidates[0]

    def _runNcbiDownload(
            self,
            arguments,
            zipFile,
            retries=3):
        """Run an NCBI Datasets download command with retries."""

        lastError = None

        for attempt in range(1, retries + 1):

            if os.path.exists(zipFile):
                os.remove(zipFile)

            try:
                Plugin.runCondaCommand(
                    self,
                    arguments,
                    RNASEQ_DIC,
                    'datasets'
                )

                if not os.path.exists(zipFile):
                    raise RuntimeError(
                        'NCBI Datasets finished without generating '
                        'the expected archive: {}'
                        .format(zipFile)
                    )

                return

            except Exception as error:
                lastError = error

                if os.path.exists(zipFile):
                    os.remove(zipFile)

                self.warning(
                    'NCBI download attempt {}/{} failed: {}'
                    .format(
                        attempt,
                        retries,
                        error
                    )
                )

        raise RuntimeError(
            'NCBI download failed after {} attempts: {}'
            .format(
                retries,
                lastError
            )
        )

    # =====================================================================
    # Common download utilities
    # =====================================================================
    @staticmethod
    def _validateEnsemblUrlComponent(value, name):
        value = str(value)

        if not re.fullmatch(r'[A-Za-z0-9_.-]+', value):
            raise RuntimeError(
                'Invalid Ensembl {}: {}'.format(
                    name,
                    value
                )
            )

        return value

    def _downloadAndUncompress(
            self,
            url,
            compressedFile,
            outputFile):

        if (
                os.path.exists(outputFile)
                and not self.overwrite.get()
        ):
            return

        if self.overwrite.get():

            for path in (
                    compressedFile,
                    outputFile
            ):
                if os.path.exists(path):
                    os.remove(path)

        request = urllib.request.Request(
            url,
            headers={
                'User-Agent': 'Scipion-Chem'
            }
        )

        try:
            with urllib.request.urlopen(
                    request,
                    timeout=30) as response:
                with open(
                        compressedFile,
                        'wb') as output:
                    shutil.copyfileobj(
                        response,
                        output
                    )

        except (
                urllib.error.URLError,
                TimeoutError
        ) as error:
            raise RuntimeError(
                'Could not download Ensembl file: {}'
                .format(error)
            )

        if not os.path.exists(compressedFile):
            raise RuntimeError(
                'The downloaded file was not created: {}'
                .format(compressedFile)
            )

        with gzip.open(
                compressedFile,
                'rb'
        ) as inputFile:
            with open(
                    outputFile,
                    'wb'
            ) as output:
                shutil.copyfileobj(
                    inputFile,
                    output
                )

        os.remove(compressedFile)

        if not os.path.exists(outputFile):
            raise RuntimeError(
                'The uncompressed file was not created: {}'
                .format(outputFile)
            )

    # ---------------------------------------------------------------------
    # Generic utilities
    # ---------------------------------------------------------------------

    @staticmethod
    def _removeGzipExtension(path):
        return (
            path[:-3]
            if path.endswith('.gz')
            else path
        )

    @staticmethod
    def _safeName(value):
        return re.sub(
            r'[^A-Za-z0-9_.-]+',
            '_',
            value
        )

    @staticmethod
    def _speciesIdentifierToScientificName(species):

        species = species.replace(
            '_',
            ' '
        )

        return ' '.join(
            word.capitalize()
            for word in species.split()
        )

    @staticmethod
    def _requestJson(url, description):

        request = urllib.request.Request(
            url,
            headers={
                'Accept': 'application/json',
                'User-Agent': 'Scipion-Chem'
            }
        )

        try:
            with urllib.request.urlopen(
                    request,
                    timeout=30) as response:

                return json.loads(
                    response.read().decode('utf-8')
                )

        except (
            urllib.error.URLError,
            TimeoutError,
            json.JSONDecodeError
        ) as error:

            raise RuntimeError(
                'Could not retrieve {}: {}'
                .format(
                    description,
                    error
                )
            )

    # =====================================================================
    # Persistent genome information
    # =====================================================================

    def _getGenomeInfoFile(self):
        return self._getExtraPath(
            'genomes.json'
        )

    def _writeGenomeInfo(self, genomes):

        with open(
                self._getGenomeInfoFile(),
                'w') as outputFile:

            json.dump(
                genomes,
                outputFile,
                indent=2,
                sort_keys=True
            )

    def _readGenomeInfo(self):

        with open(
                self._getGenomeInfoFile()
        ) as inputFile:

            return json.load(inputFile)

    # =====================================================================
    # Output
    # =====================================================================

    def createOutputStep(self):

        genomesInfo = self._readGenomeInfo()

        referenceGenomes = SetOfGenomes().create(
            outputPath=self._getPath()
        )

        referenceGenomes.setObjLabel(
            'Reference genomes'
        )

        for genomeInfo in genomesInfo:

            genome = Genome()

            genome.setScientificName(
                genomeInfo['scientificName']
            )

            if genomeInfo['source'] == 'Ensembl':
                genome.setEnsemblSpecies(
                    genomeInfo['ensemblName']
                )

            genome.setAssembly(
                genomeInfo['assembly']
            )

            genome.setSource(
                genomeInfo['source']
            )

            genome.setRelease(
                genomeInfo['release']
            )

            genome.setFastaFile(
                genomeInfo['fastaFile']
            )

            if genomeInfo.get('gtfFile'):
                genome.setGtfFile(
                    genomeInfo['gtfFile']
                )

            genome.setObjLabel(
                '{} ({})'.format(
                    genomeInfo['scientificName'],
                    genomeInfo['assembly']
                )
            )

            referenceGenomes.append(
                genome
            )

        self._defineOutputs(
            referenceGenomes=referenceGenomes
        )

    # =====================================================================
    # Validation
    # =====================================================================

    def _validate(self):

        errors = []

        try:
            genomes = self._getSelectedSpecies()

        except ValueError as error:
            errors.append(str(error))
            genomes = []

        if self.source.get() == self.SOURCE_ENSEMBL:

            self._validateEnsemblParameters(
                genomes,
                errors
            )

        else:

            self._validateNcbiParameters(
                genomes,
                errors
            )

        return errors

    # ---------------------------------------------------------------------
    # Ensembl validation
    # ---------------------------------------------------------------------

    def _validateEnsemblParameters(
            self,
            genomes,
            errors):

        assemblies = self._splitParameterValues(
            self.assemblies.get()
        )

        if not assemblies:

            errors.append(
                'Assemblies cannot be empty.'
            )

        elif (
            genomes
            and len(assemblies) not in (
                1,
                len(genomes)
            )
        ):
            errors.append(
                'Assemblies must contain one value or one value '
                'per genome.'
            )

        releases = self._splitParameterValues(
            self.releases.get()
        )

        if not releases:

            errors.append(
                'Ensembl releases cannot be empty.'
            )

        elif (
            genomes
            and len(releases) not in (
                1,
                len(genomes)
            )
        ):
            errors.append(
                'Releases must contain one value or one value '
                'per genome.'
            )

        for release in releases:

            self._validateReleaseValue(
                release,
                'Ensembl release',
                errors
            )

        for species in genomes:

            if ' ' in species:

                errors.append(
                    'Ensembl species identifiers must use '
                    'underscores: {}'.format(species)
                )

    # ---------------------------------------------------------------------
    # NCBI validation
    # ---------------------------------------------------------------------

    def _validateNcbiParameters(
            self,
            genomes,
            errors):

        assemblies = self._splitParameterValues(
            self.ncbiAssemblies.get()
        )

        if not assemblies:

            errors.append(
                'NCBI assemblies cannot be empty.'
            )

            return

        if (
            genomes
            and len(assemblies) not in (
                1,
                len(genomes)
            )
        ):
            errors.append(
                'NCBI assemblies must contain one value or '
                'one value per genome.'
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
                    'NCBI assembly must be "Latest" or a valid '
                    'versioned GCF_/GCA_ accession: {}'
                    .format(assembly)
                )

    @staticmethod
    def _validateReleaseValue(
            release,
            name,
            errors):

        if release.lower() == 'latest':
            return

        try:
            value = int(release)

            if value <= 0:
                raise ValueError

        except ValueError:

            errors.append(
                '{} must be "Latest" or a positive integer.'
                .format(name)
            )

    # =====================================================================
    # Summary
    # =====================================================================

    def _summary(self):

        if not hasattr(
                self,
                'referenceGenomes'):

            return [
                'No reference genomes have been downloaded.'
            ]

        summary = [
            'Downloaded reference genomes: {}'
            .format(
                self.referenceGenomes.getSize()
            )
        ]

        for genome in self.referenceGenomes:

            if genome.getSource() == 'Ensembl':

                text = (
                    '{} ({}, Ensembl release {})'
                    .format(
                        genome.getScientificName(),
                        genome.getAssembly(),
                        genome.getRelease()
                    )
                )

            else:

                text = (
                    '{} ({}, NCBI accession {})'
                    .format(
                        genome.getScientificName(),
                        genome.getAssembly(),
                        genome.getRelease()
                    )
                )

            if genome.hasGtfFile():
                text += ' + GTF'

            summary.append(text)

        return summary

    # =====================================================================
    # Methods
    # =====================================================================

    def _methods(self):

        if not hasattr(
                self,
                'referenceGenomes'):

            return []

        source = (
            'Ensembl'
            if self.source.get() == self.SOURCE_ENSEMBL
            else 'NCBI'
        )

        methods = [
            '{} reference genome(s) were downloaded from {}.'
            .format(
                self.referenceGenomes.getSize(),
                source
            )
        ]

        for genome in self.referenceGenomes:

            methods.append(
                '{} ({}, {})'.format(
                    genome.getScientificName(),
                    genome.getAssembly(),
                    genome.getRelease()
                )
            )

        return methods