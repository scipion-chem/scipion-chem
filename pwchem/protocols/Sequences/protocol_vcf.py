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
from pwchem.utils.sequence_utils import (
    COMMON_SPECIES,
    getProviderSpeciesName,
    parseCustomSpecies
)


class ProtVCF(EMProtocol):
    """Download known-variant VCF files from Ensembl or NCBI dbSNP."""

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

    # ---------------------------------------------------------------------
    # NCBI dbSNP
    # ---------------------------------------------------------------------

    NCBI_DBSNP_VCF_URL = (
        'https://ftp.ncbi.nih.gov/snp/latest_release/VCF/'
    )

    # ---------------------------------------------------------------------
    # Allowed remote hosts
    # ---------------------------------------------------------------------

    ALLOWED_HOSTS = {
        'rest.ensembl.org',
        'ftp.ensembl.org',
        'ftp.ncbi.nih.gov'
    }
    EMPTY_PARAMETER_ERROR = '{} cannot be empty.'
    ENSEMBL_SPECIES = 'Ensembl species'

    _label = 'download VCF'

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
                'Select the database used to download known '
                'variant VCF files.'
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
                'Use "Latest" to resolve the current RefSeq assembly '
                'for every selected species.\n\n'
                'Alternatively, provide a versioned NCBI assembly '
                'accession (GCF_ or GCA_).\n\n'
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
                'indexFile': None
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
                'source': 'NCBI',
                'database': 'dbSNP',
                'variantType': 'known',
                'vcfFile': None,
                'indexFile': None
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

        for info in vcfInfo:

            if info.get('accession') is None:
                info['accession'] = (
                    self._getNcbiReferenceAccession(
                        info['ncbiTaxon']
                    )
                )

            assemblyName = (
                self._getNcbiAssemblyName(
                    info['accession']
                )
            )

            info['assembly'] = assemblyName

    def _getNcbiReferenceAccession(
        self,
        taxon
    ):

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

        reports = self._readJsonLines(
            outputFile
        )

        if not reports:
            raise RuntimeError(
                'NCBI does not provide a current RefSeq '
                'assembly for {}.'.format(
                    taxon
                )
            )

        def getReleaseDate(report):

            assemblyInfo = report.get(
                'assemblyInfo',
                {}
            )

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
                'The NCBI assembly report for {} does not '
                'contain an accession.'.format(
                    taxon
                )
            )

        return accession

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

        assemblyInfo = report.get(
            'assemblyInfo',
            {}
        )

        return (
            assemblyInfo.get('assemblyName')
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

            release = self._validateUrlComponent(
                info['release'],
                'Ensembl release'
            )

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
                    release
                )
            )

            os.makedirs(
                outputDir,
                exist_ok=True
            )

            filename = self._findEnsemblVCF(
                info
            )

            filename = self._validateUrlComponent(
                filename,
                'Ensembl VCF filename'
            )

            baseUrl = (
                '{}release-{}/variation/vcf/{}/'
            ).format(
                self.ENSEMBL_FTP_BASE_URL,
                release,
                ensemblName
            )

            vcfUrl = (
                baseUrl
                + filename
            )

            outputVCF = os.path.join(
                outputDir,
                'variants.vcf.gz'
            )

            self._downloadFile(
                vcfUrl,
                outputVCF
            )

            indexUrl = (
                vcfUrl
                + '.tbi'
            )

            indexFile = (
                outputVCF
                + '.tbi'
            )

            if self._remoteFileExists(
                indexUrl
            ):
                self._downloadFile(
                    indexUrl,
                    indexFile
                )
            else:
                indexFile = None

            info['vcfFile'] = outputVCF
            info['indexFile'] = indexFile

    def _findEnsemblVCF(
        self,
        info
    ):

        release = self._validateUrlComponent(
            info['release'],
            'Ensembl release'
        )

        ensemblName = self._validateUrlComponent(
            info['ensemblName'],
            self.ENSEMBL_SPECIES
        )

        url = (
            '{}release-{}/variation/vcf/{}/'
        ).format(
            self.ENSEMBL_FTP_BASE_URL,
            release,
            ensemblName
        )

        html = self._readRemoteDirectory(
            url
        )

        matches = re.findall(
            r'href="([^"]+)"',
            html
        )

        candidates = [
            filename
            for filename in matches
            if filename.endswith('.vcf.gz')
        ]

        if not candidates:
            raise RuntimeError(
                'Ensembl does not provide a VCF file for '
                '{} in release {}.'.format(
                    info['scientificName'],
                    release
                )
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
                'Multiple Ensembl VCF files found for {}. '
                'Using {}.'.format(
                    info['scientificName'],
                    candidates[0]
                )
            )

        return self._validateUrlComponent(
            candidates[0],
            'Ensembl VCF filename'
        )

    # =====================================================================
    # NCBI dbSNP download
    # =====================================================================

    def _downloadNcbiVCFs(
        self,
        vcfInfo
    ):

        directoryHtml = self._readRemoteDirectory(
            self.NCBI_DBSNP_VCF_URL
        )

        availableFiles = set(
            re.findall(
                r'href="([^"]+)"',
                directoryHtml
            )
        )

        for info in vcfInfo:

            accession = self._validateUrlComponent(
                info['accession'],
                'NCBI accession'
            )

            filename = '{}.gz'.format(
                accession
            )

            indexName = (
                filename
                + '.tbi'
            )

            if filename not in availableFiles:
                raise RuntimeError(
                    'NCBI dbSNP does not provide a VCF for '
                    '{} ({}, {}).'.format(
                        info['scientificName'],
                        info['assembly'],
                        accession
                    )
                )

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
                'variants.vcf.gz'
            )

            self._downloadFile(
                self.NCBI_DBSNP_VCF_URL
                + filename,
                outputVCF
            )

            indexFile = None

            if indexName in availableFiles:

                indexFile = (
                    outputVCF
                    + '.tbi'
                )

                self._downloadFile(
                    self.NCBI_DBSNP_VCF_URL
                    + indexName,
                    indexFile
                )

            info['vcfFile'] = outputVCF
            info['indexFile'] = indexFile


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
            ) as response:

                with open(
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
                '{} ({}, {}, {})'.format(
                    vcf.getScientificName(),
                    vcf.getAssembly(),
                    vcf.getSource(),
                    vcf.getDatabase()
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

        source = (
            'Ensembl'
            if self.source.get() == self.SOURCE_ENSEMBL
            else 'NCBI dbSNP'
        )

        return [
            '{} known-variant VCF file(s) were downloaded '
            'from {}.'.format(
                self.outputVCFs.getSize(),
                source
            )
        ]