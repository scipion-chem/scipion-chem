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
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************
import json
import os
import re


COMMON_SPECIES = {
    "Homo sapiens": {
        "ensembl": "homo_sapiens",
        "ncbi": "Homo sapiens",
        "assembly": "GRCh38",
    },
    "Mus musculus": {
        "ensembl": "mus_musculus",
        "ncbi": "Mus musculus",
        "assembly": "GRCm39",
    },
    "Rattus norvegicus": {
        "ensembl": "rattus_norvegicus",
        "ncbi": "Rattus norvegicus",
        "assembly": "mRatBN7.2",
    },
    "Danio rerio": {
        "ensembl": "danio_rerio",
        "ncbi": "Danio rerio",
        "assembly": "GRCz11",
    },
    "Drosophila melanogaster": {
        "ensembl": "drosophila_melanogaster",
        "ncbi": "Drosophila melanogaster",
        "assembly": "BDGP6.46",
    },
    "Caenorhabditis elegans": {
        "ensembl": "caenorhabditis_elegans",
        "ncbi": "Caenorhabditis elegans",
        "assembly": "WBcel235",
    },
    "Saccharomyces cerevisiae": {
        "ensembl": "saccharomyces_cerevisiae",
        "ncbi": "Saccharomyces cerevisiae",
        "assembly": "R64-1-1",
    },
    "Arabidopsis thaliana": {
        "ensembl": "arabidopsis_thaliana",
        "ncbi": "Arabidopsis thaliana",
        "assembly": "TAIR10",
    },
}

def getCommonSpecies():
    """Return the scientific names of the predefined common species.

    Returns
    -------
    list
        Scientific species names available as common genomic resources.
    """
    return list(COMMON_SPECIES.keys())

def getProviderSpeciesName(species, provider):
    """Return the species identifier expected by a genomic data provider.

    Parameters
    ----------
    species : str
        Scientific species name, for example ``"Homo sapiens"``.
    provider : str
        Genomic data provider. Currently supported values are
        ``"ensembl"`` and ``"ncbi"``.

    Returns
    -------
    str
        Species identifier expected by the selected provider.

    Notes
    -----
    If the species is not included in ``COMMON_SPECIES``, a generic
    conversion is used for Ensembl and the original scientific name is
    preserved for NCBI. This allows custom species to be used without
    requiring them to be predefined. Whitespace is normalized and common
    species names are recognized regardless of capitalization. Empty names
    and unsupported providers raise ValueError.
    """
    provider = (provider or '').strip().lower()
    species = ' '.join((species or '').split())
    if not species:
        raise ValueError('A scientific species name is required.')

    if provider not in ("ensembl", "ncbi"):
        raise ValueError(
            f"Unsupported genomic data provider: {provider}"
        )

    commonName = next((name for name in COMMON_SPECIES if name.casefold() == species.casefold()), None)
    speciesInfo = COMMON_SPECIES.get(commonName)

    if speciesInfo is not None:
        return speciesInfo[provider]

    if provider == "ensembl":
        return species.lower().replace(" ", "_")

    return species


def parseCustomSpecies(speciesText):
    """Parse a comma- or semicolon-separated list of custom species.

    Parameters
    ----------
    speciesText : str
        Species entered by the user, separated by commas or semicolons.
        For example::

            Homo sapiens;Mus musculus;Canis lupus familiaris

    Returns
    -------
    list
        Cleaned species names. Empty entries are ignored and duplicate
        species are removed case-insensitively while preserving their order.
        Repeated whitespace is collapsed and common species use their
        canonical scientific-name capitalization.
    """
    if not speciesText:
        return []

    commonNames = {name.casefold(): name for name in COMMON_SPECIES}
    selected = {}
    for item in re.split(r'[;,]', speciesText):
        name = ' '.join(item.split())
        if name:
            key = name.casefold()
            selected.setdefault(key, commonNames.get(key, name))
    return list(selected.values())


class GenomicResourceMixin:
    """Shared parameter and metadata helpers for genome/variant downloads.

    Protocols supply their Ensembl REST URL and their own _requestJson method
    so provider-specific network validation remains in place.
    """

    @staticmethod
    def splitParameterValues(value):
        return [item.strip() for item in re.split('[;,]', value or '') if item.strip()]

    @classmethod
    def expandParameterValues(cls, value, count, parameterName):
        values = cls.splitParameterValues(value)
        if not values:
            raise ValueError('{} cannot be empty.'.format(parameterName.capitalize()))
        if len(values) == 1:
            return values * count
        if len(values) != count:
            raise ValueError('The number of {} must be either 1 or equal to the number of species.'
                             .format(parameterName))
        return values

    def getLatestEnsemblRelease(self):
        data = self._requestJson(self.ENSEMBL_REST_DATA_URL, 'latest Ensembl release')
        releases = data.get('releases', [])
        if not releases:
            raise RuntimeError('The Ensembl REST response does not contain releases.')
        return max(int(release) for release in releases)

    @staticmethod
    def readJsonLines(filename):
        if not os.path.exists(filename):
            return []
        with open(filename) as handle:
            return [json.loads(line) for line in handle if line.strip()]

    @staticmethod
    def isNonemptyFile(path):
        return bool(path) and os.path.isfile(path) and os.path.getsize(path) > 0

    @staticmethod
    def safeName(value):
        return re.sub('[^A-Za-z0-9_.-]+', '_', str(value))
