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

import gzip, os, re
from pwchem import Plugin
from pwchem.constants import RNASEQ_DIC

def openFastq(fn):
  """Open plain or gzip-compressed FASTQ files."""

  if fn.endswith('.gz'):
    return gzip.open(fn, 'rt')

  return open(fn, 'r')


def getFastqStats(fn):
  """Calculate number of reads and mean read length from a FASTQ file."""

  numReads = 0
  totalLength = 0

  with openFastq(fn) as f:
    while True:
      header = f.readline()

      if not header:
        break

      seq = f.readline().strip()
      f.readline()
      f.readline()

      numReads += 1
      totalLength += len(seq)

  readLength = int(round(totalLength / numReads)) if numReads else 0

  return numReads, readLength

def runFastqc(protocol, fastqFiles):
    """Run FastQC and return the generated HTML reports.

    Parameters
    ----------
    protocol
        Scipion protocol executing FastQC.
    fastqFiles : list
        FASTQ file paths to analyze.

    Returns
    -------
    list
        FastQC HTML report paths in the same order as fastqFiles.
    """
    outDir = protocol._getExtraPath('fastqc')
    os.makedirs(outDir, exist_ok=True)

    arguments = f'-o "{outDir}"'
    arguments += ''.join(f' "{fn}"' for fn in fastqFiles)

    Plugin.runCondaCommand(
        protocol,
        arguments,
        RNASEQ_DIC,
        'fastqc'
    )

    htmlFiles = []

    for fn in fastqFiles:
        baseName = os.path.basename(fn)

        for ext in ['.fastq.gz', '.fq.gz', '.fastq', '.fq']:
            if baseName.endswith(ext):
                baseName = baseName[:-len(ext)]
                break

        htmlFile = os.path.join(
            outDir,
            f'{baseName}_fastqc.html'
        )

        if not os.path.exists(htmlFile):
            raise RuntimeError(
                f'FastQC report was not generated for: {fn}'
            )

        htmlFiles.append(htmlFile)

    return htmlFiles

# ============================================================================
# Species utilities
# ============================================================================
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
    requiring them to be predefined.
    """
    provider = provider.lower()
    species = species.strip()

    if provider not in ("ensembl", "ncbi"):
        raise ValueError(
            f"Unsupported genomic data provider: {provider}"
        )

    speciesInfo = COMMON_SPECIES.get(species)

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
        species are removed while preserving their original order.
    """
    if not speciesText:
        return []

    species = [
        item.strip()
        for item in re.split(r'[;,]', speciesText)
        if item.strip()
    ]

    return list(dict.fromkeys(species))