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
from pwchem.objects import Genome, SetOfGenomes


class ProtImportGenomes(EMProtocol):
    """Import one or more reference genomes."""

    _label = 'import genomes'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'genomesData',
            params.StringParam,
            default='[]',
            label='Genomes',
            help='Configure the genomes to import.'
        )

    def _insertAllSteps(self):
        self._insertFunctionStep(
            self.createOutputStep
        )

    def createOutputStep(self):
        genomesData = json.loads(self.genomesData.get())

        referenceGenomes = SetOfGenomes().create(
            outputPath=self._getPath()
        )
        referenceGenomes.setObjLabel('Reference genomes')

        for genomeData in genomesData:
            genome = Genome()

            scientificName = genomeData.get('scientificName', '')
            assembly = genomeData.get('assembly', '')
            release = genomeData.get('release', '')
            fastaFile = genomeData.get('fastaFile', '')
            gtfFile = genomeData.get('gtfFile', '')

            if scientificName:
                genome.setScientificName(scientificName)

            if assembly:
                genome.setAssembly(assembly)

            if release:
                genome.setRelease(release)

            genome.setSource('Imported')

            # Use the scientific name to create readable and unique filenames.
            safeName = (
                scientificName.strip().replace(' ', '_')
                if scientificName
                else 'genome'
            )

            if assembly:
                safeName += '_{}'.format(
                    assembly.strip().replace(' ', '_')
                )

            if fastaFile:
                fastaExt = os.path.splitext(fastaFile)[1]

                importedFasta = self._getExtraPath(
                    '{}_genome{}'.format(
                        safeName,
                        fastaExt
                    )
                )

                copyFile(
                    fastaFile,
                    importedFasta
                )

                # Create the FASTA index required by genome browsers
                # and downstream tools.
                self._createFastaIndex(
                    importedFasta
                )

                genome.setFastaFile(
                    importedFasta
                )

            if gtfFile:
                importedGtf = self._getExtraPath(
                    '{}_annotation.gtf'.format(
                        safeName
                    )
                )

                copyFile(
                    gtfFile,
                    importedGtf
                )

                genome.setGtfFile(
                    importedGtf
                )

            label = scientificName or 'Genome'

            if assembly:
                label += ' ({})'.format(
                    assembly
                )

            genome.setObjLabel(
                label
            )

            referenceGenomes.append(
                genome
            )

        self._defineOutputs(
            referenceGenomes=referenceGenomes
        )

    def _createFastaIndex(self, fastaFile):
        """Create the samtools FASTA index (.fai)."""

        fastaIndex = fastaFile + '.fai'

        if os.path.exists(fastaIndex):
            return

        args = 'faidx {}'.format(
            self._quote(fastaFile)
        )

        Plugin.runCondaCommand(
            self,
            args,
            RNASEQ_DIC,
            'samtools'
        )

        if not os.path.isfile(fastaIndex):
            raise RuntimeError(
                'FASTA index was not created: {}'.format(
                    fastaIndex
                )
            )

    @staticmethod
    def _quote(value):
        """Quote a filesystem path for shell execution."""
        return shlex.quote(
            str(value)
        )

    def _validate(self):
        errors = []

        try:
            genomesData = json.loads(
                self.genomesData.get()
            )
        except (TypeError, ValueError):
            return [
                'Invalid genomes configuration.'
            ]

        if not genomesData:
            errors.append(
                'At least one genome must be configured.'
            )
            return errors

        for index, genomeData in enumerate(
            genomesData,
            start=1
        ):
            fastaFile = genomeData.get(
                'fastaFile',
                ''
            )

            if not fastaFile:
                errors.append(
                    'Genome {}: FASTA file is required.'.format(
                        index
                    )
                )

            elif not os.path.isfile(fastaFile):
                errors.append(
                    'Genome {}: FASTA file does not exist: {}'.format(
                        index,
                        fastaFile
                    )
                )

            gtfFile = genomeData.get(
                'gtfFile',
                ''
            )

            if (
                gtfFile
                and not os.path.isfile(gtfFile)
            ):
                errors.append(
                    'Genome {}: GTF file does not exist: {}'.format(
                        index,
                        gtfFile
                    )
                )

        return errors