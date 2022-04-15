# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import EnumParam, PointerParam, BooleanParam, LEVEL_ADVANCED
from pwchem.objects import SequencesAlignment
from pwchem.utils.utilsFasta import clustalOmegaAlignSequences, muscleAlignmentSequences

CLUSTALO, MUSCLE, MAFFT = 'Clustal_Omega', 'Muscle', 'Mafft'

EMBOSS_FORMATS = {'Fasta': 'fasta', 'Clustal': 'aln', 'Wisconsin Package GCG 9.x and 10.x': 'gcg', 'GCG 8.x': 'gcg8',
                  'SwisProt': 'sw', 'NCBI': 'ncbi',
                  'NBRF (PIR)': 'pir', 'Intelligenetics': 'ig', 'CODATA': 'codata', 'DNA strider': 'strider',
                  'ACeDB': 'acedb', '"gap" program in the Staden package': 'experiment', 'Plain sequence': 'plain',
                  'Fitch': 'fitch', 'PHYLIP interleaved': 'phylip3', 'ASN.1': 'asn1', 'Hennig86': 'hennig86',
                  'Mega': 'mega', 'Meganon': 'meganon', 'Nexus/PAUP': 'nexus', 'Nexusnon/PAUPnon': 'nexusnon',
                  'Jackknifer': 'jackknifer','Jackknifernon': 'jackknifernon', 'Treecon': 'treecon',
                  'EMBOSS sequence object report': 'debug'}


class ProtChemMultipleSequenceAlignment(EMProtocol):
    """Run multiple sequence alignment for a set of sequences"""
    _label = 'multiple sequence alignment'

    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('setOfSequences', PointerParam, pointerClass='SetOfSequences',
                      label='Input Set of Sequence File: ', allowsNull=False,
                      help="Set of sequences to be alignment ")

        form.addParam('programList', EnumParam, choices=[CLUSTALO, MUSCLE, MAFFT], default=0,
                      label='Select a multiple sequence alignment program: ',
                      help="Program selected to run the alignment")
        form.addParam('muscleAlg', EnumParam, choices=['Slow', 'Fast'], default=0,
                      label='Muscle algorithm to use: ', condition='programList==1',
                      help="Muscle algorithm to use:\n -Slow: PPP\n -Fast: Super5")

        form.addParam('additionalFormat', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      label='Additional sequence format: ', default=False,
                      help="Other standard formats from EMBOSS seqret")

        form.addParam('embossFormats', EnumParam, default=1,
                      choices=list(EMBOSS_FORMATS.keys()),
                      condition='additionalFormat', expertLevel=LEVEL_ADVANCED,
                      label='EMBOSS seq output format: ',
                      help="Sequence formats from EMBOSS seqret")

    def _insertAllSteps(self):
        self._insertFunctionStep('multipleAlignment')

    def multipleAlignment(self):
        programName = self.getEnumText('programList')
        print('Aligning with: ', programName)
        setForAlignment = self.setOfSequences.get()
        input_file = self._getPath('SequencesForAlignment.fasta')

        # Clustal Omega
        if programName == CLUSTALO:
            output_file = self._getPath('clustal_Omega.aln')
            args = ' --outfmt=clu'
            cline = clustalOmegaAlignSequences(setForAlignment, input_file, output_file)

        # Muscle
        elif programName == MUSCLE:
            output_file = self._getPath('muscle.fa')

            setForAlignment.exportToFile(input_file)
            alg = 'align' if self.muscleAlg.get() == 0 else 'super5'
            cline = 'muscle -{} {} -output {}'.format(alg, input_file, output_file)

        elif programName == MAFFT:
            output_file = self._getPath('mafft.aln')
            setForAlignment.exportToFile(input_file)
            cline = 'mafft --auto --clustalout {} > {}'.format(input_file, output_file)

        self.runJob(cline, '')

        out_fileAligned = SequencesAlignment(alignmentFileName=os.path.abspath(output_file))
        self._defineOutputs(outputAlignment=out_fileAligned)

        # EMBOSS format
        if self.additionalFormat:
            embossSelectedFormat = EMBOSS_FORMATS[self.getEnumText('embossFormats')]
            StandardFormatEmboss = self._getPath('{}.{}'.format(programName.lower(), embossSelectedFormat))
            emboss = out_fileAligned.convertEMBOSSformat(output_file, embossSelectedFormat, StandardFormatEmboss)
            self.runJob(emboss, '')

    def _validate(self):
        errors = []
        return errors
