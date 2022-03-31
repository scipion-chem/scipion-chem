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


class ProtChemMultipleSequenceAlignment(EMProtocol):
    """Run multiple sequence alignment for a set of sequences"""
    _label = 'multiple sequence alignment'

    def _defineParams(self, form):

        form.addSection(label='Input')

        form.addParam('setOfSequences', PointerParam, pointerClass='SetOfSequences',
                      label='Input Set of Sequence File: ', allowsNull=False,
                      help="Set of sequences to be alignment ")

        form.addParam('programList', EnumParam, choices=['Clustal Omega', 'Muscle'], default=0,
                      label='Select a multiple sequence alignment program: ',
                      help="Program selected to run the alignment")

        form.addParam('alignmentOutputFileClustal', EnumParam, 
                      choices=['FASTA', 'ClustalO', 'MSF', 'PHYLIP', 'SELEX', 'STOCKHOLM', 'VIENNA'],
                      condition='programList == 0', expertLevel=LEVEL_ADVANCED,
                      label='Alignment output format: ',
                      help="Writes the output in Clustal Omega format selected")

        form.addParam('alignmentOutputFileMuscle', EnumParam, 
                      choices=['FASTA', 'ClustalW', 'Clustal(strick)', 'HTML', 'MSF'],
                      condition='programList == 1', expertLevel=LEVEL_ADVANCED,
                      label='Alignment Output Format: ',
                      help="Writes the output in Muscle format selected")

        form.addParam('additionalFormat', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      label='Additional sequence format: ', default=True, 
                      help="Other standard formats from EMBOSS seqret")

        form.addParam('embossFormats', EnumParam, 
                      choices=['Wisconsin Package GCG 9.x and 10.x', 'GCG 8.x', 'SwisProt', 'NCBI', 'NBRF (PIR)', 
                               'Intelligenetics', 'CODATA', 'DNA strider', 'ACeDB', 
                               '"gap" program in the Staden package', 'Plain sequence', 'Fitch', 'PHYLIP interleaved', 
                               'ASN.1', 'Hennig86', 'Mega', 'Meganon', 'Nexus/PAUP', 'Nexusnon/PAUPnon', 'Jackknifer', 
                               'Jackknifernon', 'Treecon', 'EMBOSS sequence object report'],
                      condition='additionalFormat', expertLevel=LEVEL_ADVANCED,
                      label='EMBOSS seq output format: ',
                      help="Sequence formats from EMBOSS seqret")


    def _insertAllSteps(self):
        self._insertFunctionStep('multipleAlignment')

    def multipleAlignment(self, outputPath=None):
        programName = self.getEnumText('programList')
        print('Program Name: ', programName)
        setForAlignment = self.setOfSequences.get()
        input_file= self._getPath('SequencesForAlignment.fasta')

        # Clustal Omega
        if programName == 'Clustal Omega':
            formatSelectedClustal = self.alignmentOutputFileClustal.get()
            clustaloFormats = ['FASTA', 'ClustalO', 'MSF', 'PHYLIP', 'SELEX', 'STOCKHOLM', 'VIENNA']
            clustalFormatsArgs = ['fasta', 'clu', 'msf', 'phy', 'selex', 'stk', 'vie']
            outputFormat = clustaloFormats[formatSelectedClustal]
            print('Output Format: ', outputFormat)
            rootProgram = "clustalOmega."
            if outputFormat == 'ClustalO':
                output_file = self._getPath(rootProgram + 'aln')
            else:
                output_file = self._getPath(rootProgram + clustalFormatsArgs[formatSelectedClustal])

            # Alignment
            run_alignment = setForAlignment.clustalOmegaAlignSequences(input_file, output_file)
            cline = run_alignment
            args = ' --outfmt=' + clustalFormatsArgs[formatSelectedClustal]
            self.runJob(cline, args)

        # Muscle
        if programName == 'Muscle':
            formatSelectedMuscle = self.alignmentOutputFileMuscle.get()
            muscleFormats = ['FASTA', 'ClustalW', 'Clustal(strick)', 'HTML', 'MSF']
            muscleFormatsArgs = ['fasta', 'clw', 'clwstrict', 'html', 'msf']
            outputFormat = muscleFormats[formatSelectedMuscle]
            print('Output Format: ', outputFormat)
            rootProgram = "muscle."
            if outputFormat == 'ClustalW' or outputFormat == 'Clustal(strick)':
                output_file = self._getPath(rootProgram + 'aln')
            else:
                output_file = self._getPath(rootProgram + muscleFormatsArgs[formatSelectedMuscle])

            # Alignment
            run_alignment = setForAlignment.muscleAlignmentSequences(input_file, output_file)
            cline = run_alignment
            args = ' -' + muscleFormatsArgs[formatSelectedMuscle]
            self.runJob(cline, args)

        outputName = os.path.basename(output_file)
        out_fileAligned = SequencesAlignment(alignmentFileName=outputName)
        self._defineOutputs(outputVariants=out_fileAligned)

        # EMBOSS format
        if self.additionalFormat:
            embossAddFormats_index = self.embossFormats.get()
            embossAddFormatsArgs = ['gcg', 'gcg8', 'sw', 'ncbi', 'pir', 'ig', 'codata', 'strider', 'acedb', 
                                    'experiment', 'plain', 'fitch', 'phylip3', 'asn1', 'hennig86', 'mega', 
                                    'meganon', 'nexus', 'nexusnon', 'jackknifer', 'jackknifernon', 'treecon', 'debug']
            embossSelectedFormat = embossAddFormatsArgs[embossAddFormats_index]
            StandardFormatEmboss = self._getPath(rootProgram + embossSelectedFormat)
            emboss = out_fileAligned.convertEMBOSSformat(output_file, embossSelectedFormat,StandardFormatEmboss)
            arguments = ''
            self.runJob(emboss, arguments)


    def _validate(self):
        errors = []
        return errors
