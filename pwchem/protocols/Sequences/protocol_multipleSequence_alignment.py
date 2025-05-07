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

from pyworkflow.protocol.params import EnumParam, PointerParam, BooleanParam, LEVEL_ADVANCED, TextParam
from pwem.protocols import EMProtocol

from pwchem.utils.utilsFasta import EMBOSS_FORMATS, parseFasta, parseAlnFile, getMultipleAlignmentCline
from pwchem.objects import SetOfSequencesChem, Sequence

CLUSTALO, MUSCLE, MAFFT = 'Clustal_Omega', 'Muscle', 'Mafft'

class ProtChemMultipleSequenceAlignment(EMProtocol):
    """Runs multiple sequence alignment for a set of sequences"""
    _label = 'multiple sequence alignment'

    def _defineParams(self, form):

        form.addSection(label='Input')
        group = form.addGroup('Input')
        group.addParam('inputSequences', PointerParam, pointerClass='SetOfSequences',
                      label='Input Set of Sequence File: ', allowsNull=False,
                      help="Set of sequences to be alignment ")

        group = form.addGroup('Program')
        group.addParam('programList', EnumParam, choices=[CLUSTALO, MUSCLE, MAFFT], default=0,
                      label='Select a multiple sequence alignment program: ',
                      help="Program selected to run the alignment.\n")

        group.addParam('additionalFormat', BooleanParam, expertLevel=LEVEL_ADVANCED,
                      label='Create additional output in other EMBOSS format: ', default=False,
                      help="Other standard formats from EMBOSS seqret")

        group.addParam('embossFormats', EnumParam, default=1,
                      choices=list(EMBOSS_FORMATS.keys()),
                      condition='additionalFormat', expertLevel=LEVEL_ADVANCED,
                      label='EMBOSS output format: ',
                      help="Sequence formats from EMBOSS seqret")

        group = form.addGroup('Parameters')
        group.addParam('extraParams', TextParam, default="",
                       label='Extra parameters (Auto params by default): ',
                       help="Extra parameters to use in the alignment:\n"
                            "CLUSTALO: http://www.clustal.org/omega/README\n"
                            "MUSCLE: https://drive5.com/muscle5/manual/commands.html\n"
                            "MAFFT: https://mafft.cbrc.jp/alignment/software/manual/manual.html")



    def _insertAllSteps(self):
        self._insertFunctionStep('multipleAlignment')

    def multipleAlignment(self):
        programName = self.getEnumText('programList')
        print('Aligning with: ', programName)
        setForAlignment = self.inputSequences.get()
        inFile = os.path.abspath(self._getPath('SequencesForAlignment.fasta'))
        setForAlignment.exportToFile(inFile)

        # Clustal Omega
        if programName == CLUSTALO:
            outFile = os.path.abspath(self._getPath('clustal_Omega.aln'))
        # Muscle
        elif programName == MUSCLE:
            outFile = os.path.abspath(self._getPath('muscle.aln'))

        elif programName == MAFFT:
            outFile = os.path.abspath(self._getPath('mafft.aln'))
        
        cline = getMultipleAlignmentCline(programName, inFile, outFile, extraArgs=self.getExtraArgs())
        self.runJob(cline, '')

        if programName == MUSCLE:
            outSeqDic = parseFasta(outFile)
        else:
            outSeqDic = parseAlnFile(outFile)

        outSeqs = SetOfSequencesChem.create(outputPath=self._getPath())
        for seqId in outSeqDic:
            outSeqs.append(Sequence(name=seqId, sequence=outSeqDic[seqId], id=seqId))

        outSeqs.setAlignmentFileName(os.path.relpath(outFile))
        outSeqs.setAligned(True)
        self._defineOutputs(outputSequences=outSeqs)

        # EMBOSS format
        if self.additionalFormat:
            embossSelectedFormat = EMBOSS_FORMATS[self.getEnumText('embossFormats')]
            embossFile = self._getPath('{}.{}'.format(programName.lower(), embossSelectedFormat))
            embossFile = outSeqs.convertEMBOSSformat(embossSelectedFormat, embossFile)


    def _validate(self):
        errors = []
        return errors

    def getExtraArgs(self):
        args = self.extraParams.get()
        if args.strip() == '':
            if self.getEnumText('programList') == MUSCLE:
                args += '-align'
            elif self.getEnumText('programList') == MAFFT:
                args += '--auto'
            elif self.getEnumText('programList') == CLUSTALO:
                args += '--auto'
        return args.strip()

