# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

import glob
import os


from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.utils.path import copyFile
from pyworkflow.protocol.params import PathParam, StringParam, BooleanParam
from pwchem import Plugin
from pwem.objects import Sequence, SetOfSequences

class ProtChemImportSetOfSequences(EMProtocol):
    """Import a set of sequences either from a combined fasta or from multiple fasta files in a directory
    """
    _label = 'Import set of sequences'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('multiple', BooleanParam, default=True,
                      label='Each file is a sequence')
        form.addParam('filesPath', PathParam, condition='multiple',
                      label="Sequences files directory: ",
                      help="Directory with the files you want to import.\n\n"
                           "The path can also contain wildcards to select"
                           "from several folders. \n\n"
                           "Examples:\n"
                           "  ~/project/data/day??_files/\n"
                           "Each '?' represents one unknown character\n\n"
                           "  ~/project/data/day*_files/\n"
                           "'*' represents any number of unknown characters\n\n"
                           "  ~/project/data/day##_files/\n"
                           "'##' represents two digits that will be used as "
                           "file ID\n\n"
                           "NOTE: wildcard characters ('*', '?', '#') "
                           "cannot appear in the actual path.)")
        form.addParam('filesPattern', StringParam,  condition='multiple',
                      label='Pattern',
                      default="*",
                      help="Pattern of the files to be imported.\n\n"
                           "The pattern can contain standard wildcards such as\n"
                           "*, ?, etc, or special ones like ### to mark some\n"
                           "digits in the filename as ID.\n\n"
                           "NOTE: wildcards and special characters "
                           "('*', '?', '#', ':', '%') cannot appear in the "
                           "actual path.\n\n"
                           "You may create small molecules from Smiles (.smi), Tripos Mol2 (.mol2), "
                           "SDF (.sdf), Maestro (.mae, .maegz), or PDB blocks (.pdb)")

        form.addParam('filePath', PathParam, condition='not multiple',
                      label='Sequences file: ',
                      help='Fasta file containing a set of sequences to import')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('importStep')

    def importStep(self):
        if not self.multiple.get():
            fnFasta = self._getExtraPath(os.path.split(self.filePath.get())[1])
            copyFile(self.filePath.get(), fnFasta)

            fastaDic = self.parseFasta(fnFasta)

        else:
            fastaDic = {}
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnFasta = self._getExtraPath(os.path.split(filename)[1])
                copyFile(filename, fnFasta)
                fastaDic.update(self.parseFasta(fnFasta, combined=False))


        outputSequences = SetOfSequences().create(outputPath=self._getPath())
        for seqId in fastaDic:
            newSeq = Sequence(name=seqId, sequence=fastaDic[seqId], id=seqId)
            outputSequences.append(newSeq)

        self._defineOutputs(outputSequences=outputSequences)


    def parseFasta(self, fastaFn, combined=True):
        '''Parses a combined fasta file and returns a dictionary {seqId: sequence}'''
        fastaDic = {}
        seqId, seq = '', ''
        with open(fastaFn) as f:
            for line in f:
                if line.startswith('>'):
                    if seqId != '':
                        fastaDic[seqId] = seq
                        if not combined:
                            #Return just the first sequence if asked for a single fasta
                            return fastaDic
                    seqId = line.strip()[1:]
                elif line.strip() != '':
                    seq += line.strip()
        fastaDic[seqId] = seq
        return fastaDic
