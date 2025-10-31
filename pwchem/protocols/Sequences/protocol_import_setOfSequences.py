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

import glob, os, urllib

from pyworkflow.utils.path import copyFile
from pyworkflow.protocol import params

from pwem.protocols import EMProtocol
from pwem.objects import Sequence, SetOfSequences

class ProtChemImportSetOfSequences(EMProtocol):
    """Import a set of sequences either from a combined fasta or from multiple fasta files in a directory
    """
    _label = 'Import set of sequences'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('fromFile', params.BooleanParam, default=True, label='Import from file')
        group = form.addGroup('From file', condition='fromFile')
        group.addParam('multiple', params.BooleanParam, default=True, condition='fromFile',
                      label='Each file is a sequence')
        group.addParam('filesPath', params.PathParam, condition='fromFile and multiple',
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
        group.addParam('filesPattern', params.StringParam,  condition='fromFile and multiple',
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

        group.addParam('filePath', params.PathParam, condition='fromFile and not multiple',
                      label='Sequences file: ',
                      help='Fasta file containing a set of sequences to import')

        group.addParam('seqType', params.EnumParam, default=0, condition='fromFile',
                       choices=['Protein', 'Nucleotide'], display=params.EnumParam.DISPLAY_HLIST,
                       label='Type of sequences: ')

        group = form.addGroup('From database', condition='not fromFile')
        group.addParam('database', params.EnumParam, default=0, condition='not fromFile',
                       label='Input  database: ', choices=['Uniprot', 'ENA'])
        group.addParam('inputListID', params.TextParam, label='List of Database Ids:',
                       help="List of input Ids for downloading. One per line.")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.fromFile:
            self._insertFunctionStep('importStep')
        else:
            self._insertFunctionStep('downloadStep')

    def importStep(self):
        outputSequences = SetOfSequences().create(outputPath=self._getPath())
        if not self.multiple.get():
            fnFasta = self._getExtraPath(os.path.split(self.filePath.get())[1])
            copyFile(self.filePath.get(), fnFasta)

            outputSequences.importFromFile(fnFasta, isAmino=self.seqType.get() == 0)
        else:
            for filename in glob.glob(os.path.join(self.filesPath.get(), self.filesPattern.get())):
                fnFasta = self._getExtraPath(os.path.split(filename)[1])
                copyFile(filename, fnFasta)
                outputSequences.importFromFile(fnFasta, isAmino=self.seqType.get() == 0)

        self._defineOutputs(outputSequences=outputSequences)

    def downloadStep(self):
        outputSequences = SetOfSequences().create(outputPath=self._getPath())
        outFns = []

        for item in self.inputListID.get().split('\n'):
            inId = item.strip()
            print("Processing %s" % inId)
            fnFasta = self._getExtraPath("%s.fasta" % inId)

            url = self.getUrl(inId)
            if not os.path.exists(fnFasta):
                try:
                    urllib.request.urlretrieve(url, fnFasta)
                    outFns.append(fnFasta)
                except: # The library raises an exception when the web is not found
                    print('Could not download {} from {} database'.
                          format(inId, self.getEnumText('database')))

        isAmino = self.database.get() in [1]
        for outFn in outFns:
            newSeq = Sequence()
            newSeq.importFromFile(outFn, isAmino=isAmino)
            outputSequences.append(newSeq)

        self._defineOutputs(outputSequences=outputSequences)
        self._defineSourceRelation(self.inputListID, outputSequences)

    def getUrl(self, dbId):
        if self.database.get() == 0:
            return "https://rest.uniprot.org/uniprotkb/%s.fasta" % dbId
        elif self.database.get() == 1:
            return "https://www.ebi.ac.uk/ena/browser/api/fasta/%s" % dbId

    def _warnings(self):
        ws = []
        return ws