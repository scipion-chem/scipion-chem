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
    """
    AI Generated:

    This protocol is used to import a set of biological sequences and generate a SetOfSequences object.

    Sequences can be imported either from local FASTA files or downloaded directly from external
    biological databases such as UniProt or ENA. The protocol supports both single and multiple
    sequence imports and automatically converts them into standardized Sequence objects for
    downstream analysis.

    Overview
    --------
    The protocol provides two main import modes:
    - File-based import: load sequences from FASTA files (single or multiple files in a directory)
    - Database-based import: download sequences using database identifiers

    Input
    -----
    fromFile:
        If True, sequences are imported from local files.
        If False, sequences are downloaded from a database.

    From file
    ---------
    multiple:
        If True, each file in a directory is treated as an independent sequence.
        If False, a single FASTA file containing multiple sequences is imported.

    filesPath:
        Directory containing FASTA files (used when multiple=True).

    filesPattern:
        File pattern used to filter input files (e.g. "*.fasta").

    filePath:
        Path to a single FASTA file (used when multiple=False).

    seqType:
        Type of sequences:
        - Protein
        - Nucleotide

    From database
    --------------
    database:
        External database used for download:
        - UniProt
        - ENA

    inputListID:
        List of sequence identifiers to download (one per line).

    Workflow
    --------
    File import mode:
    1. Read FASTA file(s) from disk
    2. Copy files to working directory
    3. Parse sequences into Sequence objects
    4. Add them to a SetOfSequences container

    Database import mode:
    1. Read list of database IDs
    2. Build download URLs (UniProt or ENA REST API)
    3. Retrieve FASTA files from remote server
    4. Parse downloaded FASTA files into Sequence objects
    5. Store results in SetOfSequences

    Output
    ------
    outputSequences:
        SetOfSequences containing all imported sequences as Sequence objects.

    Key Features
    ------------
    - Supports both local and remote sequence import
    - Handles multiple FASTA files in batch mode
    - Integrates with UniProt and ENA REST APIs
    - Automatically converts FASTA into structured Sequence objects
    - Compatible with Scipion EM workflow system

    Notes
    -----
    - Protein vs nucleotide type affects internal parsing behavior
    - Download failures are handled gracefully with warning messages
    - Each sequence is stored with standardized identifiers
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
            newSeq = Sequence(fileName=outFn)
            newSeq.importFromFile(outFn, isAmino=isAmino)
            seqName = os.path.splitext(os.path.basename(outFn))[0]
            newSeq.setSeqName(seqName)
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