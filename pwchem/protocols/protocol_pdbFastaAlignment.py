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

# import lxml.etree as ET
import os
import sys
import urllib.request

import Bio
from Bio.Align.Applications import ClustalwCommandline
from Bio import AlignIO

from pwem.protocols import EMProtocol
import pyworkflow.object as pwobj
from pyworkflow.protocol.params import StringParam, EnumParam, FileParam, PointerParam
from pwchem.objects import DatabaseID, SetOfDatabaseID, SequenceFasta, ClustalAln

class ProtChemPdbFastaAlignment(EMProtocol):
    """Download a fasta from Uniprot and insert variants"""
    _label = 'pdb fasta alignment'

    def _defineParams(self, form):
        form.addSection(label='Input')

        # form.addParam('inputPDB', PointerParam, pointerClass="AtomStruct",
        #                label='PDB:', allowsNull=False)

        form.addParam('inputPDBid', StringParam,
                      label='PDB ID:', allowsNull=False,
                      help="PDB ID for the protein")

        form.addParam('inputchain', StringParam,
                      label='PDB Chain Id:', allowsNull=False,
                      help="chain ID for the query sequence in the PDB")

        form.addParam('inputProtName', StringParam,
                      label='Protein Name:', allowsNull=False,
                      help="Protein Name")

        form.addParam('fastaFile', FileParam,
                    label='FASTA sequence:', allowsNull=False,
                    help="Subject sequence")

    def _insertAllSteps(self):
        self._insertFunctionStep('downLoadPdb')
        self._insertFunctionStep('alignmentFasta')

    def downLoadPdb(self):
        outputDatabaseID = SetOfDatabaseID().create(outputPath=self._getPath())
        pdbId = self.inputPDBid.get().lower()
        urlPdb = "https://files.rcsb.org/download/%s.pdb" % pdbId
        self.fnPdb = self._getPath("%s.pdb" % pdbId)

        if not os.path.exists(self.fnPdb):
            print("Fetching pdb: %s"%urlPdb)

            try:
                urllib.request.urlretrieve(urlPdb,self.fnPdb)
            except: # The library raises an exception when the web is not found
                pass

        chain = self.inputchain.get()
        f_pdb = open(self.fnPdb, "r")

        pos_residue = 1
        aa_sequence = []
        for linePDB in f_pdb:
            rst_lineaPDB = linePDB.rstrip()
            palabras = rst_lineaPDB.split()
            for posi_palabra in range(0, len(palabras)):
                if palabras[posi_palabra] == 'ATOM':
                    if posi_palabra == 0 and palabras[posi_palabra + 4] == chain and palabras[
                        posi_palabra + 5] != pos_residue:
                        pos_residue = palabras[posi_palabra + 5]
                        aa_sequence.append(palabras[posi_palabra + 3])
                        # print(palabras[posi_palabra + 3])

        aa_code = dict()
        aa_code['ALA'] = 'A'
        aa_code['ARG'] = 'R'
        aa_code['ASN'] = 'N'
        aa_code['ASP'] = 'D'
        aa_code['CYS'] = 'C'
        aa_code['GLN'] = 'Q'
        aa_code['GLU'] = 'E'
        aa_code['GLY'] = 'G'
        aa_code['HIS'] = 'H'
        aa_code['ILE'] = 'I'
        aa_code['LEU'] = 'L'
        aa_code['LYS'] = 'K'
        aa_code['MET'] = 'M'
        aa_code['PHE'] = 'F'
        aa_code['PRO'] = 'P'
        aa_code['SER'] = 'S'
        aa_code['THR'] = 'T'
        aa_code['TRP'] = 'W'
        aa_code['TYR'] = 'Y'
        aa_code['VAL'] = 'V'

        # Make file to write sequence
        self.seqExtractedPdb = self._getPath("sequenceFromPDB_%s.fasta" % pdbId)
        file = open(self.seqExtractedPdb, "w")
        file.close()

        aa_sequence_code1 = ''
        for aminoacid in range(0, len(aa_sequence)):
            for aa3, aa1 in aa_code.items():
                if aa_sequence[aminoacid] == aa3:
                    aa_sequence_code1 = aa_sequence_code1 + aa1

        file = open(self.seqExtractedPdb, 'a')
        file.write(('>SequencePDB_%s.fasta' % pdbId) + "\n")
        file.write(aa_sequence_code1)
        file.close()

        seqPdb = SequenceFasta()
        seqPdb.setFileName(self.seqExtractedPdb)

        self._defineOutputs(outputSeqPDB=seqPdb)

    def alignmentFasta(self):
        pdbId = self.inputPDBid.get().lower()
        fastaName = self.inputProtName.get()
        fn_fasta = self.fastaFile.get()

        # Fasta file should have 2 lines, header and sequence
        pair_fasta = self._getPath('pair_' + pdbId + fastaName + '.fasta')
        total_fasta = open(pair_fasta, "w")

        fasta_1 = open(self.seqExtractedPdb, "r")
        for ele1 in fasta_1:
            linea_ele1 = ele1.rstrip()
            total_fasta.write(linea_ele1 + "\n")

        fasta_2 = open(fn_fasta, "r")
        for ele2 in fasta_2:
            linea_ele2 = ele2.rstrip()
            total_fasta.write(linea_ele2 + "\n")

        total_fasta.close()

        # Alignment
        in_file = pair_fasta
        out_file = self._getPath("clustalw_" + pdbId + fastaName + ".aln")

        clustalw_cline = ClustalwCommandline("clustalw", infile=in_file, outfile=out_file)
        print(clustalw_cline)

        stdout, stderr = clustalw_cline()
        print(stdout)
        print(stderr)

        pairFasta = SequenceFasta()
        pairFasta.setFileName(pairFasta)
        alnClustal = ClustalAln()
        alnClustal.setFileName(out_file)

        self._defineOutputs(outputPairFasta=pairFasta)
        self._defineOutputs(outputClustalw=alnClustal)

    def _validate(self):
        errors=[]
        return errors