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

from Bio import SeqIO

from pyworkflow.utils.path import createLink
from pwem.objects.data import Sequence
from pwchem.objects import SetOfDatabaseID

def copyFastaSequenceAndRead(protocol):
    outFileName = protocol._getExtraPath("sequence.fasta")
    if isinstance(protocol.inputSeq.get(), Sequence):
        fh = open(outFileName, "w")
        fh.write(">%s\n" % protocol.inputSeq.get().getSeqName())
        fh.write("%s\n" % protocol.inputSeq.get().getSequence())
        fh.close()
    elif isinstance(protocol.inputSeq.get(), SetOfDatabaseID):
        obj = protocol.inputSeq.get().getFirstItem()
        createLink(obj._uniprotFile.get(), outFileName)
    record = SeqIO.read(outFileName, "fasta")
    return str(record.seq)

def checkInputHasFasta(protocol):
    errors = []
    if isinstance(protocol.inputSeq.get(), SetOfDatabaseID):
        if len(protocol.inputSeq.get()) != 1:
            errors.append("The input list can only have a single sequence")
        obj = protocol.inputSeq.get().getFirstItem()
        if not hasattr(obj, "_uniprotFile"):
            errors.append("The input list does not have a sequence file")
    return errors