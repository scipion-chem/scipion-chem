# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os, subprocess
from Bio import SeqIO
from pyworkflow.utils.path import createLink
from pwem.objects.data import Sequence
from pwem.convert import cifToPdb, alignClustalSequences

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


def clustalOmegaAlignSequences(SetOfSequences, seqFileName, outputFileName):
  from pwem.convert.sequence import alignClustalSequences
  SetOfSequences.exportToFile(seqFileName)
  cline = alignClustalSequences(seqFileName, outputFileName)
  return cline


def muscleAlignmentSequences(SetOfSequences, seqFileName, outputFileName):
  from pwem.convert.sequence import alignMuscleSequences
  SetOfSequences.exportToFile(seqFileName)
  cline = alignMuscleSequences(seqFileName, outputFileName)
  return cline

def pairwiseAlign(seq1, seq2, outPath, seqName1=None, seqName2=None):
    if issubclass(type(seq1), Sequence):
        seqName1 = seq1.getSeqName()
        seq1 = seq1.getSequence()
    else:
        if not seqName1:
            seqName1 = 'sequence_1'

    if issubclass(type(seq2), Sequence):
        seqName2 = seq2.getSeqName()
        seq2 = seq2.getSequence()
    else:
        if not seqName2:
            seqName2 = 'sequence_2'

    outDir = os.path.dirname(outPath)
    oriFasta = os.path.join(outDir, 'original.fasta')
    with open(oriFasta, "w") as f:
      f.write(('>{}\n{}\n>{}\n{}\n'.format(seqName1, seq1, seqName2, seq2)))

    ext = os.path.splitext(outPath)[1][1:]
    if ext in ['fa', 'fasta']:
        fmt = 'fa'
    elif ext == 'aln':
        fmt = 'clu'

    # Alignment
    cline = str(alignClustalSequences(oriFasta, outPath))
    cline += ' --outfmt={}'.format(fmt)
    subprocess.check_call(cline, cwd=outDir, shell=True)









