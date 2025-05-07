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

from pwem.convert.sequence import alignClustalSequences
from pwchem.objects import SetOfDatabaseID
from pwchem.constants import BIOCONDA_DIC
from pwchem import Plugin

EMBOSS_FORMATS = {'Fasta': 'fasta', 'Clustal': 'aln', 'Wisconsin Package GCG 9.x and 10.x': 'gcg', 'GCG 8.x': 'gcg8',
                  'SwisProt': 'sw', 'NCBI': 'ncbi',
                  'NBRF (PIR)': 'pir', 'Intelligenetics': 'ig', 'CODATA': 'codata', 'DNA strider': 'strider',
                  'ACeDB': 'acedb', '"gap" program in the Staden package': 'experiment', 'Plain sequence': 'plain',
                  'Fitch': 'fitch', 'PHYLIP interleaved': 'phylip3', 'ASN.1': 'asn1', 'Hennig86': 'hennig86',
                  'Mega': 'mega', 'Meganon': 'meganon', 'Nexus/PAUP': 'nexus', 'Nexusnon/PAUPnon': 'nexusnon',
                  'Jackknifer': 'jackknifer','Jackknifernon': 'jackknifernon', 'Treecon': 'treecon',
                  'EMBOSS sequence object report': 'debug'}

def convertEMBOSSformat(inputAlignmentFile, embossFormat, outputAligmentFile):
  '''Connvert alignment files. Options in EMBOSS_FORMATS dictionary
  Returns a command line which must be executed into the scipion environment'''
  cl_run = '%s && ' % (Plugin.getEnvActivationCommand(BIOCONDA_DIC))
  cl_run += 'seqret -sequence {} -osformat2 {} {}'. \
    format(inputAlignmentFile, embossFormat, outputAligmentFile)
  return cl_run

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


def clustalOmegaAlignSequences(sequencesSet, seqFileName, outputFileName):
  sequencesSet.exportToFile(seqFileName)
  cline = alignClustalSequences(seqFileName, outputFileName)
  return cline


def muscleAlignmentSequences(sequencesSet, seqFileName, outputFileName):
  from pwem.convert.sequence import alignMuscleSequences
  sequencesSet.exportToFile(seqFileName)
  cline = alignMuscleSequences(seqFileName, outputFileName)
  return cline

def pairwiseAlign(seq1, seq2, outPath, seqName1=None, seqName2=None, force=False):
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
    oriFasta = os.path.abspath(os.path.join(outDir, 'original.fasta'))
    with open(oriFasta, "w") as f:
      f.write(('>{}\n{}\n>{}\n{}\n'.format(seqName1, seq1, seqName2, seq2)))

    ext = os.path.splitext(outPath)[1][1:]
    if ext in ['fa', 'fasta']:
        fmt = 'fa'
    elif ext == 'aln':
        fmt = 'clu'

    # Alignment
    activate_env_line = '%s && ' % (Plugin.getEnvActivationCommand(BIOCONDA_DIC))
    cline = '{} {} --outfmt={}'.format(activate_env_line, alignClustalSequences(oriFasta, outPath), fmt)
    if force:
        cline += ' --force'
    subprocess.check_call(cline, cwd=outDir, shell=True)

def parseFasta(fastaFile):
    seqDic = {}
    with open(fastaFile) as f:
        for line in f:
            if line.startswith('>'):
                seqId = line.strip()[1:]
                seqDic[seqId] = ''
            elif line.strip() != '':
                seqDic[seqId] += line.strip()
    return seqDic

def parseAlnFile(alnFile):
    seqDic = {}
    with open(alnFile) as f:
        f.readline()
        for line in f:
            if line.strip() and not line.startswith(' '):
                id, seq = line.strip().split()[:2]
                if id in seqDic:
                    seqDic[id] += seq
                else:
                    seqDic[id] = seq
    return seqDic


def calculateIdentity(alignFile):
    if os.path.splitext(alignFile)[1] in ['.fa', '.fasta']:
        seqDic = parseFasta(alignFile)
        if len(seqDic) == 2:
            hits = 0
            seqIds = list(seqDic.keys())
            for i, j in zip(seqDic[seqIds[0]], seqDic[seqIds[1]]):
                if i == j:
                    hits += 1
            ident = hits / len(seqDic[seqIds[0]])

            return round(100 * ident, 2)


def fastFastaExport(seqSet, outFasta):
  lines = []
  for seq in seqSet:
    lines.append(f'>{seq.getSeqName()}\n{seq.getSequence()}\n')

  with open(outFasta, 'w') as f:
    f.write(''.join(lines))


CLUSTALO, MUSCLE, MAFFT = 'CLUSTAL_OMEGA', 'MUSCLE', 'MAFFT'
def getMultipleAlignmentCline(programName, inputFasta, outputFile, extraArgs=None):
  programName = programName.upper()
  if extraArgs is None or extraArgs.strip() == '':
    if programName == MUSCLE:
      extraArgs = '-align'
    elif programName == MAFFT:
      extraArgs = '--auto'
    elif programName == CLUSTALO:
      extraArgs = '--auto'
    else:
      extraArgs = ''

  cline = '%s && ' % (Plugin.getEnvActivationCommand(BIOCONDA_DIC))
  # Clustal Omega
  if programName == CLUSTALO:
    cline += 'clustalo -i {} {} -o {} --outfmt=clu'.format(inputFasta, extraArgs, outputFile)

  # Muscle
  elif programName == MUSCLE:
    cline += 'muscle {} {} -output {}'.format(extraArgs, inputFasta, outputFile)

  elif programName == MAFFT:
    cline += 'mafft {} --clustalout {} > {}'.format(extraArgs, inputFasta, outputFile)
  else:
    cline += 'exit'

  return cline

