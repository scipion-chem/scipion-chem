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
from pwem.objects.data import Sequence, Object, String, Integer, Float
from pwchem.objects import SetOfDatabaseID
from ..constants import *
import random as rd

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

def writeRawPDB(pdbFile, outFile, ter=True):
    '''Creates a new pdb with only the ATOM and HETATM lines'''
    with open(outFile, 'w') as f:
        with open(pdbFile) as fIn:
            for line in fIn:
                if line.startswith('ATOM') or line.startswith('HETATM'):
                    f.write(line)
                elif ter and line.startswith('TER'):
                    f.write(line)

def writePDBLine(j):
    '''j: elements to write in the pdb'''
    j[0] = j[0].ljust(6)  # atom#6s
    j[1] = j[1].rjust(5)  # aomnum#5d
    j[2] = j[2].center(4)  # atomname$#4s
    j[3] = j[3].ljust(3)  # resname#1s
    j[4] = j[4].rjust(1)  # Astring
    j[5] = j[5].rjust(4)  # resnum
    j[6] = str('%8.3f' % (float(j[6]))).rjust(8)  # x
    j[7] = str('%8.3f' % (float(j[7]))).rjust(8)  # y
    j[8] = str('%8.3f' % (float(j[8]))).rjust(8)  # z\
    j[9] = str('%6.2f' % (float(j[9]))).rjust(6)  # occ
    j[10] = str('%6.2f' % (float(j[10]))).ljust(6)  # temp
    if j[11] != '':
        j[11] = str('%8.3f' % (float(j[11]))).rjust(10)
    else:
        j[11] = j[11].rjust(10)
    j[12] = j[12].rjust(2)  # elname
    return "\n%s%s %s %s %s%s    %s%s%s%s%s%s%s" % \
           (j[0], j[1], j[2], j[3], j[4], j[5], j[6], j[7], j[8], j[9], j[10], j[11], j[12])

def getScipionObj(value):
    if isinstance(value, Object):
        return value
    elif isinstance(value, int):
        return Integer(value)
    elif isinstance(value, float):
        return Float(value)
    elif isinstance(value, str):
        return String(value)
    else:
        return None

def setAttribute(obj, label, value):
    if value is None:
        return
    setattr(obj, label, getScipionObj(value))

def copyAttribute(src, dst, label, default=None):
    setAttribute(dst, label, getattr(src, label, default))


def createColorVectors(nColors):
    sampling = [a/10 for a in range(1, 10)]
    colors=[]
    while len(colors) < nColors:
        newColor = rd.sample(sampling, 3)
        if not newColor in colors:
            colors += [newColor]
    return colors

def createSurfacePml(pockets):
    pdbFile = pockets.getProteinFile()
    colors = createColorVectors(len(pockets))
    surfaceStr = ''
    for i, pock in enumerate(pockets):
        pId = pock.getObjId()
        surfAtomIds = str(list(map(int, pock.getDecodedCAtoms()))).replace(' ','')
        surfaceStr += PML_SURF_EACH.format(pId, colors[i], pId, surfAtomIds, pId, pId)

    return PML_SURF_STR.format(pdbFile, surfaceStr)

def writeSurfPML(pockets, pmlFileName):
    with open(pmlFileName, 'w') as f:
        f.write(createSurfacePml(pockets))

