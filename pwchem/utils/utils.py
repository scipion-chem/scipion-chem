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

from pwem.convert import AtomicStructHandler
from pwem.objects.data import Sequence, Object, String, Integer, Float
from ..constants import *
from pwchem import Plugin as pwchemPlugin
import random as rd
import os
import numpy as np

confFirstLine = {'.pdb': 'REMARK', '.pdbqt':'REMARK',
                 '.mol2': '@<TRIPOS>MOLECULE'}


def getRawPDBStr(pdbFile, ter=True):
    outStr=''
    with open(pdbFile) as fIn:
        for line in fIn:
            if line.startswith('ATOM') or line.startswith('HETATM'):
                outStr += line
            elif ter and line.startswith('TER'):
                outStr += line
    return outStr

def writeRawPDB(pdbFile, outFile, ter=True):
    '''Creates a new pdb with only the ATOM and HETATM lines'''
    with open(outFile, 'w') as f:
        getRawPDBStr(pdbFile)

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

def splitPDBLine(line, rosetta=False):
    if line.startswith(("ATOM", "HETATM")):
        atomType = line[0:6].strip()
        atomSerialNumber = line[6:11].strip()
        atomName = line[12:16].strip()
        resName = line[17:20].strip()
        chain = line[21].strip()
        resNumber = line[22:26].strip()
        coorX = line[30:38].strip()
        coorY = line[38:46].strip()
        coorZ = line[46:54].strip()
        occupancy = line[54:60].strip()
        temperatureFact = line[60:66].strip()
        if rosetta:
            segmentIdentifier = line[70:76].strip()
            elementSymbol = line[76:79].strip()
        else:
            segmentIdentifier = line[72:76].strip()
            elementSymbol = line[76:78].strip()
        return [atomType, atomSerialNumber, atomName, resName, chain, resNumber,
                coorX, coorY, coorZ, occupancy, temperatureFact, segmentIdentifier, elementSymbol]
    else:
        return None

def mergePDBs(fn1, fn2, fnOut):
    with open(fnOut, 'w') as f:
        with open(fn1) as f1:
            f.write('\n'.join(f1.readlines()[:-1]))
        with open(fn2) as f2:
            f.write(f2.read())

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


def runOpenBabel(protocol, args, cwd):
    pwchemPlugin.runOPENBABEL(protocol=protocol, args=args, cwd=cwd)


def splitConformerFile(confFile, outDir):
    _, ext = os.path.splitext(confFile)
    fnRoot = os.path.split(confFile)[1].split('_')[0]
    iConf, lastRemark, towrite = 2, True, ''
    with open(confFile) as fConf:
        for line in fConf:
            if line.startswith(confFirstLine[ext]):
                if lastRemark:
                    towrite += line
                else:
                    newFile = os.path.join(outDir, '{}-{}{}'.format(fnRoot, iConf, ext))
                    writeFile(towrite, newFile)
                    towrite, lastRemark = line, True
                    iConf += 1
            else:
                towrite += line
                lastRemark = False
    newFile = os.path.join(outDir, '{}-{}{}'.format(fnRoot, iConf, ext))
    writeFile(towrite, newFile)
    return outDir


def writeFile(towrite, file):
    with open(file, 'w') as f:
        f.write(towrite)


def getProteinMaxDiameter(protFile):
    protCoords = np.array(getPDBCoords(protFile))
    minCoords, maxCoords = protCoords.min(axis=0), protCoords.max(axis=0)
    return max(maxCoords - minCoords)

def getPDBCoords(pdbFile):
    coords = []
    with open(pdbFile) as f:
      for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
            line = splitPDBLine(line)
            coords.append((float(line[6]), float(line[7]), float(line[8])))
    return coords

##################################################
#ADT grids

def generate_gpf(protFile, spacing, xc, yc, zc, npts, outDir, ligandFns=None):
      """
      Build the GPF file that is needed for AUTOGRID to generate the electrostatic grid
      """

      protName, protExt = os.path.splitext(os.path.basename(protFile))
      gpf_file = os.path.join(outDir, protName + '.gpf')
      npts = int(round(npts))

      protAtomTypes = parseAtomTypes(protFile)
      protAtomTypes = ' '.join(sortSet(protAtomTypes))

      if ligandFns == None:
          ligAtomTypes = 'A C HD N NA OA SA'
      else:
          ligAtomTypes = set([])
          for ligFn in ligandFns:
              ligAtomTypes = ligAtomTypes | parseAtomTypes(ligFn)

          ligAtomTypes = ' '.join(sortSet(ligAtomTypes))



      with open(os.path.abspath(gpf_file), "w") as file:
        file.write("npts %s %s %s                        # num.grid points in xyz\n" % (npts, npts, npts))
        file.write("gridfld %s.maps.fld                # grid_data_file\n" % (protName))
        file.write("spacing %s                          # spacing(A)\n" % (spacing))
        file.write("receptor_types %s     # receptor atom types\n" % (protAtomTypes))
        file.write("ligand_types %s       # ligand atom types\n" % (ligAtomTypes))
        file.write("receptor %s                  # macromolecule\n" % (os.path.abspath(protFile)))
        file.write("gridcenter %s %s %s           # xyz-coordinates or auto\n" % (xc, yc, zc))
        file.write("smooth 0.5                           # store minimum energy w/in rad(A)\n")
        for ligType in ligAtomTypes.split():
            file.write("map %s.%s.map                       # atom-specific affinity map\n" % (protName, ligType))
        file.write("elecmap %s.e.map                   # electrostatic potential map\n" % (protName))
        file.write("dsolvmap %s.d.map                  # desolvation potential map\n" % (protName))
        file.write("dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant")

      return os.path.abspath(gpf_file)


def calculate_centerMass(atomStructFile):
    """
    Returns the geometric center of mass of an Entity (anything with a get_atoms function in biopython).
    Geometric assumes all masses are equal
    """

    try:
        structureHandler = AtomicStructHandler()
        structureHandler.read(atomStructFile)
        center_coord = structureHandler.centerOfMass(geometric=True)
        structure = structureHandler.getStructure()

        return structure, center_coord[0], center_coord[1], center_coord[2]  # structure, x,y,z

    except Exception as e:
        print("ERROR: ", "A pdb file was not entered in the Atomic structure field. Please enter it.", e)
        return

def parseAtomTypes(pdbFile):
    atomTypes = set([])
    with open(pdbFile) as f:
        for line in f:
            if line.startswith('ATOM') or line.startswith('HETATM'):
              pLine = line.split()
              try:
                  at = pLine[12]
              except:
                  at = splitPDBLine(line, rosetta=True)[12]

              atomTypes.add(at)
    return atomTypes


def sortSet(seti):
    seti = list(seti)
    seti.sort()
    return seti