# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Carlos Oscar Sorzano (coss@cnb.csic.es)
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

# General imports
import os, shutil, json, requests, time, subprocess, sys
# import glob
import random as rd
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select
from Bio.PDB.SASA import ShrakeRupley
#from sklearn.cluster import DBSCAN

# Scipion em imports
from pwem.convert import AtomicStructHandler
from pwem.convert.atom_struct import cifToPdb
from pwem.objects.data import Sequence, Object, String, Integer, Float

# Plugin imports
from pwchem.constants import PML_SURF_EACH, PML_SURF_STR, OPENBABEL_DIC
from pwchem import Plugin as pwchemPlugin
from pwchem.utils.scriptUtils import makeSubsets, performBatchThreading

confFirstLine = {'.pdb': 'REMARK', '.pdbqt': 'REMARK',
                 '.mol2': '@<TRIPOS>MOLECULE'}

RESIDUES3TO1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

RESIDUES1TO3 = {v: k for k, v in RESIDUES3TO1.items()}

################# Generic function utils #####################
def insistentExecution(func, *args, maxTimes=5, sleepTime=0, verbose=False):
  """
  ### This function will try to run the given function, and if it fails,
  ### it will retry it the set number of times until it works or has done it the set number of times.

  #### Params:
  - func (function): Function to execute.
  - *args: Optional. Arguments passed to the given function.
  - maxTime (int): Optional. Max number or retries. Default: 5.
  - sleepTime (int): Optional. Number of seconds to wait between each attempt.

  #### Example:
  getUrlContent = lambda url: urlopen(url).read().decode('utf-8')
  content = insistentExecution(getUrlContent, url, maxTimes=10)
  """
  # Printing message if verbose mode is active
  if verbose:
    print("Attempting to run function. Retries left: ", maxTimes-1)
    sys.stdout.flush()
    
  # Attept to run the function
  try:
    return func(*args)
  except Exception as e:
    # Saving exception message
    exception = e
  
  # If function gets to this point, execution failed
  # Show message if verbose is set
  if verbose:
    print("Execution failed with message:\n", exception)
  # Check current number of retries
  if maxTimes > 1:
    if verbose:
      print("Retrying...")
      sys.stdout.flush()
    # If there is at least one retry left, call function again with one less retry
    if sleepTime > 0:
      # Sleep sleepTime seconds if it is greater than 0 seconds
      time.sleep(sleepTime)
    return insistentExecution(func, *args, maxTimes=maxTimes-1, sleepTime=sleepTime, verbose=verbose)
  else:
    print("------------------------- RAISE ERROR -------------------------")
    if verbose:
      print("All executions failed. Re-raising exception.")
      sys.stdout.flush()
    # If max number of retries was fulfilled, raise exception
    raise exception

def insistentRun(protocol, programPath, progArgs, nMax=5, **kwargs):
  i, finished = 1, False
  while not finished and i <= nMax:
    try:
      protocol.runJob(programPath, progArgs, **kwargs)
      finished = True
    except Exception:
      i += 1
      time.sleep(1)
  if i > 1 and i <= nMax:
    print('Program {} run without error after {} trials'.format(programPath, i))
  elif i > nMax:
    print('Program {} could not be run without error after {} trials'.format(programPath, nMax))

def getVarName(var):
  return [i for i, a in locals().items() if a == var][0]

def getBaseName(file):
  return os.path.splitext(os.path.basename(file.strip()))[0]

def getLigCoords(ASFile, ligName):
  """ Return the coordinates of the ligand specified in the atomic structure file. """
  if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
    pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
    parser = PDBParser().get_structure(pdb_code, ASFile)
  elif ASFile.endswith('.cif'):
    pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
    parser = MMCIFParser().get_structure(pdb_code, ASFile)
  else:
    print('Unknown AtomStruct file format')
    return

  coords = []
  for model in parser:
    for chain in model:
      for residue in chain:
        if residue.resname == ligName:
          for atom in residue:
            coords.append(list(atom.get_coord()))
  return coords


def downloadPDB(pdbID, structureHandler=None, outDir='/tmp/'):
  if not structureHandler:
    structureHandler = AtomicStructHandler()

  url = "https://www.rcsb.org/structure/" + str(pdbID)
  try:
    response = requests.get(url)
  except Exception:
    raise Exception("Cannot connect to PDB server")
  if (response.status_code >= 400) and (response.status_code < 500):
    raise Exception("%s is a wrong PDB ID" % pdbID)
  fileName = structureHandler.readFromPDBDatabase(os.path.basename(pdbID), dir=outDir)
  return fileName


def getRawPDBStr(pdbFile, ter=True):
  outStr = ''
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
    f.write(getRawPDBStr(pdbFile, ter))
  return outFile


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
  if j[9] != '':
    j[9] = str('%6.2f' % (float(j[9]))).rjust(6)  # occ
  else:
    j[9] = j[9].rjust(6)
  if j[10] != '':
    j[10] = str('%6.2f' % (float(j[10]))).ljust(6)  # temp
  else:
    j[10] = j[10].ljust(6)
  if j[11] != '':
    j[11] = str('%8.3f' % (float(j[11]))).rjust(10)
  else:
    j[11] = j[11].rjust(10)
  j[12] = j[12].rjust(2)  # elname
  return "%s%s %s %s %s%s    %s%s%s%s%s%s%s\n" % \
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


def mergePDBs(fn1, fn2, fnOut, hetatm2=False):
  with open(fnOut, 'w') as f:
    with open(fn1) as f1:
      for line in f1:
        f.write(line)
    with open(fn2) as f2:
      for line in f2:
        if hetatm2 and line.startswith('ATOM'):
          sLine = splitPDBLine(line)
          sLine[0] = 'HETATM'
          line = writePDBLine(sLine)
        f.write(line)


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
  sampling = [a / 10 for a in range(1, 10)]
  colors = []
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
    surfAtomIds = str(list(map(int, pock.getDecodedCAtoms()))).replace(' ', '')
    surfaceStr += PML_SURF_EACH.format(pId, colors[i], pId, surfAtomIds, pId, pId)

  return PML_SURF_STR.format(pdbFile, surfaceStr)


def writeSurfPML(pockets, pmlFileName):
  with open(pmlFileName, 'w') as f:
    f.write(createSurfacePml(pockets))


def pdbqt2other(protocol, pdbqtFile, otherFile):
  '''Convert pdbqt to pdb or others using openbabel (better for AtomStruct)'''
  inName, inExt = os.path.splitext(os.path.basename(otherFile))
  if not inExt in ['.pdb', '.mol2', '.sdf', '.mol']:
    inExt, otherFile = 'pdb', otherFile.replace(inExt, '.pdb')

  args = ' -ipdbqt {} -o{} -O {}'.format(os.path.abspath(pdbqtFile), inExt[1:], otherFile)
  runOpenBabel(protocol=protocol, args=args, popen=True)
  return os.path.abspath(otherFile)


def convertToSdf(protocol, molFile, sdfFile=None, overWrite=False):
  '''Convert molecule files to sdf using openbabel'''
  if molFile.endswith('.sdf'):
    return molFile
  if not sdfFile:
    baseName = os.path.splitext(os.path.basename(molFile))[0]
    outDir = os.path.abspath(protocol._getTmpPath())
    sdfFile = os.path.abspath(os.path.join(outDir, baseName + '.sdf'))
  else:
    baseName = os.path.splitext(os.path.basename(sdfFile))[0]
    outDir = os.path.abspath(os.path.dirname(sdfFile))
  if not os.path.exists(sdfFile) or overWrite:
    args = ' -i "{}" -of sdf --outputDir "{}" --outputName {} --overWrite'.format(os.path.abspath(molFile),
                                                                                  os.path.abspath(outDir), baseName)
    pwchemPlugin.runScript(protocol, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir, popen=True)
  return sdfFile


def runOpenBabel(protocol, args, cwd='/tmp', popen=False):
  pwchemPlugin.runOPENBABEL(protocol=protocol, args=args, cwd=cwd, popen=popen)


def splitConformerFile(confFile, outDir):
  fnRoot, ext = os.path.splitext(os.path.split(confFile)[1])
  if '_prep_' in fnRoot:
    fnRoot = fnRoot.split('_prep')[0]
  elif '_conformers' in fnRoot:
    fnRoot = fnRoot.split('_conformers')[0]
  iConf, lastRemark, towrite = 1, True, ''
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


def appendToConformersFile(confFile, newFile, outConfFile=None, beginning=True):
  '''Appends a molecule to a conformers file.
    If outConfFile == None, the output conformers file path is the same as as the start'''
  if outConfFile == None:
    iName, iExt = os.path.splitext(confFile)
    outConfFile = confFile.replace(iExt, '_aux' + iExt)
    rename = True

  with open(outConfFile, 'w') as f:
    with open(newFile) as fNew:
      newStr = fNew.read()
    with open(confFile) as fConf:
      confStr = fConf.read()

    if beginning:
      f.write(newStr)
      f.write(confStr)
    else:
      f.write(confStr)
      f.write(newStr)

  if rename:
    os.rename(outConfFile, confFile)
    return confFile

  return outConfFile


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
# ADT grids

def generate_gpf(protFile, spacing, xc, yc, zc, npts, outDir, ligandFns=None, zn_ffFile=None):
  """
    Build the GPF file that is needed for AUTOGRID to generate the electrostatic grid
    """
  protName, protExt = os.path.splitext(os.path.basename(protFile))
  gpf_file = os.path.join(outDir, protName + '.gpf')
  npts = int(round(npts))

  protAtomTypes = parseAtomTypes(protFile)

  if ligandFns == None:
      ligAtomTypes = 'A C HD N NA OA SA'
  else:
      ligAtomTypes = set([])
      for ligFn in ligandFns:
          ligAtomTypes = ligAtomTypes | parseAtomTypes(ligFn)

      ligAtomTypes = protAtomTypes.union(ligAtomTypes)
      ligAtomTypes = ' '.join(sortSet(ligAtomTypes))

  protAtomTypes = ' '.join(sortSet(protAtomTypes))

  with open(os.path.abspath(gpf_file), "w") as file:
    file.write("npts %s %s %s                        # num.grid points in xyz\n" % (npts, npts, npts))
    if zn_ffFile:
        file.write("parameter_file %s                        # force field default parameter file\n" % (zn_ffFile))
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
    file.write("dielectric -0.1465                   # <0, AD4 distance-dep.diel;>0, constant\n")
    if zn_ffFile:
        file.write('''nbp_r_eps 0.25 23.2135 12 6 NA TZ\nnbp_r_eps 2.1   3.8453 12 6 OA Zn\nnbp_r_eps 2.25  7.5914 12 6 SA Zn\nnbp_r_eps 1.0   0.0    12 6 HD Zn\nnbp_r_eps 2.0   0.0060 12 6 NA Zn\nnbp_r_eps 2.0   0.2966 12 6  N Zn''')

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


def parseAtomTypes(pdbqtFile, allowed=None):
  atomTypes = set([])
  if pdbqtFile.endswith('.pdbqt'):
    with open(pdbqtFile) as f:
      for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
          pLine = line.split()
          at = pLine[-1]
          if allowed is None or at in allowed:
            atomTypes.add(at)
  else:
      struct = PDBParser().get_structure("SASAstruct", pdbqtFile)

      for atom in struct.get_atoms():
        atomId = atom.get_id()
        atomTypes.add(removeNumberFromStr(atomId))

  return atomTypes

def sortSet(seti):
  seti = list(seti)
  seti.sort()
  return seti


class CleanStructureSelect(Select):
  def __init__(self, chain_ids, rem_HETATM, rem_WATER, het2keep=[]):
    self.chain_ids = chain_ids
    self.HETATM, self.het2keep= rem_HETATM, het2keep
    self.waters = rem_WATER

  def is_het(self, residue):
    res = residue.id[0]
    return res != " " and res != "W"

  def accept_chain(self, chain):
    return not self.chain_ids or chain.id in self.chain_ids

  def accept_residue(self, residue):
    """ Recognition of heteroatoms - Remove water molecules """
    accept = True
    if self.HETATM and self.is_het(residue) and not residue.resname in self.het2keep:
      accept = False
    elif self.waters and residue.id[0] == 'W':
      accept = False
    return accept


def clean_PDB(struct_file, outFn, waters=False, HETATM=False, chain_ids=None, het2keep=[]):
  """ Extraction of the heteroatoms of .pdb files """
  struct_name = getBaseFileName(struct_file)
  if struct_file.endswith('.pdb') or struct_file.endswith('.ent'):
    struct = PDBParser().get_structure(struct_name, struct_file)
  elif struct_file.endswith('.cif'):
    struct = MMCIFParser().get_structure(struct_name, struct_file)
  else:
    print('Unknown format for file ', struct_file)
    exit()

  io = PDBIO()
  io.set_structure(struct)
  io.save(outFn, CleanStructureSelect(chain_ids, HETATM, waters, het2keep))
  return outFn


def relabelAtomsMol2(atomFile, i=''):
  '''Relabel the atom names so each atom type goes from 1 to x, so if there is only one oxygen named O7,
  it will be renamed to O1'''
  atomCount = {}
  auxFile = atomFile.replace(os.path.basename(atomFile), 'aux{}.mol2'.format(i))
  atomLines = False
  with open(auxFile, 'w') as fOut:
    with open(atomFile) as fIn:
      for line in fIn:
        if atomLines and line.startswith('@<TRIPOS>'):
          atomLines = False

        if atomLines:
          atom = removeNumberFromStr(line.split()[1])
          if atom in atomCount:
            atomCount[atom] += 1
          else:
            atomCount[atom] = 1
          numSize = len(str(atomCount[atom]))
          line = line[:8] + ' ' * (2 - len(atom)) + atom + str(atomCount[atom]).ljust(8) + line[18:]

        if line.startswith('@<TRIPOS>ATOM'):
          atomLines = True

        fOut.write(line)

  shutil.move(auxFile, atomFile)
  return atomFile

def invertDic(d):
    di = {}
    for k, v in d.items():
        di[v] = k
    return di

def relabelMapAtomsMol2(atomFile, i=''):
  '''Relabel the atom names according to a mapping dictionary.'''
  atomCount = {}
  auxFile = atomFile.replace(os.path.basename(atomFile), '{}_aux{}.mol2'.format(os.path.basename(atomFile), i))
  atomLines = False
  with open(auxFile, 'w') as fOut:
    with open(atomFile) as fIn:
      for line in fIn:
        if atomLines and line.startswith('@<TRIPOS>'):
          atomLines = False

        if atomLines:
            atomName = line.split()[1]
            atomType = removeNumberFromStr(atomName)

            atomCount[atomType] = 1 if not atomType in atomCount else atomCount[atomType] + 1
            atomName = atomName if atomType != atomName else atomName + str(atomCount[atomType])

            numSize = len(str(atomCount[atomType]))
            line = line[:8] + ' ' * (2 - len(atomType)) + atomName.ljust(8 + len(atomType)) + line[18:]

        if line.startswith('@<TRIPOS>ATOM'):
            atomLines = True

        fOut.write(line)

  shutil.move(auxFile, atomFile)
  return atomFile


def removeNumberFromStr(s):
  newS = ''
  for i in s:
    if not i.isdigit():
      newS += i
    else:
      pass
  return newS


def getNumberFromStr(s):
  num = ''
  for i in s:
    if not i.isdigit():
      pass
    else:
      num += i
  return num


def calculateDistance(coord1, coord2):
  dist = 0
  if len(coord1) != len(coord2):
    print('ERROR: coordinates of different size')
    return None
  for c1, c2 in zip(coord1, coord2):
    dist += (c1 - c2) ** 2
  return dist ** (1 / 2)


def calculateRMSD(coords1, coords2):
  rmsd = 0
  for c1, c2 in zip(coords1, coords2):
    if len(c1) != len(c2):
      print('ERROR: coordinates of different size')
      return None
    for x1, x2 in zip(c1, c2):
      rmsd += (x1 - x2) ** 2
  return (rmsd / len(coords2)) ** (1 / 2)


def calculateRMSDKeys(coordDic1, coordDic2):
  '''Calculate the RMSD from two dic containing coordinates, using the keys of the
    first one
    '''
  rmsd, count = 0, 0
  for k in coordDic1:
    if k in coordDic2:
      c1, c2 = coordDic1[k], coordDic2[k]
      if len(c1) != len(c2):
        print('ERROR: coordinates of different size')
        return None
      for x1, x2 in zip(c1, c2):
        rmsd += (x1 - x2) ** 2
      count += 1
  return (rmsd / count) ** (1 / 2)


################# UTILS Sequence Object ################

def getSequenceFastaName(sequence):
  '''Return a fasta name for the sequence.
    Priorizes the Id; if None, just "sequence"'''
  return cleanPipeIds(sequence.getId()) if sequence.getId() is not None else 'sequence'


def cleanPipeIds(idStr):
  '''Return a guessed ID when they contain "|"'''
  seqId = idStr.strip()
  if '|' in seqId:
    seqId = seqId.split('|')[1]
  return seqId

def groupConsecutiveIdxs(idxs):
    idxs.sort()
    groups, newGroup = [], []
    for idx in idxs:
      idx = int(idx)
      if not newGroup or idx - 1 == newGroup[-1]:
        newGroup.append(idx)
      else:
        groups.append(newGroup)
        newGroup = [idx]
    groups.append(newGroup)
    return groups

def natural_sort(listi, rev=False):
  import re
  convert = lambda text: int(text) if text.isdigit() else text.lower()
  alphanum_key = lambda key: [convert(c) for c in re.split('([0-9]+)', key)]
  return sorted(listi, key=alphanum_key, reverse=rev)


def fillEmptyAttributes(inputSets):
  '''Fill all items with empty attributes'''
  attributes = getAllAttributes(inputSets)
  for inSet in inputSets:
    for item in inSet.get():
      for attr in attributes:
        if not hasattr(item, attr):
          item.__setattr__(attr, attributes[attr])
  return inputSets


def getAllAttributes(inputSets):
  '''Return a dic with {attrName: ScipionObj=None}'''
  attributes = {}
  for inpSet in inputSets:
    item = inpSet.get().getFirstItem()
    attrKeys = item.getObjDict().keys()
    for attrK in attrKeys:
      if not attrK in attributes:
        value = item.__getattribute__(attrK)
        attributes[attrK] = value.clone()
        attributes[attrK].set(None)
  return attributes


def getBaseFileName(filename):
  return os.path.splitext(os.path.basename(filename))[0]


################# Wizard utils #####################

def getChainIds(chainStr):
  '''Parses a line of json with the description of a chain or chains and returns the ids'''
  chainJson = json.loads(chainStr)
  if 'chain' in chainJson:
    chain_ids = [chainJson["chain"].upper().strip()]
  elif 'model-chain' in chainJson:
    modelChains = chainJson["model-chain"].upper().strip()
    chain_ids = [x.split('-')[1] for x in modelChains.split(',')]
  return chain_ids


def calculate_SASA(structFile, outFile):
  if structFile.endswith('.pdb') or structFile.endswith('.ent'):
    p = PDBParser(QUIET=1)
  elif structFile.endswith('.cif'):
    p = MMCIFParser(QUIET=1)
  struct = p.get_structure("SASAstruct", structFile)

  sr = ShrakeRupley()
  sr.compute(struct, level="R")

  with open(outFile, 'w') as f:
    for model in struct:
      for chain in model:
        chainID = chain.get_id()
        for residue in chain:
          resId = residue.get_id()[1]
          f.write('{}:{}\t{}\n'.format(chainID, resId, residue.sasa))

################# Test utils #####################
def assertHandle(func, *args, cwd='', message=''):
  """
  ### This function runs the given assertion and handles the potential error, showing the protocol's error trace instead of a generic assertion.
  ### Note: Only designed to handle the assertions used by protocols.

  #### Params:
  - func (function): A function (assertion) to be executed.
  - *args: All the arguments passed to the function.
  - cwd (str): Optional. Current working directory. Usually protocol's cwd.
  - message (str): Optional. Custom message to display if no error messages were found in stdout/stderr.

  #### Example:
  assertHandle(self.assertIsNotNone, getattr(protocol, 'outputSet', None), cwd=protocol.getWorkingDir())
  """
  # Defining full path to error log
  stderr = os.path.abspath(os.path.join(cwd, 'logs', 'run.stderr'))
  stdout = os.path.abspath(os.path.join(cwd, 'logs', 'run.stdout'))

  # Attempt to run assertion
  try:
    return func(*args)
  except AssertionError:
    # If assertion fails, show error log
    # Getting error logs (stderr has priority over stdout)
    # Most errors are dumped on stderr, while some others on stdout
    errorMessage = ''
    for stdFile in [stderr, stdout]:
      if os.path.exists(stdFile):
        errorMessage = subprocess.run(['cat', stdFile], check=True, capture_output=True).stdout.decode()
        # Sometimes stderr file exists but it is empty, in those cases, fall back to stdout
        if errorMessage:
          break
    if not errorMessage:
      errorMessage = "Something went wrong with the protocol, but there are no stderr/stdout files right now, try manually opening the project to check it."
      if message:
        errorMessage += f"\n{message}"
    raise AssertionError(f"Assertion {func.__name__} failed for the following reasons:\n\n{errorMessage}")