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
import os, shutil, json, requests, time, subprocess, sys, multiprocessing, re
import random as rd
import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, PDBIO, Select, MMCIFIO
from Bio.PDB.SASA import ShrakeRupley

# Scipion em imports
from pwem.convert import AtomicStructHandler
from pwem.objects.data import Sequence, Object, String, Integer, Float, Pointer

# Plugin imports
from ..constants import PML_SURF_EACH, PML_SURF_STR, OPENBABEL_DIC, RDKIT_DIC
from .. import Plugin as pwchemPlugin
from ..utils.scriptUtils import makeSubsets, performBatchThreading

confFirstLine = {'.pdb': 'REMARK', '.pdbqt': 'REMARK',
                 '.mol2': '@<TRIPOS>MOLECULE'}

RESIDUES3TO1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
                'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
                'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
                'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

RESIDUES1TO3 = {v: k for k, v in RESIDUES3TO1.items()}

SEED_RANDOM = '''get_seeded_random()
{
  seed="$1"
  openssl enc -aes-256-ctr -pass pass:"$seed" -nosalt \
    </dev/zero 2>/dev/null
}'''

def checkNormalResidues(sequence):
  return all([res in RESIDUES1TO3 for res in sequence])

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
    print("Attempting to run function. Retries left: ", maxTimes)
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
    if verbose:
      print("All executions failed. Re-raising exception.")
      sys.stdout.flush()
    # If max number of retries was fulfilled, raise exception
    raise exception

def reduceNRandomLines(inFile, n, oDir='/tmp', seed=None):
  tFile = os.path.abspath(os.path.join(oDir, 'temp.txt'))
  precall = '' if seed is None else f'{SEED_RANDOM} && '
  randomArg = '' if seed is None else f' --random-source=<(get_seeded_random {seed}) '
  command = f'{precall}shuf -n {n} {os.path.abspath(inFile)}{randomArg}> {tFile}'
  subprocess.check_call(['bash', '-c', command])
  os.rename(tFile, inFile)

def gunzipFile(gzFile, oFile=None, remove=False):
  if oFile is None:
    filename = '.'.join(getBaseFileName(gzFile).split('.')[:-1])
    oFile = os.path.join(os.path.dirname(gzFile), filename)

  subprocess.check_call(f'gunzip -d {gzFile} -c > {oFile}', shell=True)
  if remove:
    os.remove(gzFile)
  return oFile

def swapColumns(smiFile, swaps=(1,2), oDir='/tmp'):
  smiFile = os.path.abspath(smiFile)
  tFile = os.path.abspath(os.path.join(oDir, 'smis.smi'))
  subprocess.check_call(f"awk '{{ t=${swaps[0]}; ${swaps[0]}=${swaps[1]}; ${swaps[1]}=t; print }}' "
                        f"{smiFile} > {tFile}", shell=True)
  os.rename(tFile, smiFile)

def concatFiles(inFiles, oFile, remove=False, skipHead=0):
  subprocess.check_call(f'awk FNR!={skipHead} {" ".join(inFiles)} > {oFile}', shell=True)
  if remove:
    [os.remove(file) for file in inFiles]

def concatGZFiles(inFiles, oFile, remove):
  cmd = f'zcat {" ".join(inFiles)} > {oFile}'
  subprocess.check_call(cmd, shell=True)
  if remove:
    [os.remove(file) for file in inFiles]

def mergeFiles(inFiles, outFile=None, oDir='', sep='', remove=False):
  inTexts = []
  for inFile in inFiles:
    with open(inFile) as f:
      inTexts.append(f.read())

  outFile = os.path.join(oDir, f'mergedFiles{os.path.splitext(inFile)[-1]}') if outFile is None else outFile
  with open(outFile, 'w') as fo:
    fo.write(sep.join(inTexts))

  if remove:
    [os.remove(inFile) for inFile in inFiles]
  return outFile

def findThreadFiles(filename, directory=None):
  '''Finds the thread files (file_i.txt) in a directory.
  filename: name of the merged file (file.txt)
  directory: directory where the thread files are. If None, will expect filename to be a path
  '''
  if not directory:
    directory, filename = os.path.dirname(filename), os.path.split(filename)[-1]

  basename, ext = os.path.splitext(getBaseFileName(filename))
  pattern = re.compile(rf"{basename}_\d+{ext}")

  matchingFiles = [os.path.join(directory, f) for f in os.listdir(directory) if pattern.match(f)]
  return matchingFiles

def concatThreadFiles(concatFile, inDir=None, remove=True):
  thFiles = findThreadFiles(concatFile, directory=inDir)
  concatFiles(thFiles, concatFile, remove=remove)

def removeThreadDirectories(pattern, inDir=None):
  '''Remove the set of thread directories: pattern_id'''
  if inDir is None:
    inDir = os.path.dirname(pattern)
  inDir = os.path.abspath(inDir)
  basePatern = os.path.basename(pattern)

  for file in os.listdir(inDir):
    if basePatern in file:
      shutil.rmtree(os.path.join(inDir, file))


def organizeThreads(nTasks, nThreads):
  if nTasks > nThreads:
    return [1] * nTasks
  else:
    subsets = [0 for _ in range(nTasks)]
    for i in range(nThreads):
      subsets[i % nTasks] += 1
  return subsets

def insistentRun(protocol, programPath, progArgs, nMax=5, sleepTime=1, popen=False, **kwargs):
  i, finished = 1, False
  while not finished and i <= nMax:
    try:
      if not popen:
        protocol.runJob(programPath, progArgs, **kwargs)
      else:
        subprocess.check_call(programPath + progArgs, shell=True, **kwargs)

      finished = True
    except Exception:
      i += 1
      time.sleep(sleepTime)

  if i > 1 and i <= nMax:
    print('Program {} run without error after {} trials'.format(programPath, i))
  elif i > nMax:
    print('Program {} could not be run without error after {} trials'.format(programPath, nMax))

def getVarName(var):
  return [i for i, a in locals().items() if a == var][0]

def getBaseFileName(file):
  '''From a file path, returns the name of the file without the directory:
  /abc/def/file.txt -> filename.txt'''
  return os.path.basename(file.strip())

def getBaseName(file):
  '''From a file path, returns the name of the file without the directory nor the extension:
    /abc/def/filename.txt -> filename'''
  return os.path.splitext(getBaseFileName(file))[0]

def parseAtomStruct(asFile):
  '''Parse an atom struct using biopython'''
  if asFile.endswith('.pdb') or asFile.endswith('.ent'):
    pdbCode = os.path.basename(os.path.splitext(asFile)[0])
    parser = PDBParser().get_structure(pdbCode, asFile)
  elif asFile.endswith('.cif'):
    pdbCode = os.path.basename(os.path.splitext(asFile)[0])
    parser = MMCIFParser().get_structure(pdbCode, asFile)
  else:
    print('Unknown AtomStruct file format')
    parser = None
  return parser

def isHet(residue):
  res = residue.id[0]
  return res != " " and res != "W"

def getLigCoords(asFile, ligName):
  """ Return the coordinates of the ligand specified in the atomic structure file. """
  parser = parseAtomStruct(asFile)
  if parser:
    coords = []
    for model in parser:
      for residue in model.get_residues():
        if residue.resname == ligName:
          for atom in residue:
            coords.append(list(atom.get_coord()))
  return coords

def downloadUrlFile(url, oDir, trials=3):
  '''Downloads a file from its url and return the local file downloaded'''
  coms = f"wget -P {oDir} {url} -q"
  oFile = os.path.join(oDir, getBaseFileName(url))
  i = 0
  while i < trials and not os.path.exists(oFile):
    subprocess.run(coms, shell=True)
    i += 1
  return oFile

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
      if line.startswith('ATOM') or line.startswith('HETATM') or (ter and line.startswith('TER')):
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

def splitFile(inFile, b=None, n=None, oDir=None, ext=None, pref=None, remove=True):
  '''Split file into a) n files or b) files of size b (MB)
  '''
  assert b is not None or n is not None
  if ext is None:
    ext = os.path.splitext(inFile)[1]
  if oDir is None:
    oDir = os.path.dirname(inFile)
  if pref is None:
    pref = getBaseName(inFile)

  if n:
    arg = f'-n l/{n}'
  elif b:
    arg = f'-C{b}m'

  subprocess.check_call(f'split {arg} {inFile} {pref}', shell=True, cwd=oDir)

  oFiles, i = [], 1
  for file in os.listdir(oDir):
    if pref in file and file != getBaseFileName(inFile):
      nBase = f'{pref}_{i}{ext}'
      file = os.path.join(oDir, file)
      oFile = os.path.join(oDir, nBase)
      os.rename(file, oFile)
      i += 1
      oFiles.append(oFile)

  if remove:
    os.remove(inFile)
  return oFiles

def mergeSDFs(sdfFiles, outFile=None, oDir=''):
  inTexts = []
  for inFile in sdfFiles:
    with open(inFile) as f:
      inTexts.append(f.read().strip())

  outFile = os.path.join(oDir, f'mergedFiles{os.path.splitext(inFile)[-1]}') if outFile is None else outFile
  with open(outFile, 'w') as fo:
    fo.write('\n'.join(inTexts))
  return outFile

def mergePDBs(fn1, fn2, fnOut, hetatm2=False):
  with open(fnOut, 'w') as f:
    with open(fn1) as f1:
      for line in f1:
        if line.strip() != 'END' and not line.startswith('CONECT'):
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
    if newColor not in colors:
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
  inExt = os.path.splitext(os.path.basename(otherFile))[1]
  if inExt not in ['.pdb', '.mol2', '.sdf', '.mol']:
    inExt, otherFile = 'pdb', otherFile.replace(inExt, '.pdb')

  args = ' -ipdbqt {} -o{} -O {}'.format(os.path.abspath(pdbqtFile), inExt[1:], otherFile)
  runOpenBabel(protocol=protocol, args=args, popen=True)
  return os.path.abspath(otherFile)

def convertToSdf(protocol, molFile, sdfFile=None, overWrite=False, addHydrogens=False):
  '''Convert molecule files to sdf using openbabel'''
  def writeAddHParamsFile(outDir, molFn):
    oFile = os.path.join(outDir, 'addHydrogensParams.txt')
    with open(oFile, 'w') as f:
      f.write(f"ligandFiles: {molFn}\n")

      f.write(f'outputDir: {outDir}\n')
      f.write(f'doHydrogens: True\n')
      f.write(f'doGasteiger: True\n')
    return oFile

  if molFile.endswith('.sdf'):
    if sdfFile:
      shutil.copy(molFile, sdfFile)
      molFile = sdfFile
    return molFile
  if not sdfFile:
    baseName = getBaseName(molFile)
    outDir = os.path.abspath(protocol._getTmpPath())
    sdfFile = os.path.abspath(os.path.join(outDir, baseName + '.sdf'))
  else:
    baseName = getBaseName(sdfFile)
    outDir = os.path.abspath(os.path.dirname(sdfFile))

  if not os.path.exists(sdfFile) or overWrite:
    args = f' -i "{os.path.abspath(molFile)}" -of sdf --outputDir "{os.path.abspath(outDir)}" ' \
           f'--outputName {baseName} --overWrite'
    pwchemPlugin.runScript(protocol, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir, popen=True)

  if addHydrogens:
    paramFile = writeAddHParamsFile(outDir, sdfFile)
    pwchemPlugin.runScript(protocol, 'rdkit_addHydrogens.py', paramFile, env=RDKIT_DIC, cwd=outDir)

  return sdfFile

def getMAEMoleculeFiles(molList):
  '''Return in different lists the mae and non-mae files'''
  maeMols, otherMols = [], []
  for mol in molList:
    molFile = os.path.abspath(mol.getPoseFile())
    if '.mae' in molFile:
      maeMols.append(mol.clone())
    else:
      otherMols.append(mol.clone())
  return maeMols, otherMols

def runOpenBabel(protocol, args, cwd='/tmp', popen=False):
  pwchemPlugin.runOPENBABEL(protocol=protocol, args=args, cwd=cwd, popen=popen)

def splitConformerFile(confFile, outDir):
  writtenFiles = []
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
          writtenFiles.append(newFile)
          writeFile(towrite, newFile)
          towrite, lastRemark = line, True
          iConf += 1
      else:
        towrite += line
        lastRemark = False
  newFile = os.path.join(outDir, '{}-{}{}'.format(fnRoot, iConf, ext))
  writeFile(towrite, newFile)
  writtenFiles.append(newFile)
  return writtenFiles

def appendToConformersFile(confFile, newFile, outConfFile=None, beginning=True):
  '''Appends a molecule to a conformers file.
    If outConfFile == None, the output conformers file path is the same as as the start'''
  rename = False
  if outConfFile == None:
    iExt = os.path.splitext(confFile)[1]
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
def generate_gpf(protFile, spacing, xc, yc, zc, npts, outDir, ligandFns=None, znFFfile=None, addLigTypes=True):
  """
    Build the GPF file that is needed for AUTOGRID to generate the electrostatic grid
    """
  protName = os.path.splitext(os.path.basename(protFile))[0]
  gpfFile = os.path.join(outDir, protName + '.gpf')
  npts = int(round(npts))

  protAtomTypes = parseAtomTypes(protFile)

  if ligandFns == None or not addLigTypes:
      ligAtomTypes = 'A C HD N NA OA SA'
  else:
      ligAtomTypes = set([])
      for ligFn in ligandFns:
          ligAtomTypes = ligAtomTypes | parseAtomTypes(ligFn)

      ligAtomTypes = protAtomTypes.union(ligAtomTypes)
      ligAtomTypes = ' '.join(sortSet(ligAtomTypes))

  protAtomTypes = ' '.join(sortSet(protAtomTypes))

  with open(os.path.abspath(gpfFile), "w") as file:
    file.write("npts %s %s %s                        # num.grid points in xyz\n" % (npts, npts, npts))
    if znFFfile:
        file.write("parameter_file %s                        # force field default parameter file\n" % (znFFfile))
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
    if znFFfile:
        file.write('''nbp_r_eps 0.25 23.2135 12 6 NA TZ\nnbp_r_eps 2.1   3.8453 12 6 OA Zn\nnbp_r_eps 2.25  7.5914 12 6 SA Zn\nnbp_r_eps 1.0   0.0    12 6 HD Zn\nnbp_r_eps 2.0   0.0060 12 6 NA Zn\nnbp_r_eps 2.0   0.2966 12 6  N Zn''')

  return os.path.abspath(gpfFile)


def calculate_centerMass(atomStructFile):
  """
    Returns the geometric center of mass of an Entity (anything with a get_atoms function in biopython).
    Geometric assumes all masses are equal
    """

  try:
    structureHandler = AtomicStructHandler()
    structureHandler.read(atomStructFile)
    centerCoord = structureHandler.centerOfMass(geometric=True)
    structure = structureHandler.getStructure()

    return structure, centerCoord[0], centerCoord[1], centerCoord[2]  # structure, x,y,z

  except Exception as e:
    print("ERROR: ", "A pdb file was not entered in the Atomic structure field. Please enter it.", e)
    return

def parseAtomTypes(pdbqtFile, allowed=None, ignore=['Si', 'B', 'G0', 'CG0', 'G1', 'CG1']):
  atomTypes = set([])
  if pdbqtFile.endswith('.pdbqt'):
    with open(pdbqtFile) as f:
      for line in f:
        if line.startswith('ATOM') or line.startswith('HETATM'):
          pLine = line.split()
          at = pLine[-1]
          if (allowed is None or at in allowed) and at not in ignore:
              atomTypes.add(at)
  else:
      struct = PDBParser().get_structure("SASAstruct", pdbqtFile)

      for atom in struct.get_atoms():
        atomId = atom.get_id()
        if removeNumberFromStr(atomId) not in ignore:
          atomTypes.add(removeNumberFromStr(atomId))

  return atomTypes

def sortSet(seti):
  seti = list(seti)
  seti.sort()
  return seti

class CleanStructureSelect(Select):
  def __init__(self, chainIds, remHETATM, remWATER, het2keep=[], het2rem=[]):
    self.chain_ids = chainIds
    self.remHETATM, self.remWATER = remHETATM, remWATER
    self.het2keep, self.het2rem = het2keep, het2rem

  def accept_chain(self, chain):
    return not self.chain_ids or chain.id in self.chain_ids

  def accept_residue(self, residue):
    """ Recognition of heteroatoms - Remove water molecules """
    accept = True
    if self.remWATER and residue.id[0] == 'W':
      accept = False
    elif isHet(residue):
      if (self.remHETATM and residue.resname not in self.het2keep) or \
              (not self.remHETATM and residue.resname in self.het2rem):
        accept = False
    return accept

def cleanPDB(structFile, outFn, waters=False, hetatm=False, chainIds=None, het2keep=[], het2rem=[]):
  """ Extraction of the heteroatoms of .pdb files """
  structName = getBaseName(structFile)
  if os.path.splitext(structFile)[1] in ['.pdb', '.ent', '.pdbqt']:
    struct = PDBParser().get_structure(structName, structFile)
  elif structFile.endswith('.cif'):
    struct = MMCIFParser().get_structure(structName, structFile)
  else:
    print('Unknown format for file ', structFile)
    exit()

  io = PDBIO() if outFn.endswith('pdb') else MMCIFIO()
  io.set_structure(struct)
  io.save(outFn, CleanStructureSelect(chainIds, hetatm, waters, het2keep, het2rem))
  return outFn

def removeStartWithLines(inFile, start):
  '''Remove the lines in a file starting with a certain string'''
  subprocess.check_call(f"sed -i '/^{start}/d' {inFile} ", shell=True)

def obabelMolConversion(mol, outFormat, outDir, pose=False):
  '''Converts a molecule into the specified format using OBabel binary'''
  outFormat = outFormat[1:] if outFormat.startswith('.') else outFormat
  molFile = os.path.abspath(mol.getFileName()) if not pose else os.path.abspath(mol.getPoseFile())
  inName, inExt = os.path.splitext(os.path.basename(molFile))
  oFile = os.path.abspath(os.path.join(outDir, f'{inName}.{outFormat}'))

  if not molFile.endswith(f'.{outFormat}'):
    args = f' -i{inExt[1:]} {os.path.abspath(molFile)} -o{outFormat} -O {oFile}'
    runOpenBabel(protocol=None, args=args, cwd=outDir, popen=True)
    if outFormat == 'mol2':
      relabelMapAtomsMol2(oFile)
  else:
    os.link(molFile, oFile)
  return oFile

def relabelAtomsMol2(atomFile, i=''):
  '''Relabel the atom names so each atom type goes from 1 to x, so if there is only one oxygen named O7,
  it will be renamed to O1'''
  atomCount = {}
  auxFile = os.path.join(os.path.dirname(atomFile), f'{getBaseName(atomFile)}_aux{i}.mol2')
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

            atomCount[atomType] = 1 if atomType not in atomCount else atomCount[atomType] + 1
            atomName = atomName if atomType != atomName else atomName + str(atomCount[atomType])

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
  return newS

def getNumberFromStr(s):
  num = ''
  for i in s:
    if i.isdigit():
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
  '''Sort a list of strs containing numbers, taking into account that number real value
  [A3, A1, B2, B4] -> [A1, A3, B2, B4]
  '''
  convert = lambda text: int(text) if text.isdigit() else text.lower()
  alphanumKey = lambda key: [convert(c) for c in re.split('(\d+)', key)]
  return sorted(listi, key=alphanumKey, reverse=rev)

def numberSort(strings, rev=False):
  '''Sort a list of strs containing numbers, taking into account only that number real value,
  ignoring the letter part
  [A3, A1, B2, B4] -> [A1, B2, A3, B4]
  '''
  def keyFunct(string):
    # Use regular expression to extract the numeric part of the string
    match = re.search(r'\d+', string)
    if match:
      return int(match.group())
    else:
      return float('inf')  # If there are no numbers, put at the end

  return sorted(strings, key=keyFunct, reverse=rev)

def unifyAttributes(itemList):
  '''Unify the attributes of the objects in a list by setting the missing ones with None so all items
  end up having all the attributes'''
  attributes = getAllItemAttributes(itemList)
  for item in itemList:
    for attr in attributes:
      if not hasattr(item, attr):
        item.__setattr__(attr, attributes[attr])
  return itemList

def getAllItemAttributes(itemList):
  attributes = {}
  for item in itemList:
    newAttrs = getItemAttributes(item)
    attributes.update(newAttrs)
  return attributes

def getItemAttributes(item):
  '''Return a dic with the attributes of an object and its values set to None in the specified type'''
  attributes = {}
  attrKeys = item.getObjDict().keys()
  for attrK in attrKeys:
    if attrK not in attributes and hasattr(item, attrK):
      value = item.__getattribute__(attrK)
      attributes[attrK] = value.clone()
      attributes[attrK].set(None)
  return attributes

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
    if isinstance(inpSet, Pointer):
      inpSet = inpSet.get()
    item = inpSet.getFirstItem()
    attributes.update(getItemAttributes(item))
  return attributes

def addToDic(dic, key, item):
  '''Add an element to a dic list creting it if the key was not there'''
  if key in dic:
    dic[key].append(item)
  else:
    dic[key] = [item]
  return dic

def clusterSurfaceCoords(surfCoords, intraDist):
  clusters = []
  for coord in surfCoords:
    newClusters = []
    newClust = [coord]
    for clust in clusters:
      merge = False
      for cCoord in clust:
        dist = calculateDistance(coord, cCoord)
        if dist < intraDist:
          merge = True
          break

      if merge:
        newClust += clust
      else:
        newClusters.append(clust)

    newClusters.append(newClust)
    clusters = newClusters.copy()
  return clusters

def createPocketFile(clust, i, outDir):
  outFile = os.path.join(outDir, f'pocketFile_{i}.pdb')
  with open(outFile, 'w') as f:
    for j, coord in enumerate(clust):
      f.write(writePDBLine(['HETATM', str(j), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
  return outFile

################# Wizard utils #####################
def getChainIds(chainStr):
  '''Parses a line of json with the description of a chain or chains and returns the ids'''
  chainJson = json.loads(chainStr)
  if 'chain' in chainJson:
    chainIds = [chainJson["chain"].upper().strip()]
  elif 'model-chain' in chainJson:
    modelChains = chainJson["model-chain"].upper().strip()
    chainIds = [x.split('-')[1] for x in modelChains.split(',')]
  return chainIds

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

def runInParallel(func, *args, paramList, jobs):
  """
  This function creates a pool of workers to run the given function in parallel.
  Also returns a list with the failed commands.
  """
  # Create a pool of worker processes
  nJobs = len(paramList) if len(paramList) < jobs else jobs
  pool = multiprocessing.Pool(processes=nJobs)

  # Apply the given function to the given param list using the pool
  results = [pool.apply_async(func, args=(param, *args,)) for param in paramList]

  # Initializing list of failed commands
  failedCommands = []

  # Check if any process encountered an error
  for result in results:
    if result.get():
      failedCommands.append(result.get())

  # Close the pool to release resources
  pool.close()

  # Join the pool, waiting for all processes to complete
  pool.join()

  # Return list of failed commands
  return failedCommands

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

################# File management utils #####################
def removeElements(elements):
  """ This function removes all given files and directories. """
  # Removing selected elements
  for item in elements:
    if os.path.exists(item):
      if os.path.isdir(item):
        shutil.rmtree(item)
      else:
        os.remove(item)

def createMSJDic(protocol):
  msjDic = {}
  for pName in protocol.getStageParamsDic(type='Normal').keys():
    if hasattr(protocol, pName):
      msjDic[pName] = getattr(protocol, pName).get()
    else:
      print('Something is wrong with parameter ', pName)

  for pName in protocol.getStageParamsDic(type='Enum').keys():
    if hasattr(protocol, pName):
      msjDic[pName] = protocol.getEnumText(pName)
    else:
      print('Something is wrong with parameter ', pName)
  return msjDic

########### SEQUENCE INTERACTING MOL UTILS ########
def getFilteredOutput(inSeqs, filtSeqNames, filtMolNames, scThres):
  '''Filters the setofsequences (inSeqs) to return an array with the interacting molecules scores
  The array is formed only by filt(Seq/Mol)Names and over the scThres score threshold'''
  intDic = inSeqs.getInteractScoresDic()

  seqNames, molNames = inSeqs.getSequenceNames(), inSeqs.getInteractMolNames()
  seqNames, molNames = filterNames(seqNames, molNames, filtSeqNames, filtMolNames)

  intAr = formatInteractionsArray(intDic, seqNames, molNames)
  intAr, seqNames, molNames = filterScores(intAr, seqNames, molNames, scThres)
  return intAr, seqNames, molNames

def filterNames(seqNames, molNames, filtSeqNames, filtMolNames):
  if 'All' not in filtSeqNames:
    seqNames = [seqName for seqName in seqNames if seqName in filtSeqNames]

  if 'All' not in filtMolNames:
    molNames = [molName for molName in molNames if molName in filtMolNames]

  return seqNames, molNames

def filterScores(intAr, seqNames, molNames, scThres):
  ips, ims = [], []

  for ip, seqName in enumerate(seqNames):
    if any(intAr[ip, :] >= scThres):
      ips.append(ip)

  for im, molName in enumerate(molNames):
    if any(intAr[:, im] >= scThres):
      ims.append(im)

  if len(seqNames) != len(ips):
    seqNames = list(np.array(seqNames)[ips])
    intAr = intAr[ips, :]

  if len(molNames) != len(ims):
    molNames = list(np.array(molNames)[ims])
    intAr = intAr[:, ims]

  return intAr, seqNames, molNames

def formatInteractionsArray(intDic, seqNames, molNames):
  intAr = np.zeros((len(seqNames), len(molNames)))
  for i, seqName in enumerate(seqNames):
    for j, molName in enumerate(molNames):
      intAr[i, j] = intDic[seqName][molName]
  return intAr
  
def normalizeToRange(iterable, normRange=[0, 1]):
  maxIt, minIt = max(iterable), min(iterable)
  return [((normRange[1] - normRange[0]) * (i - minIt)) / (maxIt - minIt) + normRange[0] for i in iterable]

def replaceInFiles(directory, old_string, new_string, file_extension="*"):
  command = f"find {directory} -type f -name '*{file_extension}' -exec sed -i 's/{old_string}/{new_string}/g' {{}} +"
  subprocess.run(command, shell=True, check=True)

def replaceInFile(file, inStr, repStr):
  inStr, repStr = inStr.replace('\n', '\\n'), repStr.replace('\n', '\\n')
  subprocess.check_call(f'''sed -i -z 's/{inStr}/{repStr}/g' {file}''', shell=True)
  return file

def getReplaceCommand(file, inStr, repStr):
  inStr, repStr = inStr.replace('\n', '\\n'), repStr.replace('\n', '\\n')
  quote = "'" if '"' in inStr or '"' in repStr else '"'
  return f'''sed -i -z {quote}s/{inStr}/{repStr}/g{quote} {file}'''

def flipDic(dic):
  return {v:k for k, v in dic.items()}
