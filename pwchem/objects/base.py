# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

import enum, io, subprocess, pickle

import pyworkflow.object as pwobj
import pwem.objects.data as data

from scipy import spatial
from pwchem.utils import *
from pwchem.constants import *


class DatabaseID(data.EMObject):
  """ Database identifier """

  def __init__(self, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.database = pwobj.String(kwargs.get('database', None))
    self.dbId = pwobj.String(kwargs.get('dbId', None))

  def getDbId(self):
    return self.dbId.get()

  def setDbId(self, value):
    self.dbId.set(value)

  def getDatabase(self):
    return self.database.get()

  def setDatabase(self, value):
    """Valid databases: pdb, uniprot, ... """
    self.database.set(value)

  def copyInfo(self, other, copyId=False):
    self.copyAttributes(other, 'database', 'dbId')
    if copyId:
      self.copyObjId(other)


class SetOfDatabaseID(data.EMSet):
  """ Set of DatabaseIDs """
  ITEM_TYPE = DatabaseID
  FILE_TEMPLATE_NAME = 'setOfDatabaseIds%s.sqlite'

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)


class SequenceChem(data.Sequence):
  def __init__(self, **kwargs):
    data.Sequence.__init__(self, **kwargs)
    self._attrFile = pwobj.String(kwargs.get('attributesFile', None))

    # Dictionary that contains the interacting score of each SequenceChem with each SmallMolecule,
    # so there is no need to create 1 SetOfSmallMolecules per SequenceChem {molId: score}
    self._interactScoresFile = pwobj.String(kwargs.get('interactScoreFile', None))

  def __str__(self):
    return "SequenceChem (name = {})\n".format(self.getSeqName())

  def getAttrFile(self):
    return self._attrFile.get()

  def setAttrFile(self, value):
    self._attrFile.set(value)

  def addAttributes(self, attrDic):
    attrFile = self.getAttrFile()
    with open(attrFile, 'a') as f:
      for key, values in attrDic.items():
        f.write('{}: {}\n'.format(key, values))

  def getAttributesDic(self):
    attrDic = {}
    with open(self.getAttrFile()) as f:
      for line in f:
        key, values = line.split(':')
        attrDic[key.strip()] = eval(values.strip())
    return attrDic

  def getInteractScoresDic(self):
    '''Returns a dictionary of the form {molName: score}, from the file where the interaction scores are stored.
        '''
    with open(self.getInteractScoresFile(), 'rb') as f:
      intDic = pickle.load(f)
    return intDic

  def setInteractScoresDic(self, intDic, outFile):
    '''From a dictionary of the form {molName: score}, writes a file in outFile
        containing that information to be stored in the object'''
    with open(outFile, 'wb') as f:
      pickle.dump(intDic, f)
    self.setInteractScoresFile(outFile)

  def getInteractScoresFile(self):
    return self._interactScoresFile.get()

  def setInteractScoresFile(self, intFile):
    self._interactScoresFile.set(intFile)


class SetOfSequencesChem(data.SetOfSequences):
  def __init__(self, **kwargs):
    data.SetOfSequences.__init__(self, **kwargs)
    self._aligned = pwobj.Boolean(kwargs.get('aligned', False))
    self._alignFile = pwobj.String(kwargs.get('alignFile', None))

    self._interactMols = pwobj.Pointer()
    self._interactScoresFile = pwobj.String(kwargs.get('interactScoreFile', None))

  def copyInfo(self, other):
    """ Copy basic information from other set of classes to current one"""
    self.copyAttributes(other, '_aligned', '_alignFile', '_interactMols', '_interactScoresFile')

  def getSetPath(self):
    return os.path.abspath(self._mapperPath[0])

  def getSetDir(self):
    return os.path.dirname(self.getSetPath())

  def getAligned(self):
    return self._aligned.get()

  def setAligned(self, value):
    self._aligned.set(value)

  def getAlignmentFileName(self):
    return self._alignFile.get()

  def setAlignmentFileName(self, value):
    return self._alignFile.set(value)

  def convertEMBOSSformat(self, embossFormat, outputFile):
    from pwchem import Plugin
    clRun = '%s && ' % (Plugin.getEnvActivationCommand(BIOCONDA_DIC))
    clRun += 'seqret -sequence {} -osformat2 {} {}'. \
      format(self.getAlignmentFileName(), embossFormat, outputFile)

    subprocess.check_call(clRun, shell=True)
    return outputFile

  def __str__(self):
    alignStr = super().__str__()
    alignStr += ', aligned={}'.format(self._aligned.get())
    return alignStr

  def getInteractMols(self):
    return self._interactMols.get()

  def getInteractMolsPointer(self):
    return self._interactMols

  def hasInteractMols(self):
    return self.getInteractMolsPointer() != pwobj.Pointer()

  def setInteractMols(self, mols=None):
    if mols.isPointer():
      self._interactMols.copy(mols)
    else:
      self._interactMols.set(mols)

  def getInteractScoresDic(self, calculate=False):
    '''Returns a dictionary of the form {seqName: {molName: score}},
        from the files where the interaction scores are stored.
        '''
    if not calculate and self.getInteractScoresFile() and os.path.getsize(self.getInteractScoresFile()) > 0:
      with open(self.getInteractScoresFile(), 'rb') as f:
        intDic = pickle.load(f)
    else:
      intDic = {}
      for seq in self:
        with open(seq.getInteractScoresFile(), 'rb') as f:
          intDic[seq.getSeqName()] = pickle.load(f)
    return intDic

  def setInteractScoresDic(self, intDic=None, outFile=None):
    '''From a dictionary of the form {seqName: {molName: score}}, writes a file in outFile
        containing that information to be stored in the object.
        If intDic=None, it is build from the elements of the set
        If outFile=None, the file is saved in the same directory as the set mapper'''
    if not intDic:
      intDic = self.getInteractScoresDic()
    if not outFile:
      outFile = os.path.join(self.getSetDir(), f'{super().__str__()}_interactions.pickle')
    with open(outFile, 'wb') as f:
      pickle.dump(intDic, f)
    self.setInteractScoresFile(outFile)

  def getInteractScoresFile(self):
    return self._interactScoresFile.get()

  def setInteractScoresFile(self, intFile):
    self._interactScoresFile.set(intFile)

  def getSequenceNames(self):
    return [seq.getSeqName() for seq in self]

  def getInteractMolNames(self):
    intDic = self.getInteractScoresDic()
    return list(set([molName for molDic in intDic.values() for molName in molDic]))

  def getInteractMolsNumber(self):
    n = 0
    if self.hasInteractMols():
      mols = self.getInteractMols()
      if isinstance(mols, SetOfSmallMolecules):
        n = len(mols)
      else:
        n = mols.getLength()
    return n


class SequenceVariants(data.EMFile):
  """A fasta file for a Protein ID"""

  def __init__(self, filename=None, **kwargs):
    data.EMFile.__init__(self, filename, **kwargs)
    self._sequence = None

  def __str__(self):
    return ("{} (id={}, numberOfVariants={})".format(self.getClassName(), self.getId(),
                                                     len(self.getMutationsInLineage())))

  def setSequence(self, sequenceObj):
    self._sequence = sequenceObj

  def getId(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getId()

  def getSeqName(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getSeqName()

  def setVariantsFileName(self, fnVars):
    self.setFileName(fnVars)

  def getVariantsFileName(self):
    return self.getFileName()

  def getSequence(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getSequence()

  def getMutationsInLineage(self):
    fnVars = self.getVariantsFileName()
    var2mutDic = {}
    with open(fnVars) as f:
      for line in f:
        lineInfo = line.split(';')[0].strip()
        mut, varGroups = lineInfo.split()[0], lineInfo.split()[1:]
        for vGroup in varGroups:
          if vGroup.endswith(',') or vGroup.endswith('.'):
            vGroup = vGroup[:-1]
          for variant in vGroup.split('/'):
            if variant != 'and':
              if variant in var2mutDic:
                var2mutDic[variant] += [mut]
              else:
                var2mutDic[variant] = [mut]

    return var2mutDic

  def generateVariantLineage(self, selectedVariant):
    '''Perform the substitutions related to a variant'''
    varDic = self.getMutationsInLineage()
    mutList = varDic[selectedVariant]
    return self.performSubstitutions(mutList)

  def performSubstitutions(self, mutList):
    '''Perform substitutions in the sequence from a list as: ["L5F", "M89C", ...]'''
    sequence = self.getSequence()
    subsDic = {}
    for mutation in mutList:
      mutation = mutation.strip()
      mutIdx = getNumberFromStr(mutation)
      subsDic[int(mutIdx)] = mutation.split(mutIdx)[1]

    lettersSequence = list(sequence)
    for mutPos in subsDic:
      lettersSequence[mutPos - 1] = subsDic[mutPos]

    mutatedSequence = ''.join(lettersSequence)
    return mutatedSequence

  def exportToFile(self, outPath):
    if os.path.exists(outPath):
      os.remove(outPath)
    wholeSeq = self.getSequence()
    self._sequence.exportToFile(outPath)

    var2Mut = self.getMutationsInLineage()
    for var in var2Mut:
      tmpSeq = ['-'] * len(wholeSeq)
      for mut in var2Mut[var]:
        mut = mut.strip()
        mutIdx = getNumberFromStr(mut)
        tmpSeq[int(mutIdx) - 1] = mut.split(mutIdx)[1]
      tmpSeq = ''.join(tmpSeq)

      tmpSeqObj = Sequence(sequence=tmpSeq, id=var)
      tmpSeqObj.appendToFile(outPath, doClean=False)


class SmallMolecule(data.EMObject):
  """ Small molecule """

  def __init__(self, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.smallMoleculeFile = pwobj.String(kwargs.get('smallMolFilename', None))
    self.poseFile = pwobj.String(kwargs.get('poseFile', None))  # File of position
    self._mappingFile = pwobj.String(kwargs.get('mappingFile', None))
    self.confId = pwobj.Integer(kwargs.get('confId', None))  # pocketID
    self.molName = pwobj.String(kwargs.get('molName', None))
    if self.molName.get() and self.molName.get().lower() == 'guess':
      self.molName.set(self.guessMolName())

    self.gridId = pwobj.Integer(kwargs.get('gridId', None))  # pocketID
    self.poseId = pwobj.Integer(kwargs.get('poseId', None))
    self.dockId = pwobj.Integer(kwargs.get('dockId', None))  # dockProtocol ID
    self._type = pwobj.String(kwargs.get('type', 'Standard'))

    self.proteinFile = pwobj.String(kwargs.get('proteinFile', None))  # to be used when each mol has diff receptor

  def __str__(self):
    s = '{} ({} molecule)'.format(self.getClassName(), self.getUniqueName())
    return s

  def getFileName(self):
    '''Original filename of the molecule prior to any docking'''
    return self.smallMoleculeFile.get()

  def setFileName(self, value):
    self.smallMoleculeFile.set(value)

  def getMolName(self):
    molName = self.molName.get()
    if not molName:
      molName = self.guessMolName()
    return molName

  def setMolName(self, value):
    self.molName.set(value)

  def guessMolName(self):
    fBase = self.getFileName().split('/')[-1].split('.')[0]
    if self.getConfId():
      return '-'.join(fBase.split('-')[:-1])
    else:
      return fBase

  def getMolBase(self):
    return self.getMolName()

  def getPoseFile(self):
    '''Filename of the molecule after docking'''
    return self.poseFile.get()

  def setPoseFile(self, value):
    return self.poseFile.set(value)

  def getMappingFile(self):
    '''Filename of the molecule after docking'''
    return self._mappingFile.get()

  def setMappingFile(self, value):
    return self._mappingFile.set(value)

  def setPoseId(self, value):
    return self.poseId.set(value)

  def getPoseId(self):
    return self.poseId.get()

  def getConfId(self):
    return self.confId.get()

  def setConfId(self, confId):
    self.confId.set(confId)

  def getGridId(self):
    return self.gridId.get()

  def setGridId(self, gridId):
    self.gridId.set(gridId)

  def getDockId(self):
    return self.dockId.get()

  def setDockId(self, value):
    self.dockId.set(value)

  def getEnergy(self):
    if hasattr(self, '_energy'):
      return self._energy.get()

  def setEnergy(self, value):
    if hasattr(self, '_energy'):
      self._energy.set(value)
    else:
      self._energy = pwobj.Float(value)

  def getConformersFileName(self):
    if hasattr(self, '_ConformersFile'):
      return self._ConformersFile.get()

  def getParamsFileName(self):
    if hasattr(self, '_ParamsFile'):
      return self._ParamsFile.get()

  def getPDBFileName(self):
    if hasattr(self, '_PDBFile'):
      return self._PDBFile.get()

  def setProteinFile(self, value):
    self.proteinFile.set(value)

  def getProteinFile(self):
    return self.proteinFile.get()

  def getMolClass(self):
    return self._type

  def setMolClass(self, value):
    self._type.set(value)

  def getUniqueName(self, grid=True, conf=True, pose=True, dock=True):
    name = self.getMolName()
    if self.getGridId() and grid:
      name = 'g{}_'.format(self.getGridId()) + name
    if self.getConfId() and conf:
      name += '-{}'.format(self.getConfId())
    if self.getPoseId() and pose:
      name += '_{}'.format(self.getPoseId())
    if self.getDockId() and dock:
      name += '_{}'.format(self.getDockId())
    return name

  def getAtomsPosDic(self, onlyHeavy=True):
    '''Returns a dictionary with the atoms coordinates:
        {atomId: [c1, c2, c3], ...}'''
    molFile = self.getFileName()
    if self.getPoseFile() != None:
      molFile = self.getPoseFile()

    posDic = {}
    if '.pdb' in molFile:
      with open(molFile) as fIn:
        for line in fIn:
          if line.startswith('ATOM') or line.startswith('HETATM'):
            elements = splitPDBLine(line)
            atomId, coords = elements[2], elements[6:9]
            atomType = removeNumberFromStr(atomId)
            if atomType != 'H' or not onlyHeavy:
              posDic[atomId] = list(map(float, coords))

    elif molFile.endswith('.mol2'):
      with open(molFile) as fIn:
        parse = False
        for line in fIn:
          if parse and line.startswith('@'):
            # finished
            return posDic

          if parse:
            elements = line.split()
            atomId, atomType, coords = elements[1], elements[5], elements[2:5]
            if atomType != 'H' or not onlyHeavy:
              posDic[atomId] = list(map(float, coords))

          elif line.startswith('@<TRIPOS>ATOM'):
            parse = True

    elif molFile.endswith('.sdf'):
      with io.open(molFile, encoding='cp1252') as fIn:
        numType = {}
        for i, line in enumerate(fIn):
          if i >= 4:
            elements = line.split()
            if len(elements) > 7:
              atomType, coords = elements[3], elements[:3]
              if atomType in numType:
                numType[atomType] += 1
              else:
                numType[atomType] = 1

              atomId = '{}{}'.format(atomType, numType[atomType])
              if atomType != 'H' or not onlyHeavy:
                posDic[atomId] = list(map(float, coords))
            else:
              break

    return posDic

  def getMapDic(self):
    mapDic = {}
    if self.getMappingFile():
      with open(self.getMappingFile()) as fIn:
        for line in fIn:
          sline = line.split()
          mapDic[sline[1]] = sline[0]
    return mapDic

  def mapLabels(self, probMol, mapBy='coords', distThres=0.1):
    '''Maps the labels of the reference molecule to other probe molecule, which must be the same and
        have the same positions. Returns a dic with the mapping'''
    mapDic = {}
    refCoords, probCoords = self.getAtomsPosDic(), probMol.getAtomsPosDic()
    if mapBy == 'coords':
      for refLabel in refCoords:
        for probLabel in probCoords:
          dist = calculateDistance(refCoords[refLabel], probCoords[probLabel])
          if dist <= distThres:
            mapDic[refLabel.upper()] = probLabel.upper()

    elif mapBy == 'order':
      for refLabel, probLabel in zip(refCoords, probCoords):
        mapDic[refLabel.upper()] = probLabel.upper()

    if len(mapDic) == len(refCoords):
      return mapDic
    print('Atom label mapping could not be completed correctly')

  def writeMapFile(self, probMol, mapBy='coords', overWrite=False, outFile=None, outDir=None):
    if not outFile:
      if not outDir:
        outDir = os.path.dirname(probMol.getFileName())
      outFile = os.path.join(outDir, 'mapping_{}.tsv'.format(self.getMolName()))

    if not os.path.exists(outFile) or overWrite:
      mapDic = self.mapLabels(probMol, mapBy)
      if mapDic:
        with open(outFile, 'w') as f:
          for refLabel in mapDic:
            f.write('{}\t{}\n'.format(refLabel, mapDic[refLabel]))
    return outFile


class SetOfSmallMolecules(data.EMSet):
  """ Set of Small molecules """
  ITEM_TYPE = SmallMolecule
  FILE_TEMPLATE_NAME = 'setOfSmallMolecules%s.sqlite'

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)
    self._molClass = pwobj.String('Standard')
    self.proteinFile = pwobj.String(kwargs.get('proteinFile', None))
    self._docked = pwobj.Boolean(False)

  def __str__(self):
    s = '{} ({} items, {} class)'.format(self.getClassName(), self.getSize(), self.getMolClass())
    return s

  def clone(self):
    clone = self.getClass()()
    clone.copy(self, ignoreAttrs=[])
    return clone

  def copyInfo(self, other):
    self._molClass = other._molClass
    self.proteinFile = other.proteinFile
    self._docked = other._docked

  def getMolClass(self):
    return self._molClass.get()

  def updateMolClass(self):
    mClass = self.getMolClass()
    for mol in self:
      if mol.getMolClass() != mClass:
        if mClass == 'Standard':
          mClass = mol.getMolClass()
        else:
          mClass = 'Mixed'
          break
    self._molClass.set(mClass)

  def append(self, item, update=False):
    super().append(item)
    if update:
      self.updateMolClass()

  def getSetPath(self):
    return os.path.abspath(self._mapperPath[0])

  def getSetDir(self):
    return '/'.join(self.getSetPath().split('/')[:-1])

  def setProteinFile(self, value):
    self.proteinFile.set(value)

  def getProteinFile(self):
    return self.proteinFile.get()

  def getOriginalReceptorFile(self):
    return self.proteinFile.get()

  def isDocked(self):
    return self._docked.get()

  def setDocked(self, value=True):
    self._docked.set(value)

  def getUniqueMolNames(self):
    names = []
    for mol in self:
      names.append(mol.getUniqueName())
    return names

  def updateLigandsDic(self, outputLigandsDic, molSet, vType, oLabelSet='outputSmallMolecules'):
    '''Generates and updates a dictionary containing the indexes of the objects for each of the subsets found
        in the set for the viewer: based on set, molecule, pocket and single'''
    SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'
    if vType == SET:
      oLabel = oLabelSet
      outputLigandsDic[oLabel] = [mol.getObjId() for mol in molSet]
    else:
      for mol in molSet:
        curMol = mol.clone()
        if vType == SINGLE:
          oLabel = curMol.getUniqueName()
        elif vType == MOLECULE:
          oLabel = curMol.getMolName()
        elif vType == POCKET:
          oLabel = 'g_{}'.format(curMol.getGridId())

        if oLabel not in outputLigandsDic:
          outputLigandsDic[oLabel] = [curMol.getObjId()]
        else:
          outputLigandsDic[oLabel] += [curMol.getObjId()]
    return outputLigandsDic

  def getGroupIndexesPath(self, setLabel=None):
    setLabel = getBaseName(self.getSetPath()) if setLabel is None else setLabel
    return os.path.join(os.path.dirname(self.getSetPath()), f'{setLabel}_indexGroups.pickle')

  def saveGroupIndexes(self, setLabel=None, outputLigandsDic=None):
    '''Computes the subgroups dictionary and saves it into the corresponding file or saves it directly if passed'''
    setLabel = getBaseName(self.getSetPath()) if setLabel is None else setLabel
    outputLigandsDic = self.getGroupIndexes() if outputLigandsDic is None else outputLigandsDic
    indexFile = open(self.getGroupIndexesPath(setLabel), 'wb')
    pickle.dump(outputLigandsDic, indexFile)

  def getGroupIndexes(self, setLabel=None):
    '''Return the dictionary containing the subsets information as:
            {subGroupLabel: [subGroupIndexes]}
        Check if the pickle file containing it exists and loads it (if not compute) else, it generates it
        '''
    setLabel = getBaseName(self.getSetPath()) if setLabel is None else setLabel
    indexFile = self.getGroupIndexesPath(setLabel)
    if os.path.exists(indexFile):
      with open(indexFile, 'rb') as f:
        outputLigandsDic = pickle.load(f)
    else:
      outputLigandsDic = self.computeGroupIndexes(setLabel)
      self.saveGroupIndexes(setLabel, outputLigandsDic)
    return outputLigandsDic

  def computeGroupIndexes(self, setLabel=None):
    '''Generete the subgroup index dictionaries'''
    setLabel = getBaseName(self.getSetPath()) if setLabel is None else setLabel
    outputLigandsDic = {}
    for vType in ['single', 'molecule', 'pocket', 'set']:
      outputLigandsDic[vType] = self.updateLigandsDic({}, self, vType, oLabelSet=setLabel)
    return outputLigandsDic

  def getMolsFromIds(self, itemIds):
    '''Return a list of molecules '''
    mols, mapper = [], self._getMapper()
    for itemId in itemIds:
      mols.append(mapper.selectById(itemId).clone())
    self.close()
    return mols


class SmallMoleculesLibrary(data.EMObject):
  '''Small molecules library object to wrap a single small molecules file containing the whole library in SMI
    '''

  def __init__(self, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.libraryFile = pwobj.String(kwargs.get('libraryFilename', None))
    self.origin = pwobj.String(kwargs.get('origin', None))
    self.length = pwobj.Integer(kwargs.get('length', None))
    self.headers = pwobj.List(objName='_headers')
    inHeaders = kwargs.get('headers', ['SMI', 'Name'])
    self.setHeaders(inHeaders)

  def __str__(self):
    length = self.getLength()
    if not length:
      length = self.calculateLength()

    s = f'{self.getClassName()} ({length} items, headers: {self.getHeaders()})'
    return s

  def calculateLength(self):
    inFile = self.getFileName()
    count = sum(1 for _ in open(inFile, 'rb'))
    self.setLength(count)
    return count

  def getLength(self):
    return self.length.get()

  def setLength(self, value):
    self.length.set(value)

  def getFileName(self):
    return self.libraryFile.get()

  def setFileName(self, value):
    self.libraryFile.set(value)
    self.calculateLength()

  def getLibraryMap(self, inverted=False, fullLine=False, lineDic=False):
    '''Returns a map dictionary as: {smi: name} or {name: smi} if inverted
        '''
    mapDic = {}
    for key, val in self.yieldLibraryMapItems(inverted, fullLine, lineDic):
        mapDic[key] = val
    return mapDic

  def yieldLibraryMapItems(self, inverted=False, fullLine=False, lineDic=False):
    with open(self.getFileName()) as f:
      for line in f:
        smi, name = line.split()[0].strip(), line.split()[1].strip()
        key = name if inverted else smi
        if not fullLine:
          val = smi if inverted else name
        else:
          if lineDic:
            val = {ki: vi for ki, vi in zip(self.getHeaders(), line.strip().split())}
          else:
            val = line.strip()
        yield key, val

  def getHeaders(self):
    hs = [h.get() for h in self.headers]
    return hs

  def setHeaders(self, headers):
    self.headers.clear()
    for head in headers:
      self.headers.append(pwobj.String(head))

  def validateSplit(self):
    maxMols = pwchemPlugin.getVar(MAX_MOLS_SET)
    return self.length <= maxMols

  def splitInFiles(self, outDir, col=1):
    assert self.validateSplit(), WARNLIBBIG
    oFiles = []
    inFile = self.getFileName()
    with open(inFile) as f:
      for line in f:
        molName = line.split()[col]

        oFile = os.path.join(outDir, f'{molName}.smi')
        with open(oFile, 'w') as fO:
          fO.write(line)
        oFiles.append(oFile)
    return oFiles

  def yieldLibraryValues(self, colIndexes, batchSize=1024):
    batch = []
    with open(self.getFileName()) as f:
      for line in f:
        row = line.split()
        try:
          vals = [row[ci].strip() for ci in colIndexes]
          batch.append(vals)
        except (IndexError):
          continue

        if len(batch) >= batchSize:
          yield batch
          batch = []

      if batch:
        yield batch

class BindingSite(data.EMObject):
  """ Binding site """

  def __init__(self, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.bindingSiteFile = pwobj.String(kwargs.get('bindingSiteFilename', None))
    self.structurePtr = None

  def getFileName(self):
    return self.bindingSiteFile.get()


class SetOfBindingSites(data.EMSet):
  """ Set of Binding sites """
  ITEM_TYPE = BindingSite
  FILE_TEMPLATE_NAME = 'setOfBindingSites%s.sqlite'
  EXPOSE_ITEMS = True

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)


class SequenceROI(data.EMObject):
  """ Represent a sequence ROI"""

  def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self._sequence = kwargs.get('sequence', None)  # Sequence object
    self._ROISequence = kwargs.get('seqROI', None)  # Sequence object
    self._roiIdx = Integer(kwargs.get('roiIdx', None))
    self._roiIdx2 = Integer(kwargs.get('roiIdx2', None))
    self._type = String(kwargs.get('type', None))

  def __str__(self):
    s = '{} (Idx: {}-{}, ROI: {})'.format(self.getClassName(), self.getROIIdx(), self.getROIIdx2(),
                                          self.getROISequence())
    return s

  def setSequence(self, sequenceObj):
    self._sequence = sequenceObj

  def getSequence(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getSequence()

  def getId(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getId()

  def getSeqName(self):
    if type(self._sequence) == Sequence:
      return self._sequence.getSeqName()

  def getROIId(self):
    return self._ROISequence.getId()

  def setROIId(self, roiId):
    self._ROISequence.setId(roiId)

  def getROIName(self):
    return self._ROISequence.getSeqName()

  def setROIName(self, roiName):
    self._ROISequence.setSeqName(roiName)

  def getROISequence(self):
    return self._ROISequence.getSequence()

  def setROISequence(self, seqObj):
    self._ROISequence.set(seqObj)

  def findROIIdx(self):
    # python idx starts with 0, sequences with 1
    return self.getSequence().find(self.getROISequence()) + 1

  def getROIIdx(self):
    return self._roiIdx.get()

  def setROIIdx(self, idx):
    self._roiIdx.set(idx)

  def getROIIdx2(self):
    return self._roiIdx2.get()

  def setROIIdx2(self, idx):
    self._roiIdx2.set(idx)

  def setROIIdxs(self, idxs):
    self.setROIIdx(idxs[0])
    self.setROIIdx2(idxs[1])

  def getROIIdxs(self):
    return [self._roiIdx.get(), self._roiIdx2.get()]

  def getROILength(self):
    return len(self.getROISequence())

  def setType(self, type):
    self._type.set(type)

  def getType(self):
    return self._type.get()


class SetOfSequenceROIs(data.EMSet):
  ITEM_TYPE = SequenceROI

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)

  def __str__(self):
    s = '{} ({} items, sequence={})'.format(self.getClassName(), self.getSize(),
                                            self.getSequenceObj().getId())
    return s

  def getSetPath(self):
    return os.path.abspath(self._mapperPath[0])

  def getSetDir(self):
    return '/'.join(self.getSetPath().split('/')[:-1])

  def getSequenceObj(self):
    return self.getFirstItem()._sequence

  def getSequence(self):
    return self.getSequenceObj().getSequence()

  def exportToFile(self, outPath, mainSeq=True, minLen=None, maxLen=None):
    '''Export the whole set of sequence ROIs to a fasta file.
        If mainSeq, the whole sequence of the ROIs will first be written and the ROIs will be aligned to it.
        Else, only the ROIs, unaligned, will be written'''
    if os.path.exists(outPath):
      os.remove(outPath)
    if mainSeq:
      wholeSeqObj = self.getSequenceObj()
      wholeSeq = wholeSeqObj.getSequence()
      wholeSeqObj.exportToFile(outPath)
    for roi in self:
      roiSeq, roiIdx = roi.getROISequence(), roi.getROIIdx()
      if (not minLen or len(roiSeq) >= minLen) and (not maxLen or len(roiSeq) <= maxLen):
        if mainSeq:
          tmpSeq = ['-'] * len(wholeSeq)
          r1, r2 = roiIdx - 1, roiIdx - 1 + len(roiSeq)
          tmpSeq[r1:r2] = roiSeq
        else:
          tmpSeq = roiSeq

        tmpSeqObj = Sequence(sequence=''.join(tmpSeq), id=roi._ROISequence.getId())
        tmpSeqObj.appendToFile(outPath, doClean=False)


class MultiEpitope(SetOfSequenceROIs):
  ITEM_TYPE = SequenceROI

  def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self._sequence = kwargs.get('sequence', None)

  def __str__(self):
    s = '{} ({} epitopes)'.format(self.getClassName(), self.countEpitopes())
    return s

  def countEpitopes(self):
    return len(self.getEpitopes())

  def getEpitopes(self):
    eps = []
    for item in self:
      if hasattr(item, '_type') and getattr(item, '_type') == 'Epitope':
        eps.append(item)
    return eps

  def getLinkers(self):
    links = []
    for item in self:
      if hasattr(item, '_type') and getattr(item, '_type') == 'Linker':
        links.append(item)
    return links

  def getMultiEpitopeSequence(self):
    s = ''
    for item in self:
      s += item.getROISequence()
    return s


class StructROI(data.EMFile):
  """ Represent a structural region of interest"""

  def __init__(self, filename=None, proteinFile=None, extraFile=None, pClass='Standard', **kwargs):
    self._class = String(pClass)
    if filename != None and extraFile != None:
      self.properties, self.pocketId = self.parseFile(extraFile, filename)
      kwargs.update(self.getKwargs(self.properties, POCKET_ATTRIBUTES_MAPPING))

    data.EMFile.__init__(self, filename, **kwargs)
    if hasattr(self, 'pocketId') and self.pocketId:
      self.setObjId(self.pocketId)
    self._proteinFile = String(proteinFile)
    self._extraFile = String(extraFile)
    self._nPoints = Integer(kwargs.get('nPoints', None))
    self._contactResidues = String(kwargs.get('contactResidues', None))
    self._contactAtoms = String(kwargs.get('contactAtoms', None))
    self._volume = Float(kwargs.get('volume', None))
    self._score = Float(kwargs.get('score', None))
    self._energy = Float(kwargs.get('energy', None))

    if proteinFile != None:
      self.calculateContacts()

  def __str__(self):
    s = 'Structural ROI {}, {} class'.format(self.getObjId(), self.getPocketClass())
    return s

  # Attributes functions
  def getPocketClass(self):
    return str(self._class)

  def getNumberOfPoints(self):
    return self._nPoints

  def setNumberOfPoints(self, value):
    self._nPoints.set(value)

  def getContactResidues(self):
    return self._contactResidues

  def getDecodedCResidues(self):
    return self.decodeIds(self.getContactResidues())

  def setContactResidues(self, values):
    self._contactResidues.set(values)

  def getNumberOfContactResidues(self):
    return len(self.decodeIds(str(self._contactResidues)))

  def getContactAtoms(self):
    return self._contactAtoms

  def getDecodedCAtoms(self):
    return self.decodeIds(self.getContactAtoms())

  def setContactAtoms(self, values):
    self._contactAtoms.set(values)

  def getNumberOfContactAtoms(self):
    return len(self.decodeIds(str(self._contactAtoms)))

  def getVolume(self):
    return self._volume

  def setVolume(self, value):
    self._volume.set(value)

  def getScore(self):
    return self._score

  def setScore(self, value):
    self._score.set(value)

  def getEnergy(self):
    return self._score

  def setEnergy(self, value):
    self._score.set(value)

  def setProteinFile(self, value):
    self._proteinFile.set(value)

  def getProteinFile(self):
    return str(self._proteinFile)

  def getKwargs(self, props, am):
    nkwargs = {}
    for k in props:
      if k in am:
        nkwargs[am[k]] = props[k]
    return nkwargs

  def calculateContacts(self):
    cAtoms = self.buildContactAtoms(calculate=True)
    self.setContactAtoms(self.encodeIds(self.getAtomsIds(cAtoms)))
    cResidues = self.getResiduesFromAtoms(cAtoms)
    self.setContactResidues(self.encodeIds(self.getResiduesIds(cResidues)))

  # Complex pocket attributes functions
  def buildContactAtoms(self, calculate=False, maxDistance=5):
    '''Return the reported proteins atoms in contact with the structural ROI.
        If not reported, returns the protein atoms at less than 5A than any ROI point'''
    contactCodes = self.getContactAtoms()
    contactAtoms = []
    if str(contactCodes) != 'None' and not calculate:
      contactsIds = self.decodeIds(contactCodes)
      with open(self.getProteinFile()) as f:
        for line in f:
          if line.startswith('ATOM'):
            atomId = splitPDBLine(line)[1]
            if atomId in contactsIds:
              contactAtoms.append(ProteinAtom(line))
    else:
      # Manually calculate the contacts
      pocketCoords = self.getPointsCoords()
      proteinAtoms = self.getProteinAtoms()
      proteinCoords = self.getAtomsCoords(proteinAtoms)
      dists = spatial.distance.cdist(proteinCoords, pocketCoords)
      for i in range(len(proteinCoords)):
        if min(dists[i, :]) < maxDistance:
          contactAtoms.append(proteinAtoms[i])
    return contactAtoms

  def calculateMassCenter(self, weights=None):
    '''Calculates the center of mass of a set of points: [(x,y,z), (x,y,z),...]
        A weight for each point can be specified'''
    if weights == None:
      weights = self.getSpheresRadius()
      if weights == []:
        weights = None
    coords = self.getPointsCoords()
    return list(np.average(coords, axis=0, weights=weights))

  def getDiameter(self, radius=[]):
    '''Returning max distance of points found in the convex hull
        Possibility of adding radius to each row and column to calculate the diameter
        on spheres instead of points'''
    if radius == []:
      radius = self.getSpheresRadius()
    coords = np.array(self.getPointsCoords())
    if len(coords) > 3:
      cHullIndex = spatial.ConvexHull(coords).vertices
    else:
      cHullIndex = np.array(list(range(len(coords))))
    candidates = coords[cHullIndex]
    distMat = spatial.distance_matrix(candidates, candidates)
    if radius != []:
      distMat = self.addRadius(distMat, np.array(radius)[cHullIndex])
    i, j = np.unravel_index(distMat.argmax(), distMat.shape)

    return distMat[i, j]

  def addRadius(self, dMat, radius):
    '''Add the radius of each alpha sphere to their corresponding row and column in the distances
        matrix'''
    for i in range(dMat.shape[0]):
      dMat[i, :] += radius[i]
      dMat[:, i] += radius[i]
    return dMat

  def getPocketVolume(self, mode=0):
    '''Return the pocket volume:
        0: getSurfaceConvexVolume (over protein contact points)
        1: reported volume
        2: getConvexVolume (over pocket points)
        '''
    if mode == 2:
      return self.getConvexVolume()
    elif mode == 1 and self.getVolume() != None:
      return self.getVolume()
    else:
      vol = self.getSurfaceConvexVolume()
      if vol:
        return vol
      else:
        return self.getConvexVolume()

  def getConvexVolume(self):
    '''Calculate the convex volume of the points forming the pocket'''
    coords = np.array(self.getPointsCoords())
    cHull = spatial.ConvexHull(coords)
    return cHull.volume

  def getSurfaceConvexVolume(self):
    '''Calculate the convex volume of the protein contact atoms'''
    cAtoms = self.buildContactAtoms()
    cCoords = self.getAtomsCoords(cAtoms)
    if len(cCoords) >= 4:
      cHull = spatial.ConvexHull(cCoords)
      return cHull.volume
    else:
      return len(cCoords)

  def getSurfaceConvexArea(self):
    '''Calculate the convex area of the protein contact atoms'''
    cAtoms = self.buildContactAtoms()
    cCoords = self.getAtomsCoords(cAtoms)
    if len(cCoords) >= 3:
      cHull = spatial.ConvexHull(cCoords)
      return cHull.area
    else:
      return len(cCoords)

  def getPocketBox(self):
    '''Return the coordinates of the 2 corners determining the box (ortogonal to axis) where the pocket fits in
        For example: ([-1,0,2], [2,5,4])'''
    coords = np.array(self.getPointsCoords())
    return coords.min(axis=0), coords.max(axis=0)

  # Utils functions
  def encodeIds(self, idList):
    return '-'.join(idList)

  def decodeIds(self, idStr):
    return str(idStr).split('-')

  def getProteinAtoms(self):
    atoms = []
    with open(self.getProteinFile()) as f:
      for line in f:
        if line.startswith('ATOM'):
          atoms.append(ProteinAtom(line))
    return atoms

  def getProteinCoords(self):
    coords = []
    with open(self.getProteinFile()) as f:
      for line in f:
        if line.startswith('ATOM'):
          coords.append(tuple(line.split()[6:9]))
    return coords

  def getAtomsCoords(self, atoms):
    coords = []
    for atom in atoms:
      coords.append(atom.getCoords())
    return coords

  def getAtomsIds(self, atoms):
    ids = []
    for atom in atoms:
      ids.append(atom.atomId)
    ids = natural_sort(ids)
    return ids

  def getResiduesIds(self, residues):
    ids = set([])
    for res in residues:
      ids.add(res.residueId)
    ids = natural_sort(list(ids))
    return ids

  def getAtomsResidues(self, atoms):
    res = set([])
    for atom in atoms:
      res.add(atom.get)

  def buildPocketPoints(self):
    pocketPoints = []
    with open(self.getFileName()) as f:
      for line in f:
        if line.startswith('HETATM') or line.startswith('ATOM'):
          pocketPoints += [ProteinAtom(line)]
    return pocketPoints

  def getPointsCoords(self):
    coords = []
    for point in self.buildPocketPoints():
      coords.append(point.getCoords())
    return coords

  def getResiduesFromAtoms(self, atoms):
    res = []
    for atom in atoms:
      res.append(ProteinResidue(atom.line))
    return res

  def getMostCentralResidues(self, n=2):
    cMass = self.calculateMassCenter()
    cAtoms = self.buildContactAtoms()

    closestResidues = self.getCloserResidues(cMass, cAtoms, n)
    return closestResidues

  def getCloserResidues(self, refCoord, atoms, n=2):
    '''Returns the atoms sorted as they are close to the reference coordinate'''
    dists = []
    for at in atoms:
      dists += [calculateDistance(refCoord, at.getCoords())]

    zippedLists = sorted(zip(dists, atoms))
    dists, atoms = zip(*zippedLists)
    residues = self.getAtomResidues(atoms)

    closestResidues = []
    for r in residues:
      if r not in closestResidues:
        closestResidues.append(r)
        if len(closestResidues) == n:
          return closestResidues
    return closestResidues

  def getAtomResidues(self, atoms):
    residues = []
    for at in atoms:
      residues.append(at.residueId)
    return residues

  def parseFile(self, extraFile, filename):
    if self.getPocketClass() == 'FPocket':
      props, atoms, residues = {}, [], []
      atomsIds, residuesIds = [], []
      ini, parse = 'HEADER Information', False
      with open(extraFile) as f:
        for line in f:
          if line.startswith(ini):
            parse = True
            pocketId = int(line.split()[-1].replace(':', ''))
          elif line.startswith('HEADER') and parse:
            name = line.split('-')[1].split(':')[0]
            val = line.split(':')[-1]
            props[name.strip()] = float(val.strip())

          elif line.startswith('ATOM') and parse:
            atoms.append(ProteinAtom(line))
            atomsIds.append(atoms[-1].atomId)
            newResidue = ProteinResidue(line)
            if newResidue.residueId not in residuesIds:
              residues.append(newResidue)
              residuesIds.append(newResidue.residueId)
      props['contactAtoms'] = self.encodeIds(atomsIds)
      props['contactResidues'] = self.encodeIds(residuesIds)

    elif self.getPocketClass() == 'P2Rank':
      props = {}
      pocketId = int(filename.split('/')[-1].split('_')[1].split('.')[0])
      with open(extraFile) as f:
        keys = f.readline().split(',')
        for line in f:
          if int(line.split(',')[1]) == pocketId:
            values = line.split(',')

      for i, k in enumerate(keys):
        props[k.strip()] = values[i]

      props['residue_ids'] = props['residue_ids'].strip().replace(' ', '-')
      props['surf_atom_ids'] = props['surf_atom_ids'].strip().replace(' ', '-')
      pocketId = int(values[1])

    elif self.getPocketClass() == 'AutoLigand':
      props, i = {}, 1
      pocketId = int(filename.split('out')[1].split('.')[0])
      with open(extraFile) as f:
        for line in f:
          if i == pocketId:
            sline = line.split(',')
            # Volume
            props[sline[1].split('=')[0].strip()] = float(sline[1].split('=')[1].strip())
            # Energy/vol
            props[sline[2].strip()] = float(sline[3].split('=')[1].strip())
            if 'NumberOfClusters' in line:
              print(int(sline[4].split('=')[1].strip()))
              self.setNClusters(int(sline[4].split('=')[1].strip()))
          i += 1
      with open(filename) as f:
        npts, points = 0, []
        for line in f:
          line = line.split()
          points += [tuple(map(float, line[5:8]))]
          npts += 1
        props['nPoints'] = npts
        self.setAutoLigandPoints(points)

    elif self.getPocketClass() == 'SiteMap':
      props, pId = {}, 1
      pocketId = int(filename.split('-')[-1].split('.')[0])
      with open(extraFile) as fh:
        for line in fh:
          if line.startswith("SiteScore"):
            # SiteScore size   Dscore  volume  exposure enclosure contact  phobic   philic   balance  don/acc
            if pocketId == pId:
              keys = line.split()
              values = [float(x) for x in fh.readline().split()]
              props = dict(zip(keys, values))
            pId += 1

    else:
      props, pocketId = {}, None
    return props, pocketId

  def getStructureMaeFile(self):
    if hasattr(self, '_maeFile'):
      return self._maeFile.get()
    else:
      return None

  def setStructureMaeFile(self, value):
    if hasattr(self, '_maeFile'):
      self._maeFile.set(value)
    else:
      self._maeFile = String(value)

  def getAutoLigandPoints(self):
    if hasattr(self, '_adPoints'):
      return self._adPoints
    else:
      return None

  def setAutoLigandPoints(self, value):
    if hasattr(self, '_adPoints'):
      self._adPoints.set(value)
    else:
      self._adPoints = String(value)

  def incrNClusters(self):
    '''Increase in 1 the number of clusters which support a autoligand pocket'''
    self._nClusters.set(self.getNClusters() + 1)

  def getNClusters(self):
    if hasattr(self, '_nClusters'):
      return self._nClusters.get()
    else:
      return None

  def setNClusters(self, value):
    if hasattr(self, '_nClusters'):
      self._nClusters.set(value)
    else:
      self._nClusters = Integer(value)

  def getSpheresRadius(self):
    radius = []
    if self.getPocketClass() == 'FPocket':
      with open(str(self.getFileName())) as f:
        for line in f:
          if line.startswith('ATOM'):
            radius.append(float(line.split()[-1]))
    return radius


class SetOfStructROIs(data.EMSet):
  ITEM_TYPE = StructROI

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)
    self._pocketsClass = String(kwargs.get('pocketsClass', None))
    self._hetatmFile = String(kwargs.get('hetatmFile', None))

  def __str__(self):
    s = '{} ({} items, {} class)'.format(self.getClassName(), self.getSize(), self.getPocketsClass())
    return s

  def copyInfo(self, other):
    self._hetatmFile = other._hetatmFile
    self._pocketsClass = other._pocketsClass

  def getSetPath(self):
    return os.path.abspath(self._mapperPath[0])

  def getSetDir(self):
    return '/'.join(self.getSetPath().split('/')[:-1])

  def getProteinName(self):
    return self.getProteinFile().split('/')[-1].split('.')[0]

  def getProteinFile(self):
    return self.getFirstItem().getProteinFile()

  def getProteinSequencesDic(self):
      '''Returns the chains sequences for the protein file
      '''
      seqDic = {}
      from pwem.convert.atom_struct import AtomicStructHandler
      handler = AtomicStructHandler(self.getProteinFile())
      listChains, listRes = handler.getModelsChains()
      for chain in listChains[0]:
          seqDic[chain] = str(handler.getSequenceFromChain(0, chain))
      return seqDic

  def getProteinHetatmFile(self):
      return self._hetatmFile.get()

  def setProteinHetatmFile(self, value):
    self._hetatmFile.set(value)

  def getPocketsClass(self):
    return self._pocketsClass.get()

  def updatePocketsClass(self):
    pClass = self.getPocketsClass()
    for pocket in self:
      if pocket.getPocketClass() != pClass:
        if pClass == None:
          pClass = pocket.getPocketClass()
        else:
          pClass = 'Mixed'
          break
    self._pocketsClass.set(pClass)

  def append(self, item):
    super().append(item)
    self.updatePocketsClass()

  def buildPDBhetatmFile(self, suffix=''):
    protName = self.getProteinName()
    atmFile = self.getProteinFile()
    atmExt = os.path.splitext(atmFile)[1]
    outDir = self.getSetDir()
    outFile = os.path.join(outDir, protName + '{}_out{}'.format(suffix, atmExt))

    with open(outFile, 'w') as f:
      f.write(getRawPDBStr(atmFile, ter=False))
      f.write(self.getPocketsPDBStr())

    self.setProteinHetatmFile(outFile)
    return outFile

  def createPML(self, bBox=False):
    outHETMFile = os.path.abspath(self.getProteinHetatmFile())
    outExt = os.path.splitext(outHETMFile)[1]
    pmlFile = outHETMFile.replace('_out{}'.format(outExt), '.pml')

    # Creates the pml for pymol visualization
    idList = [pock.getObjId() for pock in self]
    with open(pmlFile, 'w') as f:
      if bBox:
        toWrite = FUNCTION_BOUNDING_BOX
        for pocket in self:
          pDia = pocket.getDiameter()
          toWrite += PML_BBOX_STR_EACH.format([0, 1, 0], pocket.calculateMassCenter(),
                                              [pDia * bBox] * 3,
                                              'BoundingBox_' + str(pocket.getObjId()))
        f.write(PML_BBOX_STR_POCK.format(outHETMFile, outHETMFile, idList, toWrite))
      else:
        f.write(PML_STR.format(outHETMFile, idList))

    return pmlFile

  def createSurfacePml(self, bBox=False):
    outHETMFile = os.path.abspath(self.getProteinHetatmFile())
    outExt = os.path.splitext(outHETMFile)[1]
    pmlFile = outHETMFile.replace('_out{}'.format(outExt), '_surf.pml')
    colors = createColorVectors(len(self))
    surfaceStr = ''
    for i, pock in enumerate(self):
      pId = pock.getObjId()
      surfAtomIds = str(list(map(int, pock.getDecodedCAtoms()))).replace(' ', '')
      surfaceStr += PML_SURF_EACH.format(pId, colors[i], pId, surfAtomIds, pId, pId)

    # Creates the pml for pymol visualization
    with open(pmlFile, 'w') as f:
      if bBox:
        toWrite = FUNCTION_BOUNDING_BOX
        for pocket in self:
          pDia = pocket.getDiameter()
          toWrite += PML_BBOX_STR_EACH.format([0, 1, 0], pocket.calculateMassCenter(),
                                              [pDia * bBox] * 3,
                                              'BoundingBox_' + str(pocket.getObjId()))
        f.write(PML_BBOX_STR_POCKSURF.format(outHETMFile, surfaceStr, toWrite))
      else:
        f.write(PML_SURF_STR.format(outHETMFile, surfaceStr))
    return pmlFile

  def createTCL(self):
    outHETMFile = os.path.abspath(self.getProteinHetatmFile())
    pqrFile = outHETMFile.replace('_out.pdb', '.pqr')
    tclFile = outHETMFile.replace('_out.pdb', '.tcl')
    with open(pqrFile, 'w') as f:
      for pocket in self:
        pqrPocket = getRawPDBStr(pocket.getFileName(), ter=False).strip()
        f.write(pqrPocket + '\n')
      f.write('TER\nEND')
    tclStr = TCL_STR % (outHETMFile, pqrFile)
    with open(tclFile, 'w') as f:
      f.write(tclStr)
    return tclFile

  def buildPocketsFiles(self, suffix='', tcl=False, bBox=False):
    outHETMFile = self.buildPDBhetatmFile(suffix)
    pmlFile = self.createPML(bBox=bBox)
    pmlFileSurf = self.createSurfacePml(bBox=bBox)
    if self.getPocketsClass() == 'FPocket' and tcl == True:
      self.createTCL()
    return outHETMFile, pmlFile, pmlFileSurf

  def getMAEFile(self):
    pock = self.getFirstItem()
    if hasattr(pock, '_maeFile'):
      return pock._maeFile.get()

  ######### Utils

  def getPocketsPDBStr(self):
    outStr = ''
    for i, pocket in enumerate(self):
      pocket.setObjId(i + 1)
      outStr += self.formatPocketStr(pocket)
    return outStr

  def formatPocketStr(self, pocket):
    outStr = ''
    numId, pocketFile = str(pocket.getObjId()), pocket.getFileName()
    rawStr = getRawPDBStr(pocketFile, ter=False).strip()
    if pocket.getPocketClass() in ['AutoLigand', 'AutoSite']:
      for line in rawStr.split('\n'):
        sline = splitPDBLine(line)
        replacements = ['HETATM', sline[1], 'APOL', 'STP', 'C', numId, *sline[6:-1], 'Ve']
        pdbLine = writePDBLine(replacements)
        outStr += pdbLine

    elif pocket.getPocketClass() == 'FPocket':
      for line in rawStr.split('\n'):
        sline = line.split()
        replacements = ['HETATM', sline[1], 'APOL', 'STP', 'C', numId, *sline[5:], '', 'Ve']
        try:
          pdbLine = writePDBLine(replacements)
        except:
          sline = splitPDBLine(line)
          replacements = ['HETATM', sline[1], 'APOL', 'STP', 'C', numId, *sline[6:], '', 'Ve']
          pdbLine = writePDBLine(replacements)
        outStr += pdbLine

    elif pocket.getPocketClass() == 'SiteMap':
      for line in rawStr.split('\n'):
        line = line.split()
        replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], '', 'Ve']
        pdbLine = writePDBLine(replacements)
        outStr += pdbLine

    else:  # (P2Rank, ElliPro, Standard...)
      for line in rawStr.split('\n'):
        line = splitPDBLine(line)
        line[5] = numId
        pdbLine = writePDBLine(line)
        outStr += pdbLine

    return outStr


class ProteinAtom(data.EMObject):
  def __init__(self, pdbLine, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.parseLine(pdbLine)

  def __str__(self):
    return 'Atom: {}. Type: {}'.format(self.atomId, self.atomType)

  def parseLine(self, pdbLine):
    if not pdbLine.startswith('ATOM') and not pdbLine.startswith('HETATM'):
      print('Passed line does not seems like an pdb ATOM line')
    else:
      self.line = pdbLine
      line = splitPDBLine(pdbLine)
      self.atomId = line[1]
      self.atomType = line[2]
      self.residueType = line[3]
      self.proteinChain = line[4]
      self.residueNumber = line[5]
      self.residueId = '{}.{}'.format(self.proteinChain, self.residueNumber)
      self.xCoord = line[6]
      self.yCoord = line[7]
      self.zCoord = line[8]
      self.atomLetter = line[-1]

  def getCoords(self):
    return tuple(map(float, (self.xCoord, self.yCoord, self.zCoord)))


class ProteinResidue(data.EMObject):
  def __init__(self, pdbLine, **kwargs):
    data.EMObject.__init__(self, **kwargs)
    self.parseLine(pdbLine)

  def __str__(self):
    return 'Residue: {}. Type: {}'.format(self.residueId, self.residueType)

  def parseLine(self, pdbLine):
    if not pdbLine.startswith('ATOM'):
      print('Passed line does not seems like an pdb ATOM line')
    else:
      line = splitPDBLine(pdbLine)
      self.residueType = line[3]
      self.proteinChain = line[4]
      self.residueNumber = line[5]
      self.residueId = '{}_{}'.format(self.proteinChain, self.residueNumber)


class MDSystem(data.EMFile):
  """A system atom structure (prepared for MD). Base class for Gromacs, Amber
        _topoFile: topology file
        _trjFile: trajectory file
        _ff: main force field
        _wff: water force field model"""

  def __init__(self, filename=None, **kwargs):
    super().__init__(filename=filename, **kwargs)
    self._oriStructFile = pwobj.String(kwargs.get('oriStructFile', None))
    self._topoFile = pwobj.String(kwargs.get('topoFile', None))
    self._trjFile = pwobj.String(kwargs.get('trjFile', None))
    self._ff = pwobj.String(kwargs.get('ff', None))
    self._wff = pwobj.String(kwargs.get('wff', None))

    self._ligName = pwobj.String(kwargs.get('ligName', None))
    self._topoLigFile = pwobj.String(kwargs.get('topoLigFile', None))

  def __str__(self):
    return '{} ({}, hasTrj={})'.format(self.getClassName(), os.path.basename(self.getSystemFile()),
                                       self.hasTrajectory())

  def getSystemFile(self):
    return self.getFileName()

  def setSystemFile(self, value):
    self.setFileName(value)

  def getTopologyFile(self):
    return self._topoFile.get()

  def setTopologyFile(self, value):
    self._topoFile.set(value)

  def getLigandName(self):
    return self._ligName.get()

  def setLigandName(self, value):
    self._ligName.set(value)

  def getLigandID(self):
    return 'LIG'

  def getLigandTopologyFile(self):
    return self._topoLigFile.get()

  def setLigandTopologyFile(self, value):
    self._topoLigFile.set(value)

  def hasTopology(self):
    if self.getTopologyFile():
      return True
    else:
      return False

  def hasTrajectory(self):
    if self.getTrajectoryFile():
      return True
    else:
      return False

  def getTrajectoryFile(self):
    return self._trjFile.get()

  def setTrajectoryFile(self, value):
    self._trjFile.set(value)

  def getOriStructFile(self):
    return self._oriStructFile.get()

  def setOriStructFile(self, value):
    self._oriStructFile.set(value)

  def getForceField(self):
    return self._ff.get()

  def setForceField(self, value):
    self._ff.set(value)

  def getWaterForceField(self):
    return self._wff.get()

  def setWaterForceField(self, value):
    self._wff.set(value)

  def getSystemName(self):
    return getBaseName(self.getSystemFile())


class PharmFeature(data.EMObject):
  """ Represent a pharmacophore feature """

  def __init__(self, **kwargs):
    super().__init__(**kwargs)
    self._type = String(kwargs.get('type', None))
    self._X, self._Y, self._Z = Float(kwargs.get('x', None)), Float(kwargs.get('y', None)), \
                                Float(kwargs.get('z', None))
    self._radius = Float(kwargs.get('radius', 1.0))

  def __str__(self):
    s = '{} {} (Type: {}. Coords: ({:.2f}, {:.2f}, {:.2f}). Radius: {:.2f})'. \
      format(self.getClassName(), self.getObjId(), self.getType(), *self.getCoords(), self.getRadius())
    return s

  def getType(self):
    return self._type.get()

  def setType(self, value):
    self._type.set(value)

  def getRadius(self):
    return self._radius.get()

  def setRadius(self, value):
    self._radius.set(value)

  def getCoords(self):
    return self._X.get(), self._Y.get(), self._Z.get()

  def setCoords(self, values):
    self._X.set(values[0]), self._Y.set(values[1]), self._Z.set(values[2])

  def feat2Dic(self):
    return {'type': self.getType(), 'coords': self.getCoords(), 'radius': self.getRadius()}

  def setFeatFromDic(self, featDic):
    if 'type' in featDic:
      self.setType(featDic['type'])
    else:
      print('Type for {} has not been specified'.format(self))

    if 'coords' in featDic:
      self.setCoords(featDic['coords'])
    else:
      print('Coordinates for {} has not been specified'.format(self))

    if 'radius' in featDic:
      self.setRadius(featDic['radius'])
    else:
      self.setRadius(1.0)


class PharmacophoreChem(data.EMSet):
  """ Pharmacophore (built as a set of PharmFeature) """
  ITEM_TYPE = PharmFeature

  def __init__(self, **kwargs):
    data.EMSet.__init__(self, **kwargs)
    self._proteinFile = pwobj.String(kwargs.get('proteinFile', None))

  def __str__(self):
    s = '{} ({} features)'.format(self.getClassName(), self.getSize())
    return s

  def clone(self):
    clone = self.getClass()()
    clone.copy(self)
    clone.copyInfo(self)
    return clone

  def copyInfo(self, other):
    self._proteinFile = other._proteinFile

  def getSetPath(self):
    return os.path.abspath(self._mapperPath[0])

  def getSetDir(self, path=''):
    return os.path.join('/'.join(self.getSetPath().split('/')[:-1]), path)

  def setProteinFile(self, value):
    self._proteinFile.set(value)

  def getProteinFile(self):
    return self._proteinFile.get()

  def pharm2Dic(self):
    pDic = {}
    for feat in self:
      pDic[feat.getObjId()] = feat.feat2Dic()
    return pDic

  def setPharmFromDic(self, pDic):
    for featId in pDic:
      feat = PharmFeature().setFeatFromDic(pDic[featId])
      feat.setObjId(featId)
      self.append(feat)


############################################################
##############  POSSIBLE OUTPUTS OBJECTS ###################
############################################################

class PredictStructROIsOutput(enum.Enum):
  outputStructROIs = SetOfStructROIs