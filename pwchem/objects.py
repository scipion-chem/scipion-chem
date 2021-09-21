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

import pyworkflow.object as pwobj
import pwem.objects.data as data
from pyworkflow.object import (Float, Integer, List, String)
import numpy as np
import os
from scipy import spatial
from .utils import *
from .constants import *

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

class ProteinSequenceFile(data.EMFile):
    """A file with a list of protein sequences"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

class NucleotideSequenceFile(data.EMFile):
    """A file with a list of nucleotide sequences"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

class SmallMolecule(data.EMObject):
    """ Small molecule """
    def __init__(self, **kwargs):
        data.EMObject.__init__(self, **kwargs)
        self.smallMoleculeFile = pwobj.String(kwargs.get('smallMolFilename', None))

    def getFileName(self):
        return self.smallMoleculeFile.get()

    def getConformersFileName(self):
        return self._ConformersFile.get()

    def getParamsFileName(self):
        return self._ParamsFile.get()

    def getPDBFileName(self):
        return self._PDBFile.get()

class SetOfSmallMolecules(data.EMSet):
    """ Set of Small molecules """
    ITEM_TYPE = SmallMolecule
    FILE_TEMPLATE_NAME = 'setOfSmallMolecules%s.sqlite'

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)

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

class ProteinPocket(data.EMFile):
    """ Represent a pocket file """

    def __init__(self, filename=None, proteinFile=None, extraFile=None, **kwargs):
        data.EMFile.__init__(self, filename, **kwargs)
        self._proteinFile = String(proteinFile)
        self._extraFile = String(extraFile)
        self._nPoints = Integer(kwargs.get('nPoints', None))
        self._contactResidues = String(kwargs.get('contactResidues', None))
        self._contactAtoms = String(kwargs.get('contactAtoms', None))
        self._volume = Float(kwargs.get('volume', None))
        self._score = Float(kwargs.get('score', None))
        self._energy = Float(kwargs.get('energy', None))
        self._class = String(kwargs.get('class', None))

    #Attributes functions
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

    def getKwargs(self, props, AM):
        nkwargs = {}
        for k in props:
            if k in AM:
                nkwargs[AM[k]] = props[k]
        return nkwargs

    #Complex pocket attributes functions
    def buildContactAtoms(self, calculate=False, maxDistance=4):
        '''Return the reported proteins atoms in contact with the pocket.
        If not reported, returns the protein atoms at less than 4A than any pocket point'''
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
            #Manually calculate the contacts
            pocketCoords = self.getPointsCoords()
            proteinAtoms = self.getProteinAtoms()
            proteinCoords = self.getAtomsCoords(proteinAtoms)
            dists = spatial.distance.cdist(proteinCoords, pocketCoords)
            for i in range(len(proteinCoords)):
                if min(dists[i,:]) < maxDistance:
                    contactAtoms.append(proteinAtoms[i])
        return contactAtoms

    def calculateMassCenter(self, weights=None):
        '''Calculates the center of mass of a set of points: [(x,y,z), (x,y,z),...]
        A weight for each point can be specified'''
        coords = self.getPointsCoords()
        return list(np.average(coords, axis=0, weights=weights))

    def getDiameter(self, radius=[]):
        '''Returning max distance of points found in the convex hull
        Possibility of adding radius to each row and column to calculate the diameter
        on spheres instead of points'''
        coords = np.array(self.getPointsCoords())
        cHullIndex = spatial.ConvexHull(coords).vertices
        candidates = coords[cHullIndex]
        dist_mat = spatial.distance_matrix(candidates, candidates)
        if radius!=[]:
            dist_mat = self.addRadius(dist_mat, radius[cHullIndex])
        i, j = np.unravel_index(dist_mat.argmax(), dist_mat.shape)

        return dist_mat[i, j]

    def addRadius(self, dMat, radius):
        '''Add the radius of each alpha sphere to their corresponding row and column in the distances
        matrix'''
        for i in range(dMat.shape[0]):
            dMat[i,:] += radius[i]
            dMat[:, i] += radius[i]
        return dMat

    def getPocketVolume(self, mode=0):
        '''Return the pocket volume:
        0: getSurfaceConvexVolume (over protein contact points)
        1: reported volume
        2: getConvexVolume (over pocket points)
        '''
        if mode==2:
            return self.getConvexVolume()
        elif mode==1 and float(self.getVolume()) != None:
            return self.getVolume()
        else:
            return self.getSurfaceConvexVolume()

    def getConvexVolume(self):
        '''Calculate the convex volume of the points forming the pocket'''
        coords = np.array(self.getPointsCoords())
        cHull = spatial.ConvexHull(coords)
        return cHull.volume

    def getSurfaceConvexVolume(self):
        '''Calculate the convex volume of the protein contact atoms'''
        cAtoms = self.buildContactAtoms()
        cCoords = self.getAtomsCoords(cAtoms)
        cHull = spatial.ConvexHull(cCoords)
        return cHull.volume

    def getPocketBox(self):
        '''Return the coordinates of the 2 corners determining the box (ortogonal to axis) where the pocket fits in
        For example: ([-1,0,2], [2,5,4])'''
        coords = np.array(self.getPointsCoords())
        return coords.min(axis=0), coords.max(axis=0)

    #Utils functions
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
        coords=[]
        for atom in atoms:
            coords.append(atom.getCoords())
        return coords

    def getAtomsIds(self, atoms):
        ids = []
        for atom in atoms:
            ids.append(atom.atomId)
        return ids

    def getResiduesIds(self, residues):
        ids = set([])
        for res in residues:
            ids.add(res.residueId)
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


class SetOfPockets(data.EMSet):
    ITEM_TYPE = ProteinPocket

    def __init__(self, **kwargs):
        data.EMSet.__init__(self, **kwargs)
        self._pocketsClass = String(kwargs.get('pocketsClass', None))
        self._hetatmFile = String(kwargs.get('hetatmFile', None))
        self._pmlFile = String(kwargs.get('pmlFile', None))
        self._pmlFileSurf = String(kwargs.get('pmlFileSurf', None))

    def __str__(self):
      s = '{} ({} items, {} class)'.format(self.getClassName(), self.getSize(), self.getPocketsClass())
      return s

    def getSetPath(self):
        return os.path.abspath(self._mapperPath[0])

    def getSetDir(self):
        return '/'.join(self.getSetPath().split('/')[:-1])

    def getProteinName(self):
        return self.getProteinFile().split('/')[-1].split('.')[0]

    def getPmlFile(self):
        return self._pmlFile.get()

    def setPmlFile(self, value):
        self._pmlFile.set(value)

    def getPmlFileSurf(self):
        return self._pmlFileSurf.get()

    def setPmlFileSurf(self, value):
        self._pmlFileSurf.set(value)

    def getProteinFile(self):
        return self.getFirstItem().getProteinFile()

    def getProteinHetatmFile(self):
        return self._hetatmFile.get()

    def setProteinHetatmFile(self, value):
        self._hetatmFile.set(value)

    def getPocketsClass(self):
        return self._pocketsClass.get()

    def updatePocketsClass(self):
        pClass = self.getPocketsClass()
        for pocket in self:
            if not pocket.getPocketClass() == pClass:
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
        outDir = self.getSetDir()
        outFile = os.path.join(outDir, protName+'{}_out.pdb'.format(suffix))

        with open(outFile, 'w') as f:
            f.write(getRawPDBStr(atmFile, ter=False))
            f.write(self.getPocketsPDBStr())

        self.setProteinHetatmFile(outFile)
        return outFile

    def createPML(self, outHETMFile):
        outHETMFile = os.path.abspath(outHETMFile)
        pmlFile = outHETMFile.replace('_out.pdb', '.pml')

        # Creates the pml for pymol visualization
        with open(pmlFile, 'w') as f:
            f.write(PML_STR.format(outHETMFile))
        self.setPmlFile(pmlFile)
        return pmlFile

    def createSurfacePml(self, outHETMFile):
        outHETMFile = os.path.abspath(outHETMFile)
        pmlFile = outHETMFile.replace('_out.pdb', '_surf.pml')
        colors = createColorVectors(len(self))
        surfaceStr = ''
        for i, pock in enumerate(self):
            pId = pock.getObjId()
            surfAtomIds = str(list(map(int, pock.getDecodedCAtoms()))).replace(' ', '')
            surfaceStr += PML_SURF_EACH.format(pId, colors[i], pId, surfAtomIds, pId, pId)

        # Creates the pml for pymol visualization
        with open(pmlFile, 'w') as f:
            f.write(PML_SURF_STR.format(outHETMFile, surfaceStr))
        self.setPmlFileSurf(pmlFile)
        return pmlFile

    def createTCL(self, outHETMFile):
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


    def buildPocketsFiles(self, suffix='', tcl=False):
        outHETMFile = self.buildPDBhetatmFile(suffix)
        self.createPML(outHETMFile)
        self.createSurfacePml(outHETMFile)
        if self.getPocketsClass() == 'FPocket' and tcl==True:
            self.createTCL(outHETMFile)
        return outHETMFile

    ######### Utils

    def getPocketsPDBStr(self):
        outStr = ''
        for i, pocket in enumerate(self):
            pocket.setObjId(i+1)
            outStr += self.formatPocketStr(pocket)
        return outStr

    def formatPocketStr(self, pocket):
        outStr = ''
        numId, pocketFile = str(pocket.getObjId()), pocket.getFileName()
        rawStr = getRawPDBStr(pocketFile, ter=False).strip()
        if pocket.getPocketClass() == 'AutoLigand':
            for line in rawStr.split('\n'):
                line = line.split()
                replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], 'Ve']
                pdbLine = writePDBLine(replacements)
                outStr += pdbLine

        elif pocket.getPocketClass() == 'P2Rank':
            for line in rawStr.split('\n'):
                line = splitPDBLine(line)
                line[5] = numId
                pdbLine = writePDBLine(line)
                outStr += pdbLine

        elif pocket.getPocketClass() == 'FPocket':
            for line in rawStr.split('\n'):
                line = line.split()
                replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:], '', 'Ve']
                pdbLine = writePDBLine(replacements)
                outStr += pdbLine

        elif pocket.getPocketClass() == 'SiteMap':
            for line in rawStr.split('\n'):
                line = line.split()
                replacements = ['HETATM', line[1], 'APOL', 'STP', 'C', numId, *line[5:-1], '', 'Ve']
                pdbLine = writePDBLine(replacements)
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
