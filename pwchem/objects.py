# -*- coding: utf-8 -*-
#  **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

class AutodockGrid(data.EMFile):
    """A search grid in the file format of Autodock"""
    def __init__(self, **kwargs):
        data.EMFile.__init__(self, **kwargs)

class proteinPocket(data.EMFile):
  """ Represent a pocket file """

  def __init__(self, filename=None, **kwargs):
    data.EMFile.__init__(self, filename, **kwargs)

class fpocketPocket(proteinPocket):
  """ Represent a pocket file from fpocket"""
  def __init__(self, filename=None, **kwargs):
    proteinPocket.__init__(self, filename, **kwargs)
    self.properties, self.pocketId = self.parseFile(filename)
    self.setObjId(self.pocketId)

  def __str__(self):
    s = 'Fpocket pocket {}\nFile: {}'.format(self.pocketId, self.getFileName())
    return s

  def getVolume(self):
    return self.properties['Real volume (approximation)']

  def getPocketScore(self):
    return self.properties['Pocket Score']

  def getNumberOfVertices(self):
    return self.properties['Number of V. Vertices']

  def parseFile(self, filename):
    props = {}
    ini, parse = 'HEADER Information', False
    with open(filename) as f:
      for line in f:
        if line.startswith(ini):
          parse=True
          pocketId = int(line.split()[-1].replace(':', ''))
        elif line.startswith('HEADER') and parse:
          name = line.split('-')[1].split(':')[0]
          val = line.split(':')[-1]
          props[name.strip()] = float(val.strip())
    return props, pocketId

class autoLigandPocket(proteinPocket):
    """ Represent a pocket file from autoLigand"""

    def __init__(self, filename=None, resultsFile=None, **kwargs):
        proteinPocket.__init__(self, filename, **kwargs)
        print(filename)
        print(filename.split('out'))
        print(filename.split('out')[1].split('.'))
        self.pocketId = int(filename.split('out')[1].split('.')[0])
        self.properties = self.parseFile(resultsFile)
        self.setObjId(self.pocketId)

    def __str__(self):
        s = 'Fpocket pocket {}\nFile: {}'.format(self.pocketId, self.getFileName())
        return s

    def getVolume(self):
        return self.properties['Total Volume']

    def getEnergyPerVol(self):
        return self.properties['Total Energy per Vol']

    def parseFile(self, filename):
        props, i = {}, 1
        with open(filename) as f:
            for line in f:
                if i==self.pocketId:
                    line = line.split(',')
                    print(line)
                    #Volume
                    props[line[1].split('=')[0].strip()] = float(line[1].split('=')[1].strip())
                    #Energy/vol
                    props[line[2].strip()] = float(line[3].split('=')[1].strip())

        return props

class SetOfPockets(data.EMSet):
    ITEM_TYPE = proteinPocket

    def __str__(self):
      s = '{} ({} items)'.format(self.getClassName(), self.getSize())
      return s
