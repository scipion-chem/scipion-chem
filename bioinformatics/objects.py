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
