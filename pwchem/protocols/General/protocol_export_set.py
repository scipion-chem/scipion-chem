# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, importlib

from pwem.protocols import EMProtocol
from pwem.objects import EMSet, String

from pyworkflow.protocol import params

IMP, EXP = 0, 1
SET_NOCOPY = ['_size', '_streamState', '_mapperPath']
TYPE_STR = 'obj_type'

class ProtChemImportExportSet(EMProtocol):
    """Export a set by copying the elements in a desired directory and import from those created directories
    """
    _label = 'Import/export output'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('mode', params.EnumParam, choices=['Import', 'Export'], label='Protocol mode: ',
                      display=params.EnumParam.DISPLAY_HLIST, default=0,
                      help='Whether to use this protocol to import or to export an object via files')
        form.addParam('input', params.PointerParam, pointerClass="EMObject", label='Input object: ',
                      condition='mode==1',
                      help='Input object whose attributes will be copied to the desired directory')

        form.addParam('directory', params.FileParam, label='Import/export directory: ', default='',
                      help='Directory where the information of the object will be saved '
                           'By default, this protocol directory/exportedObject_objId')


    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        if self.mode.get() == IMP:
            self._insertFunctionStep(self.importStep)
        else:
            self._insertFunctionStep(self.exportStep)

    def importStep(self):
        iDir = self.directory.get()
        if os.path.exists(iDir):
            oTypeFile, oAttrFile, oItemFile = self.getObjectsFile(iDir)
            with open(oTypeFile) as f:
                typePath = f.read().strip()
            oClass = self.importObjectType(typePath)
            attrDic = self.parseObjAttributes(oAttrFile)

            if issubclass(oClass, EMSet):
                oSet = self.buildSet(oClass, attrDic)
                oSet = self.appendImportedItems(oSet, oItemFile)
                self._defineOutputs(importedSet=oSet)

            else:
                oObj = oClass()
                for k, val in attrDic.items():
                    setattr(oObj, k, String(val))

                self._defineOutputs(importedObject=oObj)

    def exportStep(self):
        oDir = self.getOutputDirectory()
        self.exportObjectInfo(oDir)
        if isinstance(self.input.get(), EMSet):
            self.exportItems(oDir)


    ##################### UTILS FUNCTIONS ###################

    def buildSet(self, oClass, attrDic):
        oSet = oClass().create(outputPath=self._getPath())
        for k, val in attrDic.items():
            if k not in SET_NOCOPY:
                setattr(oSet, k, String(val))
        return oSet

    def appendImportedItems(self, oSet, oItemFile):
        attrKeys, attrValues, outObjs = self.parseItems(oItemFile)
        for i, obj in enumerate(outObjs):
            objValues = attrValues[i]
            for k, v in zip(attrKeys, objValues):
                setattr(obj, k, v)
            oSet.append(obj)
        return oSet

    def parseItems(self, itemFile):
        with open(itemFile) as f:
            attrKeys = [k.strip() for k in f.readline().strip().split(';')]
            attrTypes = [t for t in f.readline().strip().split(';')]

            attributeValues, outObjs = [], []
            for line in f:
                attributeValues.append([])
                values = line.strip().split(';')
                for k, t, v in zip(attrKeys, attrTypes, values):
                    k, v = k.strip(), v.strip()
                    if k != TYPE_STR:
                        t = self.importObjectType(t)
                        if v == 'None':
                            v = None
                        attributeValues[-1].append(t(v))
                    else:
                        outObjs.append(self.importObjectType(v)())
        return attrKeys[:-1], attributeValues, outObjs

    def parseObjAttributes(self, attrFile):
        with open(attrFile) as f:
            keys = f.readline().strip().split(';')
            values = f.readline().strip().split(';')
        return {key: val for key, val in zip(keys, values)}

    def importObjectType(self, pathString):
        # Parse the class path from the string
        classPath = pathString.strip("<>").split(" '")[1].strip("'")
        modulePath, className = classPath.rsplit('.', 1)

        # Import the module dynamically
        module = importlib.import_module(modulePath)

        # Get the class
        oClass = getattr(module, className)
        return oClass

    def exportItems(self, oDir):
        oFile = self.getObjectsFile(oDir)[2]
        with open(oFile, "w") as fh:
            lineHdr, lineTypes = "", ""
            for i, entry in enumerate(self.input.get()):
                line = ""
                for key, value in entry.getAttributes():
                    if i == 0:
                        if lineHdr != "":
                            lineHdr += "; "
                            lineTypes += "; "
                        lineHdr += key
                        lineTypes += str(type(value))

                    if line != "":
                        line += "; "
                    line += str(value.get())
                if i == 0:
                    fh.write(lineHdr + f"; {TYPE_STR}\n")
                    fh.write(lineTypes + f"; ...\n")
                fh.write(line + f"; {type(entry)}\n")

    def exportObjectInfo(self, oDir):
        oFile = self.getObjectsFile(oDir)[0]
        with open(oFile, 'w') as f:
            f.write(str(type(self.input.get())))

        oFile = self.getObjectsFile(oDir)[1]
        heads, lines = [], []
        for key, value in self.input.get().getAttributes():
            heads.append(key), lines.append(str(value.get()))

        with open(oFile, 'w') as f:
            f.write(';'.join(heads) + '\n')
            f.write(';'.join(lines))

    def getOutputDirectory(self):
        oDir = self.directory.get()
        if not oDir.strip():
            oDir = self._getPath()
        if os.path.exists(oDir):
            oDir = os.path.join(oDir, f'exportedObject_{self.input.get().getObjId()}')

        if not os.path.exists(oDir):
            os.mkdir(oDir)
        return oDir

    def getObjectsFile(self, oDir):
        objTypeFile = os.path.join(oDir, 'outputType.txt')
        objAttrFile = os.path.join(oDir, 'outputAttributes.csv')
        objItemFile = os.path.join(oDir, 'items.csv')
        return objTypeFile, objAttrFile, objItemFile