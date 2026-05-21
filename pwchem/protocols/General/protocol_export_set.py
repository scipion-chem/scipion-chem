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

import os, importlib, shutil

from pwem.protocols import EMProtocol
from pwem.objects import EMSet, String

from pyworkflow.protocol import params

from pwchem.utils import getBaseFileName

IMP, EXP = 0, 1
SET_NOCOPY = ['_size', '_streamState', '_mapperPath']
TYPE_STR = 'obj_type'

class ProtChemImportExportSet(EMProtocol):
    """
    Protocol to export and import Scipion objects and datasets via filesystem serialization.

    This protocol allows saving a Scipion EMObject (including EMSet collections)
    into a directory structure and later reconstructing it from those files.

    It supports two modes:
    - Export: serialize object attributes and items into files
    - Import: reconstruct objects from previously exported files

    Inputs
    ------
    mode:
        Operation mode:
        - Import: reconstruct object from directory
        - Export: serialize object into directory

    input:
        Input EMObject to export (only required in Export mode).

    directory:
        Target directory for export or source directory for import.
        If not specified during export, a default directory is created
        inside the protocol working directory.

    Export mode
    -----------
    The protocol serializes the object into three files:

    1. outputType.txt
       - Stores the full class path of the object.

    2. outputAttributes.csv
       - First line: attribute names
       - Second line: attribute values

    3. items.csv (only for EMSet objects)
       - First line: attribute names (+ obj_type)
       - Second line: attribute types
       - Remaining lines: attribute values for each item

    Workflow (Export)
    -----------------
    1. Create output directory.
    2. Save object type.
    3. Save object attributes.
    4. If object is a set:
       - Iterate over items
       - Save attributes and types
       - Store file paths as absolute paths

    Import mode
    -----------
    The protocol reconstructs the object from exported files.

    Workflow (Import)
    -----------------
    1. Read object type from outputType.txt.
    2. Dynamically import class using importlib.
    3. Parse object attributes from outputAttributes.csv.

    4. If object is EMSet:
       - Create empty set
       - Parse items.csv
       - Reconstruct each item:
         * Restore attribute types
         * Instantiate correct object class
         * Assign attributes
         * Copy referenced files if needed

    5. If object is not a set:
       - Instantiate object
       - Assign attributes directly

    Output
    ------
    importedSet:
        Reconstructed dataset (if imported object is EMSet)

    importedObject:
        Reconstructed object (if not EMSet)

    Error handling
    --------------
    - Skips missing directories silently.
    - Copies external files into project if needed.
    - Assumes exported files follow expected format.

    Validation
    ----------
    - Requires valid directory structure for import.
    - Requires valid class paths for dynamic import.
    - Attributes must be compatible with String casting or original types.

    Summary
    -------
    This protocol enables:
    - Exporting Scipion objects for portability
    - Sharing datasets across projects
    - Reconstructing workflows from serialized outputs

    It provides a lightweight alternative to full project export/import
    by focusing on individual objects and datasets.

    Notes
    -----
    - Uses dynamic class loading (importlib).
    - File paths are preserved and copied if necessary.
    - Some internal attributes (e.g. _size, _mapperPath) are not copied.
    - Attribute types are reconstructed based on stored type metadata.
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
                val = v.get()
                if type(v.get()) == str and os.path.isfile(v.get()) and os.getcwd() not in v.get():
                    # Copy files if import in different project
                    oDir = self._getExtraPath(k)
                    if not os.path.exists(oDir):
                        os.mkdir(oDir)
                    nVal = os.path.join(oDir, getBaseFileName(val))
                    shutil.copy(val, nVal)
                    v = String(nVal)
                setattr(obj, k, v)
            oSet.append(obj)
        return oSet

    def parseItems(self, itemFile):
        with open(itemFile) as f:
            attrKeys = [k.strip() for k in f.readline().strip().split(';')]
            attrTypes = list(f.readline().strip().split(';'))

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
        return dict(zip(keys, values))

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
                    val = str(value.get())
                    if os.path.isfile(val):
                        val = os.path.abspath(val)
                    line += val
                if i == 0:
                    fh.write(lineHdr + f"; {TYPE_STR}\n")
                    fh.write(lineTypes + "; ...\n")
                fh.write(line + f"; {type(entry)}\n")

    def exportObjectInfo(self, oDir):
        oFile = self.getObjectsFile(oDir)[0]
        with open(oFile, 'w') as f:
            f.write(str(type(self.input.get())))

        oFile = self.getObjectsFile(oDir)[1]
        heads, lines = [], []
        for key, value in self.input.get().getAttributes():
            val = str(value.get())
            if os.path.isfile(val):
                val = os.path.abspath(val)
            heads.append(key), lines.append(val)

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