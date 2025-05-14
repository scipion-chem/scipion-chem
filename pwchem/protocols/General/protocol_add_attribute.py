# **************************************************************************
# *
# * Authors:		Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
This module will extract the ligand from a complex pdb.
"""

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.utils import Message

# Scipion chem imports
from ...utils import String

OBJECT, ITEM = 0, 1

class ProtAddAttribute(EMProtocol):
  _label = 'Add attribute'

  def _defineParams(self, form):
    # From input condition string
    fromInputCondition = 'fromInput=='

    form.addSection(label=Message.LABEL_INPUT)
    group = form.addGroup('Attribute')
    group.addParam('atKey', params.StringParam, label="Attribute name: ", default='',
                    help='Name of the new attribute of the set')
    group.addParam('atType', params.EnumParam, label="Attribute type: ", default=0,
                    choices=['String', 'Integer', 'Float'], expertLevel=params.LEVEL_ADVANCED,
                    help='Type of the attribute to add')

    group = form.addGroup('Input')
    group.addParam('fromInput', params.EnumParam, label="Input to label: ", default=0, choices=['Object', 'Items'],
                    help='Whether to add the attribute to the input object or to the items inside it '
                      '(in case it is a set)')

    group.addParam('inputObject', params.PointerParam, label="Input object: ",
                    important=True, pointerClass='EMObject', condition=f'{fromInputCondition}0',
                    help='Object you want to label with a new attribute')
    group.addParam('atValue', params.StringParam, label="Attribute value: ", default='', condition=f'{fromInputCondition}0',
                    help='Value of the new attribute of the object')

    group.addParam('inputSet', params.PointerParam, label="Input set: ",
                    important=True, pointerClass='EMSet', condition=f'{fromInputCondition}1',
                    help='Set containing the items you want to label with a new attribute')
    group.addParam('mapKey', params.StringParam, label="Mapping key: ", default='', condition=f'{fromInputCondition}1',
                    help='Key of each item attribute to map the correspondance. '
                      'E.g: use the item ID to select the new attribute value for each item')

    group.addParam('inputFile', params.FileParam, label="Map attribute file: ", condition=f'{fromInputCondition}1',
                    help='File containing the mapping of the items to the attribute values. It must be a csv '
                      'file where each line contains the "key, value"')

  # --------------------------- Steps functions --------------------
  def _insertAllSteps(self):
    self._insertFunctionStep('createOutputStep')

  def createOutputStep(self):
    atKey, atValue, atType = self.atKey.get(), self.atValue.get(), self.getEnumText('atType')

    if self.fromInput.get() == OBJECT:
      inpObj = self.inputObject.get()
      outObj = inpObj.clone()
      outObj = self.setAttribute(outObj, atKey, atValue, atType)

      self._defineOutputs(outputObject=outObj)

    elif self.fromInput.get() == ITEM:
      inpSet = self.inputSet.get()
      mapDic = self.getMapDic()

      outSet = inpSet.createCopy(self._getPath(), copyInfo=True)
      for item in inpSet:
        itemKey = str(getattr(item, self.mapKey.get()))
        if itemKey in mapDic:
          atValue = mapDic[itemKey]
          newItem = self.setAttribute(item, atKey, atValue, atType)
          outSet.append(newItem)
        else:
          print('Mapping key {} was not found in the mapping file for attribute {}.'
                .format(itemKey, self.mapKey.get()))

      self._defineOutputs(outputSet=outSet)

  # --------------------------- INFO functions -----------------------------------
  def _validate(self):
    vals = []
    if self.fromInput.get() == ITEM:
      for item in self.inputSet.get():
        if not hasattr(item, self.mapKey.get()):
          vals.append('Mapping attribute {} is not found in the items of the input set'.format(self.mapKey.get()))
          break
    return vals

  def setAttribute(self, obj, atKey, atValue, atType):
    obj.__setattr__(atKey, eval(atType)(atValue))
    return obj

  def getMapDic(self):
    mapDic = {}
    with open(self.inputFile.get()) as f:
      for line in f:
        sline = line.split(',')
        mapDic[sline[0]] = sline[1].strip()
    return mapDic