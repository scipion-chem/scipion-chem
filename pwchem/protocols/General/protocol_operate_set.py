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

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.utils import fillEmptyAttributes

UNIQUE, UNION, INTERSECTION, DIFFERENCE, FILTER, REMCOl, RANK = list(range(7))

class ProtChemOperateSet(EMProtocol):
    """Filter a set by a column value or keep just a few columns"""
    _label = 'operate set'

    def _defineParams(self, form):
        form.addSection(label='Input')
        group = form.addGroup('Operation')
        group.addParam('operation', params.EnumParam, label='Operation: ', default=0,
                      choices=['Unique', 'Union', 'Intersection', 'Difference', 'Filter', 'Remove columns',  'Ranking'],
                      help='-Sets operations: Duplicates share the same reference column value\n'
                           '\tUnique: keep just one item with the same reference column value. Similar to remove '
                           'duplicates.\n\tUnion: merges two or more sets.\n\tIntersection: keep only the items '
                           'repeated in all the input sets.\n\tDifference: keep only the items in the first set that '
                           'are not present in the second.\n\n-Modification operations:\n\tFilter: outputs only '
                           'those items passing the filter\n\tRemove columns: remove the specified columns\n\t'
                           'Ranking: outputs only the top/bottom elements for the specified column.')

        group.addParam('removeDuplicates', params.BooleanParam, default=False,
                       label='Remove duplicates: ', condition='not operation in [0]',
                       help='Remove elements with the reference column value repeated')
        group.addParam('refColumn', params.StringParam, label='Reference column: ', default='',
                       condition='removeDuplicates or operation in [0, 1, 2, 3]',
                       help='Reference attribute for the set operation.')

        group.addParam('filterColumn', params.StringParam, label='Filter column: ', condition='(operation in [4, 6])',
                       help='Attribute for the set filtering.')
        group.addParam('filterOp', params.EnumParam, label='Filter operation: ',
                      condition='(operation==4)', default=0,
                      choices=['==', '>', '>=', '<', '<=', '!=', 'between', 'startswith', 'endswith', 'contains', 
                               'does not startwith', 'does not end with', 'does not contain'])
        group.addParam('filterValue', params.StringParam,
                       label='Value: ', condition='(operation==4)',
                       help='Value to use in the filter')
        group.addParam('filterValue2', params.StringParam,
                       label='Lower Value: ', condition='(operation==4 and filterOp==6)',
                       help='Value to use in the filter')
        
        group.addParam('remColumns', params.StringParam, label='Remove columns: ', condition='operation==5',
                       help='They must exist in the input database list. Separated by semicolons '
                            '(e.g. column1 ; column2 ; ...)')
        group.addParam('threshold', params.StringParam, label="Threshold: ", condition="operation==6",
                      help='Number/proportion of items to keep:\n\tNumber: n>=1 \n\tProportion: 0<n<1\n\tPercentage: n%\n\n'
                           'Higher/lower values of the attribute: \n\tHigher: positive number\n\tLower: negative number\n\n'
                           'e.g: "-10%" == "-0.1" == 10% of the items with lower values\n'
                           'e.g: "5" == 5 items with higher values')

        group = form.addGroup('Input')
        group.addParam('inputSet', params.PointerParam, pointerClass="EMSet", label='Set to filter: ',
                       condition='not operation in [1, 2]', help='Principal set to operate')
        group.addParam('inputMultiSet', params.MultiPointerParam, pointerClass="EMSet", allowsNull=True,
                      label='Input sets: ', condition='operation in [1, 2]')
        group.addParam('secondSet', params.PointerParam, pointerClass="EMSet", condition="operation==3",
                      label='Set with items to remove: ',
                      help='Secondary set to operate')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('operateStep')

    def operateStep(self):
        outputDict, tmpDict = {}, {}
        opAttr = self.refColumn.get() if self.refColumn.get() else '_objId'
        
        if self.operation.get() == UNIQUE:
            for item in self.inputSet.get():
                opId = self.getAttrValue(item, opAttr)
                self.addItem(outputDict, opId, item, allowDup=False)
                
        elif self.operation.get() == UNION:
            inputSets = fillEmptyAttributes(self.inputMultiSet)
            for inSet in inputSets:
                for item in inSet.get():
                    opId = self.getAttrValue(item, opAttr)
                    self.addItem(outputDict, opId, item)
                  
        elif self.operation.get() == INTERSECTION:
            idsLists = []
            inputSets = fillEmptyAttributes(self.inputMultiSet)
            for inSet in inputSets:
                curList = []
                for item in inSet.get():
                    opId = self.getAttrValue(item, opAttr)
                    tmpDict[opId] = [] if not opId in tmpDict else tmpDict[opId] + [item.clone()]
                    curList.append(opId)
                idsLists.append(curList)

            interIds = set(idsLists[0])
            for idList in idsLists[1:]:
                interIds = interIds.intersection(set(idList))

            for interId in interIds:
                for item in tmpDict[interId]:
                    self.addItem(outputDict, interId, item)

        elif self.operation.get() == DIFFERENCE:
            removeIds = []
            for item in self.secondSet.get():
                opId = self.getAttrValue(item, opAttr)
                removeIds.append(opId)

            for item in self.inputSet.get():
                opId = self.getAttrValue(item, opAttr)
                if not opId in removeIds:
                    self.addItem(outputDict, opId, item)
        
        elif self.operation.get()==FILTER:
            filterOp = self.filterOp.get()
            referenceValue = self.getRefValue()
            
            if filterOp == 6:
                referenceValue2 = self.getRefValue(2)

            for item in self.inputSet.get():
                opId = self.getAttrValue(item, opAttr)
                value = item.getAttributeValue(self.filterColumn.get())

                add = False
                if filterOp in [0, 1, 2, 3, 4, 5]: # ==, >, ...
                    add = eval('{}{}{}'.format(value, self.getEnumText('filterOp'), referenceValue))

                elif filterOp == 6:  # between
                    add = (value <= referenceValue and value >= referenceValue2)
                elif filterOp == 7:  #startswith
                    add = value.startswith(referenceValue)
                elif filterOp == 8:  # endswith
                    add = value.endswith(referenceValue)
                elif filterOp == 9:  # contains
                    add = referenceValue in value
                elif filterOp == 10:  # does not startswith
                    add = not (value.startswith(referenceValue))
                elif filterOp == 11:  # does not endswith
                    add = not (value.endswith(referenceValue))
                elif filterOp == 12:  # does not contains
                    add = not (referenceValue in value)

                if add:
                    self.addItem(outputDict, opId, item)

        elif self.operation.get()==REMCOl:
            remList=[x.strip() for x in self.remColumns.get().split(';')]

            ignoreList=[]
            for name, _ in self.inputSet.get().getFirstItem().getAttributes():
                if name in remList:
                    ignoreList.append(name)

            for item in self.inputSet.get():
                newEntry = self.inputSet.get().ITEM_TYPE()
                newEntry.copy(item, ignoreAttrs=ignoreList)

                opId = self.getAttrValue(newEntry, opAttr)
                self.addItem(outputDict, opId, newEntry)

        elif self.operation.get() == RANK:
            outNumber, descending = self.parseTopRankParam()

            V = []
            for entry in self.inputSet.get():
                V.append(entry.getAttributeValue(self.filterColumn.get()))
            V.sort(reverse=descending)

            threshold = V[outNumber]

            for item in self.inputSet.get():
                opId = self.getAttrValue(item, opAttr)
                value = item.getAttributeValue(self.filterColumn.get())
                if descending and value>threshold:
                    self.addItem(outputDict, opId, item)
                elif not descending and value<threshold:
                    self.addItem(outputDict, opId, item)

        if len(outputDict)>0:
            outputSet = self.getRepInputSet().createCopy(self._getPath(), copyInfo=True)
            i = 1
            for itemId in outputDict:
                for item in outputDict[itemId]:
                    if self.operation.get() in [UNION, INTERSECTION]:
                        item.setObjId(i)
                        i += 1
                    outputSet.append(item)
            self._defineOutputs(outputSet=outputSet)


    ######################## UTILS functions #################

    def getRepInputSet(self):
        if self.operation.get() in [1,2]:
            return self.inputMultiSet[0].get()
        else:
            return self.inputSet.get()

    def addItem(self, dic, itemId, item, allowDup=None):
        allowDup = not self.removeDuplicates.get() if allowDup == None else allowDup

        if itemId in dic and allowDup:
            dic[itemId] += [item.clone()]
        else:
            dic[itemId] = [item.clone()]
        return dic

    def getRefValue(self, idx=1):
        referenceValue = self.filterValue.get() if idx == 1 else self.filterValue2.get()
        value = self.inputSet.get().getFirstItem().getAttributeValue(self.filterColumn.get())

        if isinstance(value, float):
            referenceValue = float(referenceValue)
        elif isinstance(value, int):
            referenceValue = int(referenceValue)
        return referenceValue

    def parseTopRankParam(self):
        descending = True
        inputSet, threshold = self.inputSet.get(), self.threshold.get().strip()
        # Filed with % at the end
        if threshold.endswith('%'):
            perc = float(threshold[:-1])
            finalNumber = round(perc * len(inputSet) / 100)

        # percentage specified in decimal format.
        elif -1 < float(threshold) < 1:
            prop = float(threshold)
            finalNumber = round(prop * len(inputSet))

        # Else integers (positive or negative).
        else:
            finalNumber = int(float(threshold))

        # If negative
        if finalNumber < 0:
            finalNumber = abs(finalNumber)
            descending = False

        return finalNumber, descending

    def getAttrValue(self, item, opAttr):
        try:
            opId = getattr(item, opAttr).get()
        except:
            opId = getattr(item, opAttr)
        return opId
