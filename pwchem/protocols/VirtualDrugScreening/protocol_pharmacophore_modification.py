# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# # # *
# # # *
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import json
from os.path import abspath

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import PharmFeature, PharmacophoreChem
from pwchem.constants import *
from pwchem.utils import *


ADD, REM, MOD = 'ADD', 'REMOVE', 'MODIFY'

class ProtocolPharmacophoreModification(EMProtocol):
    """
    Perform the modification of a pharmacophore.
    Allows to add, remove or modify its chemical features
    """
    _label = 'Pharmacophore modification'

    ##### -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form.addParam('inputPharmacophores', params.MultiPointerParam,
                      pointerClass='PharmacophoreChem', allowsNull=True,
                      label="Input pharmacophore: ",
                      help='Select the pharmacophore to modify.\n'
                           'If empty, you can only add features (which will become a new pharmacophore)')
        form.addParam('inputAtomStruct', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input pharmacophore related protein: ",
                      help='Select the protein structure where the pharmacophore will be placed.\n'
                           'If empty, the input pharmacophore related protein will be kept as it is.')

        group = form.addGroup('Define operation')
        group.addParam('operation', params.EnumParam, label='Feature operation: ',
                       choices=[ADD, REM, MOD], default=0,
                       help="Chose whether to add, remove or modify a feature for the pharmacophore."
                            "Once defined, the operation must be stored using the wizard so several modifications "
                            "can be defined")

        group.addParam('currentFeatures', params.StringParam, label='Select input feature: ',
                       condition='operation!=0', default='',
                       help="Chose the input present feature to remove / modify using the wizard")

        group.addParam('showCurrent', params.LabelParam, label='Show current attributes: ',
                       condition='operation==2',
                       help="Display the current attributes of the selected feature")

        group.addParam('featType', params.EnumParam, label='Feature type: ', condition='operation!=1',
                       choices=FEATURE_LABELS_SIMPLE + FEATURE_LABELS_ADVANCED, default=0,
                       help="Chose the type of feature to add.")

        line = group.addLine('Coordinates: ', condition='operation!=1',
                             help='Coordinates of the feature too add')
        line.addParam('featX', params.FloatParam, label='X: ', condition='operation!=1', default=0.0)
        line.addParam('featY', params.FloatParam, label='Y: ', condition='operation!=1', default=0.0)
        line.addParam('featZ', params.FloatParam, label='Z: ', condition='operation!=1', default=0.0)

        group.addParam('featRadius', params.FloatParam, label='Radius: ', condition='operation!=1',
                       default=1.0, help="Choose the radius for the feature to add")

        group.addParam('addOperation', params.LabelParam, label='Add operation: ',
                       help="Click on the wizard to save the defined operation in the list of operations")

        group = form.addGroup('List of operations')
        group.addParam('operationList', params.TextParam, width=100, height=10,
                       default='', label='List of operations:',
                       help='Defines the list of operations that will be performed on the pharmacophore.'
                            'If several operations are defined for the same feature, only the last one '
                            'will be performed.\nManual modification of this parameter is strongly not recommended '
                            'except for deleting entire lines.')


        # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')

    def createOutputStep(self):
        if len(self.inputPharmacophores) > 0:
            outPharm = PharmacophoreChem.createCopy(self.inputPharmacophores[0].get(), self._getPath(),
                                                    copyInfo=True, copyItems=False)
        else:
            outPharm = PharmacophoreChem().create(outputPath=self._getPath())

        if self.inputAtomStruct.get():
            outPharm.setProteinFile(self.inputAtomStruct.get().getFileName())

        operDic, addList = self.getOperationDic()
        featsDic = self.getCurrentFeaturesDic()

        for featId in featsDic:
            if featId in operDic:
                jDic = operDic[featId]
                if jDic:
                    pharmFeat = self.buildFeatureFromJDic(jDic)
                    outPharm.append(pharmFeat)
            else:
                pharmFeat = featsDic[featId]
                pharmFeat.cleanObjId()
                outPharm.append(pharmFeat)

        for jDic in addList:
            pharmFeat = self.buildFeatureFromJDic(jDic)
            outPharm.append(pharmFeat)

        self._defineOutputs(outputPharmacophore=outPharm)

    # --------------------------- INFO functions -----------------------------------

    def _validate(self):
        vals = []
        operDic, _ = self.getOperationDic()
        if len(operDic) > 0 and len(self.inputPharmacophores) == 0:
            vals.append('If no input pharmacophore is specified, you can only add new features.\n'
                        'No modification or deletion can be performed (over which features would it be?)')
        return vals

    # --------------------------- UTILS functions -----------------------------------

    def getPresentFeatures(self):
        if hasattr(self, 'inputPharmacophores') and len(self.inputPharmacophores) > 0:
            featList = []
            for pharm in self.inputPharmacophores:
              for feat in pharm.get():
                  featList.append(str(feat))
        else:
            return []

    def getOperationDic(self):
        operDic, addList = {}, []
        for operation in self.operationList.get().split('\n'):
            if operation.strip():
              oper = operation.split('|')[0]
              if oper.strip() == ADD:
                  rest = operation.split('|')[1]
                  feat = json.loads(rest.strip())
                  addList += [feat]

              elif oper.strip() == REM:
                  setStr, rest = operation.split('|')[1:]
                  setId = int(setStr.split()[1])
                  featId = int(rest.split()[1])
                  operDic[(setId, featId)] = ''

              elif oper.strip() == MOD:
                  setStr, rest = operation.split('|')[1:]
                  setId = int(setStr.split()[1])
                  featId = int(rest.split()[1])
                  feat = json.loads(rest.split('TO:')[1].strip())
                  operDic[(setId, featId)] = feat
        return operDic, addList

    def getCurrentFeaturesDic(self):
        featDic = {}
        if len(self.inputPharmacophores) > 0:
          for i, pharm in enumerate(self.inputPharmacophores):
            for feat in pharm.get():
                nFeat = feat.clone()
                featDic[(i, nFeat.getObjId())] = nFeat
        return featDic

    def buildFeatureFromJDic(self, jDic):
        loc = eval(jDic['Coords'])
        pharmFeat = PharmFeature(type=jDic['Type'], radius=jDic['Radius'],
                               x=loc[0], y=loc[1], z=loc[2])
        pharmFeat.cleanObjId()
        return pharmFeat
