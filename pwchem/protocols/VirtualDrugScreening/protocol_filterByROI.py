# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to filter a set of small molecules by the ROIs they are docked to
(number of  dockings per ROI, ROI ID, ...)

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os, re

ID, POPU = 0, 1

class ProtocolFilterByROI(EMProtocol):
    """
    Executes the filtering on a set of docked small molecules based on the ROIs they are
    docked to
    """
    _label = 'Filter by ROIs'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputMolecules', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Docked Molecules: ",
                       help='Select the SetOfSmallMolecules to filter')
        form.addParam('operation', params.EnumParam, default=0, label='Filter criteria: ',
                      choices=['ID', 'Population'],
                      help="Specify the filtering criteria")
        form.addParam('opNumber', params.IntParam, default=1,
                      label='Criteria number: ',
                      help="Number of the defined criteria. \n"
                           "For ID: id(s) of the kept ROIs (Comma separated: id1,id2,id3,...)\n"
                           "For population: number of most/less populated ROIs to keep")
        form.addParam('maxmin', params.BooleanParam, default=True,
                      label='Keep maximum values: ', condition='operation==1',
                      help='True to keep the maximum values. False to get the minimum')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        self.fMols = []
        if self.operation == ID:
            ids = self.opNumber.get().split(',')
            for mol in self.inputMolecules.get():
                if mol.getGridId() in ids:
                    self.fMols.append(mol.clone())

        elif self.operation == POPU:
            pocketCount = {}
            for mol in self.inputMolecules.get():
                pocketId = mol.getGridId()
                if pocketId in pocketCount:
                    pocketCount[pocketId] += 1
                else:
                    pocketCount[pocketId] = 1

            pKeys, pValues = list(pocketCount.keys()), list(pocketCount.values())
            pValues, pKeys = (list(t) for t in zip(*sorted(zip(pValues, pKeys), reverse=self.maxmin.get())))

            for mol in self.inputMolecules.get():
                if mol.getGridId() in pKeys[:self.opNumber.get()]:
                    self.fMols.append(mol.clone())



    def createOutputStep(self):
        inputProteinFile = self.inputMolecules.get().getProteinFile()
        outDocked = SetOfSmallMolecules(filename=self._getPath('outputSmallMolecules.sqlite'))
        outDocked.setDocked(True)
        outDocked.setProteinFile(inputProteinFile)
        for outDock in self.fMols:
            outDocked.append(outDock)
        self._defineOutputs(outputSmallMolecules=outDocked)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        """ Try to find warnings on define params. """
        validations = []
        return validations

    # --------------------------- UTILS functions -----------------------------------