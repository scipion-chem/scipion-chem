# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import re

import os

import pyworkflow.object as pwobj
from pwem.wizards import VariableWizard
from pyworkflow.gui import ListTreeProviderString
from pyworkflow.gui import dialog

from pwchem.protocols.VirtualDrugScreening.protocol_consensus_structROIs import ProtocolConsensusStructROIs


class ExtractObjectsWizard(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getAllItems(self, form, inputParamName):
        protocol = form.protocol
        allItems = []

        multiPointer = getattr(protocol, inputParamName)

        for roiSetPointer in multiPointer:
            roiSet = roiSetPointer.get()
            allItems.append(roiSet)

        return allItems

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)

        allSets = self.getAllItems(form, inputParam[0])

        displayList = [pwobj.String(str(s)) for s in allSets]

        provider = ListTreeProviderString(displayList)
        dlg = dialog.ListDialog(
            form.root,
            "Input ROI Sets",
            provider,
            "Select one input set"
        )

        if dlg.values:
            selectedStr = dlg.values[0].get()
            form.setVar(outputParam[0], selectedStr)


ExtractObjectsWizard().addTarget(protocol=ProtocolConsensusStructROIs,
                             targets=['fromSet'],
                             inputs=['inputStructROIsSets'],
                             outputs=['fromSet'])