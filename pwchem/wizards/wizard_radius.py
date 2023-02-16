# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""
This wizard will show the structure of the pdb using a matplotlib viewer 
to select the radius of the sphere that contains the protein or a desired zone.
"""

# Imports
from pyworkflow.object import Pointer, PointerList
from pwem.objects import AtomStruct
from pwem.wizards.wizard import EmWizard, VariableWizard
from pwem.wizards.wizards_3d.mask_structure_wizard import MaskStructureWizard
from pwchem.utils import getProteinMaxDiameter

from pwchem.objects import SetOfSmallMolecules

class GetRadiusProtein(VariableWizard):
    _targets, _inputs, _outputs = [], {}, {}

    def getASFile(self, inputObj):
        if issubclass(type(inputObj), AtomStruct):
            return inputObj.getFileName()
        elif issubclass(type(inputObj), SetOfSmallMolecules):
            return inputObj.getProteinFile()

    def getInputObj(self, protocol, inputParams):
        inputPointer = getattr(protocol, inputParams[0])
        if type(inputPointer) == PointerList:
            inpObj = inputPointer[0].get()
        elif type(inputPointer) == Pointer:
            inpObj = inputPointer.get()
        return inpObj

    def show(self, form):
        inputParams, outputParam = self.getInputOutput(form)
        protocol = form.protocol

        inputObj = self.getInputObj(protocol, inputParams)
        if not inputObj:
            print('You must specify input structure')
            return

        strFile = self.getASFile(inputObj)
        try:
            plt = MaskStructureWizard(strFile)
            plt.initializePlot()
            form.setVar(outputParam[0], plt.radius)
        except:
            print('Cannot open the file with MaskStructureWizard.'
                  'Trying to assign the radius directly')
            radius = getProteinMaxDiameter(strFile) / 2
            form.setVar(outputParam[0], int(radius + 1))


