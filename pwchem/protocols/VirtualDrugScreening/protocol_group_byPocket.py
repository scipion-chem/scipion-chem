# **************************************************************************
# *
# * Authors:     Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params

from pwchem.objects import SetOfSmallMolecules, SmallMolecule


class ProtChemGroupByAtt(EMProtocol):
    """Group small molecules by attribute."""
    _label = 'group by pocket id'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputMols', params.PointerParam, pointerClass="SetOfSmallMolecules", label='Input molecules: ',
                      help='Set of small molecules to group by pocket id.')
        form.addParam('refColumn', params.StringParam, label='Grouping column: ', default='',
                       help='Attribute for the grouping operation.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('defineOutputStep')

    # --------------------------- STEPS subfunctions ------------------------------
    def defineOutputStep(self):
        inputMols = self.inputMols.get()
        numPockets = self.getNumPockets()

        for i in (numPockets):
            outputPath = self._getPath(f'pocket{i}')
            os.makedirs(outputPath, exist_ok=True)

            outputMols = SetOfSmallMolecules().create(outputPath=outputPath)
            for mol in inputMols:
                newMol = SmallMolecule()
                newMol.copy(mol)
                if newMol.getGridId() == i:
                    outputMols.append(newMol)
            outputName = f"outputSmallMolecules{i}"
            outputMols.setDocked()
            outputMols.proteinFile.set(inputMols.getProteinFile())

            self._defineOutputs(**{outputName: outputMols})
            self._defineSourceRelation(inputMols, outputMols)



    # --------------------------- INFO functions -----------------------------------

    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        molSet = self.inputMols.get()
        if not molSet.isDocked():
            validations += ['{} is not docked yet'.format(molSet)]

        return validations

    # --------------------------- UTILS functions -----------------------------------
    def getNumPockets(self): #todo change and get the number of indiv values of the att
        uniquePockets = {mol.getGridId() for mol in self.inputMols.get()}
        return sorted(uniquePockets)
