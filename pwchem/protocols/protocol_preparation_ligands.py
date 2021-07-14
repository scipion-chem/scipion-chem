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
import glob
import os

from pyworkflow.protocol.params import PointerParam
from .protocol_preparation_receptor import ProtBioinformaticsADTPrepare
from pwchem.objects import SmallMolecule, SetOfSmallMolecules

class ProtBioinformaticsADTPrepareLigands(ProtBioinformaticsADTPrepare):
        """Prepare ligands using Autodocking Tools from MGL"""
        _label = 'ligand preparation ADT'
        _program = ""

        def _defineParams(self, form):
            self.typeRL="ligand"
            form.addSection(label='Input')
            form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules",
                          label='Set of small molecules:', allowsNull=False,
                          help='It must be in pdb or mol2 format, you may use Schrodinger convert to change it')
            ProtBioinformaticsADTPrepare._defineParamsBasic(self, form)

        def preparationStep(self):
            for mol in self.inputSmallMols.get():
                fnSmall = mol.smallMoleculeFile.get()
                fnMol = os.path.split(fnSmall)[1]
                fnRoot = os.path.splitext(fnMol)[0]
                fnOut = self._getExtraPath(fnRoot+".pdbqt")

                args = ' -v -l %s -o %s' % (fnSmall, fnOut)
                ProtBioinformaticsADTPrepare.callPrepare(self, "prepare_ligand4", args)

        def createOutput(self):
            outputSmallMolecules = SetOfSmallMolecules().create(path=self._getPath(), suffix='SmallMols')
            for fnOut in glob.glob(self._getExtraPath("*.pdbqt")):
                smallMolecule = SmallMolecule(smallMolFilename=fnOut)
                outputSmallMolecules.append(smallMolecule)
            if len(outputSmallMolecules)>0:
                self._defineOutputs(outputSmallMols=outputSmallMolecules)
                self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)
