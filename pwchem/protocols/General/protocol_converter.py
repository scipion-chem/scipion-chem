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


import os, shutil

from pyworkflow.protocol.params import PointerParam, EnumParam
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toCIF, toPdb

from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel

RDKIT, OBABEL = 'RDKit', 'OpenBabel'
extDic = {'PDB': '.pdb', 'cif': '.cif', 'Mol2': '.mol2', 'SDF': '.sdf', 'Smiles': '.smi'}

class ConvertStructures(EMProtocol):
    """
    Convert a set of input ligands or a protein structure to a specific file format
    """

    _label = 'Convert format'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputType', EnumParam, default=0,
                       choices=["Small molecules", 'Target structure'],
                       label='Input type')

        form.addParam('inputSmallMols', PointerParam, pointerClass="SetOfSmallMolecules", condition='inputType==0',
                      label='Set of small molecules:', allowsNull=False,
                      help="The allowed format are pdb, cif, mol2, sdf and smi")

        form.addParam('outputFormatSmall', EnumParam, default=2, condition='inputType==0',
                       choices=['PDB', 'Mol2', 'SDF', 'Smiles'],
                       label='Output format',
                       help="Output format for the converted molecules")

        form.addParam('useManager', EnumParam, default=0, label='Convert using: ',
                      condition='inputType==0', choices=[RDKIT, OBABEL],
                      help='Whether to convert the input molecules using RDKit or OpenBabel')

        form.addParam('inputStructure', PointerParam, pointerClass= "AtomStruct",
                      condition='inputType==1', label='Input structure:', allowsNull=False,
                      help="The allowed format are pdb or cif")

        form.addParam('outputFormatTarget', EnumParam, default=0,
                      condition='inputType==1', choices=['PDB', 'cif'], label='Output format')



    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')

    def convertStep(self):

        if self.inputType==0:  # Small molecules
            outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')

            self.convErrors = []  # Save the file paths that could not be transformed
            for mol in self.inputSmallMols.get():
                fnSmall = os.path.abspath(mol.smallMoleculeFile.get())
                fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]

                outFormat = extDic[self.getEnumText('outputFormatSmall')]
                outDir = os.path.abspath(self._getExtraPath())
                fnOut = os.path.join(outDir, fnRoot + outFormat)

                args = ' -i "{}" -of {} -o {} --outputDir {}'.format(fnSmall, outFormat, fnOut, outDir)
                if self.getEnumText('useManager') == OBABEL:
                    Plugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)
                else:
                    Plugin.runScript(self, 'rdkit_IO.py', args, env='rdkit', cwd=outDir)

                if os.path.exists(fnOut):
                    smallMolecule = SmallMolecule(smallMolFilename=fnOut, molName='guess')
                    outputSmallMolecules.append(smallMolecule)
                else:
                    self.convErrors.append(fnRoot)

            if len(outputSmallMolecules) > 0:
                self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
                self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)

            if len(self.convErrors) > 0:
                print("The following entries could not be converted: %s" % ','.join(self.convErrors))

        else:
            fnStructure = os.path.abspath(self.inputStructure.get().getFileName())
            args = self.inputArg(os.path.abspath(fnStructure))
            fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]

            outFormat = extDic[self.getEnumText('outputFormatTarget')]
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.join(outDir, fnRoot + outFormat)

            if outFormat == '.pdb':
                convFn = toPdb(fnStructure, fnOut)
            elif outFormat == '.cif':
                convFn = toCIF(fnStructure, fnOut)

            if convFn == fnStructure:
                shutil.copy(convFn, fnOut)

            if os.path.exists(fnOut):
                target = AtomStruct(filename=fnOut)
                self._defineOutputs(outputStructure=target)
                self._defineSourceRelation(self.inputStructure, target)

    def inputArg(self, fn):  # Input format file (fn)

        if fn.endswith('.pdb'):  # Protein Data Bank
            args = "-ipdb"
        elif fn.endswith('.cif'):  # cif (crystallography information)
            args = "-icif"

        elif fn.endswith('.sdf'):  # MDL MOL FORMAT
            args = "-isdf"
        elif fn.endswith('.sd'):  # MDL MOL FORMAT
            args = "-isd"

        elif fn.endswith('.mol2'):  # Sybyl Mol2 format (3D)
            args = "-imol2"
        elif fn.endswith('.smi') or fn.endswith('.smiles'):  # Smiles format (2D)
            args = "-ismi"
        else:
            error = " Input format was not recognize (The allowed format are pdb, cif, mol2, sdf and smi)"
            raise Exception(error)
        return args + " %s" % os.path.abspath(fn)

    # --------------------------- Summary functions --------------------
    def _summary(self):
        summary=[]

        if self.inputType.get()==0:
            if self.outputFormatSmall.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatSmall.get() == 2:
                summary.append('Converted to Mol2')
            elif self.outputFormatSmall.get() == 3:
                summary.append('Converted to SDF')
            elif self.outputFormatSmall.get() == 4:
                summary.append('Converted to Smiles')

        else:
            if self.outputFormatTarget.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatTarget.get() == 2:
                summary.append('Converted to Mol2')
        return summary