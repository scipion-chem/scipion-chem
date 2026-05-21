#Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
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

# General imports
import os

# Scipion em imports
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import getBaseName


_chargeMethods = ['gas', 'bcc']
_atomTypes = ['gaff', 'amber', 'gaff2', 'amber2']
_qProgs = ['mopac', 'sqm', 'divcon']

def performACPYPE(protocol, molFile, outDir, kwargs):
  args = f'-i {molFile} -b {kwargs["molName"]} -c {kwargs["chargeMethod"]} ' \
         f'-m {kwargs["multip"]} -a {kwargs["atomType"]} -q {kwargs["qprog"]} '
  if 'netCharge' in kwargs:
    args += f' -n {kwargs["netCharge"]}'
  if 'outType' in kwargs:
    args += f' -o {kwargs["outType"]}'
  Plugin.runACPYPE(protocol, args=args, cwd=outDir)

class ProtocolLigandParametrization(EMProtocol):
    """
    Protocol for ligand parametrization using ACPYPE.

    This protocol prepares a selected ligand from a SetOfSmallMolecules for
    molecular dynamics simulations by generating force-field compatible files
    (e.g. GAFF/AMBER formats) through the ACPYPE pipeline.

    It acts as a bridge between small molecule libraries and MD-ready ligand
    parameters.

    Overview
    --------
    The protocol performs the following steps:

    1. Input selection
       - Takes a SetOfSmallMolecules as input
       - Selects a specific ligand by name

    2. Parameter setup
       - Defines charge model, atom type, multiplicity, and quantum method
       - Optionally sets net charge

    3. Execution of ACPYPE
       - Builds command-line arguments
       - Runs ACPYPE via pwchem Plugin

    4. Output extraction
       - Retrieves generated .mol2 parametrized ligand
       - Updates molecule object inside a new set

    Input
    -----
    inputSmallMolecules:
        SetOfSmallMolecules containing candidate ligands.

    inputLigand:
        Name of the ligand to be parametrized (must match molecule identifier).

    chargeMethod:
        Method used to compute atomic charges:
        - gas
        - bcc (AM1-BCC)

    netCharge:
        Optional total charge of the ligand (if not provided, guessed automatically).

    multip:
        Spin multiplicity (2S+1), default = 1.

    atomType:
        Force field atom type definition:
        - gaff
        - amber
        - gaff2 (default)
        - amber2

    qprog:
        Quantum chemistry backend used for charge calculation:
        - mopac
        - sqm (default)
        - divcon

    Workflow
    --------
    1. Identify ligand inside input set:
       - Matches molecule by string comparison
       - Clones selected molecule

    2. Build ACPYPE command:
       - Input molecule file
       - Ligand name
       - Charge method
       - Multiplicity
       - Atom type
       - Quantum program
       - Optional net charge and output type

    3. Run ACPYPE:
       - Executed via pwchem Plugin wrapper
       - Outputs force-field parameter files in .acpype directory

    4. Locate output:
       - Finds generated .mol2 file
       - Updates molecule file path

    Output
    ------
    outputSmallMolecules:
        New SetOfSmallMolecules containing:
        - original molecules
        - updated ligand with parametrized .mol2 file

    Internal utilities
    ------------------
    getSpecifiedMol():
        Selects and clones the ligand matching inputLigand name.

    getOutDir():
        Finds ACPYPE output directory (*.acpype).

    getOutMol2():
        Retrieves generated .mol2 file from ACPYPE output folder.

    performACPYPE():
        Builds and executes ACPYPE command with selected parameters.

    Summary
    -------
    This protocol automates ligand preparation for MD simulations by wrapping
    ACPYPE execution inside Scipion, ensuring:
    - consistent ligand selection
    - reproducible parametrization
    - integration with molecular workflows

    Key features:
    - ACPYPE integration
    - Flexible charge and force field options
    - Automatic ligand extraction from molecule sets
    - MD-ready output generation (.mol2)
    """
    _label = 'Ligand parametrization'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        iGroup = form.addGroup('Input')
        iGroup.addParam('inputSmallMolecules', params.PointerParam,
                        pointerClass='SetOfSmallMolecules', allowsNull=False,
                        label="Input  Small Molecules: ",
                        help='Select the molecules to be filtered')
        iGroup.addParam('inputLigand', params.StringParam,
                        label='Ligand to prepare: ',
                        help='Specific ligand to prepare in the system')

        self._defineACPYPEparams(form)


    def _defineACPYPEparams(self, form, condition=True):
      pGroup = form.addGroup('Ligand parametrization', condition=condition)
      pGroup.addParam('chargeMethod', params.EnumParam, default=1, label='Charge methods: ',
                      choices=_chargeMethods, condition=condition,
                      help="Method used to calculate charges.\n")
      pGroup.addParam('netCharge', params.StringParam, default='', label='Net charge: ',
                      expertLevel=params.LEVEL_ADVANCED, condition=condition,
                      help="Net molecular charge (int), it tries to guess it if not not declared")
      pGroup.addParam('multip', params.StringParam, default='1', label='Multiplicity: ',
                      expertLevel=params.LEVEL_ADVANCED, condition=condition,
                      help="Multiplicity (2S+1), default is 1")

      pGroup.addParam('atomType', params.EnumParam, default=2, label='Charge methods: ', choices=_atomTypes,
                      condition=condition,
                      help="Atom type, can be gaff, gaff2 (default), amber (AMBER14SB) or amber2 (AMBER14SB + GAFF2)")
      pGroup.addParam('qprog', params.EnumParam, default=1, label='Quantum program: ', choices=_qProgs,
                      condition=condition,
                      help="Quantum program to be used fro parametrization\n"
                           "am1-bcc flag, sqm (default), divcon, mopac")
      return pGroup

    
    def getParameters(self):
      parDic = {}
      enums, ints = ['chargeMethod', 'atomType', 'qprog'], ['netCharge', 'multip']
      for enum in enums:
        parDic[enum] = self.getEnumText(enum)
      for i in ints:
        val = getattr(self, i).get()
        if val:
          parDic[i] = val
      return parDic
    
    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.parametrizeStep)
        self._insertFunctionStep(self.createOutputStep)

    def parametrizeStep(self):
        mol = self.getSpecifiedMol()
        molFile = os.path.abspath(mol.getFileName())
        outDir = os.path.abspath(self._getExtraPath())
        kwargs = self.getParameters()
        kwargs['molName'] = mol.getMolName()
        
        performACPYPE(self, molFile, outDir, kwargs)


    def createOutputStep(self):
        outDir = self.getOutDir()
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        mol = self.getSpecifiedMol()

        outFile = self.getOutMol2(outDir)
        mol.setFileName(outFile)
        newMols.append(mol)

        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- UTILS functions -----------------------------------
    def getSpecifiedMol(self):
        myMol = None
        for mol in self.inputSmallMolecules.get():
          if mol.__str__() == self.inputLigand.get():
            myMol = mol.clone()
            break
        if myMol == None:
            print('The input ligand is not found')
            return None
        else:
            return myMol

    def getOutDir(self):
      for d in os.listdir(self._getExtraPath()):
        d = self._getExtraPath(d)
        if os.path.isdir(d) and d.endswith('.acpype'):
          return d
      return None

    def getOutMol2(self, outDir):
      for d in os.listdir(outDir):
        d = os.path.join(outDir, d)
        if d.endswith('.mol2'):
          return d
      return None

