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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to obtain peptides from sequences.

"""
import pwchem.constants as constants
import os

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwchem import ESM_DIC, Plugin


class ProtPeptideFromSequence(EMProtocol):
    """
    Generates a pdb peptide file from a sequence.
    """
    _label = 'Get peptide from sequence'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        #todo start with one then do multiple to get a SetOfAtomStructs
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('name', params.StringParam, default='',
                       label='Input sequence name: ',
                       help='Input sequence name.')
        group.addParam('inputSeq', params.StringParam, default='',
                       label='Input sequence: ',
                       help='Input sequence to get peptide pdb.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        #self._insertFunctionStep('createEnvStep')
        self._insertFunctionStep('createOutputStep')

    def createEnvStep(self):
        Plugin.addESMPackage(constants.ESM_DIC)

    def createOutputStep(self):

        pdbFile = self.getESMpdb()
        print(pdbFile)


    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions ------------------------
    def getESMpdb(self):
        """
        Retrieve a predicted peptide PDB structure from ESM Atlas web endpoint.
        """
        sequence = self.inputSeq.get()
        sequence = sequence.strip().upper()
        print(sequence)

        esmHome = Plugin.getVar(ESM_DIC["home"])

        outputDir = os.path.abspath(self._getExtraPath())
        fastaPath = os.path.abspath(self._getTmpPath('sequence.fasta'))

        header = ">unnamed"
        with open(fastaPath, "w") as f:
            f.write(f"{header}\n{sequence}\n")
        pdbName = "unnamed.pdb"
        pdbPath = os.path.join(outputDir,pdbName)
        scriptName = "fold.py"
        script_dir = os.path.join(esmHome, "scripts")

        args = f"-i {fastaPath} -o {outputDir}"
        print(args)

        self.info(f"Running ESMFold for sequence: {sequence}")

        Plugin.runScript(protocol=self,
                       scriptName=scriptName,
                       args=args,
                       env=ESM_DIC,
                       cwd=str(outputDir),
                       popen=False,
                       scriptDir=script_dir)

        if not pdbPath.exists():
            raise FileNotFoundError(f"ESMFold did not produce a PDB file at {pdbPath}")

        self.info(f"ESMFold completed successfully: {pdbPath}")
        return str(pdbPath)

