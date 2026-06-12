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
This protocol is used to draw a set of small molecules using the JChemPaint software

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem import Plugin as pwchem_plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
import os


class ProtDrawMolecules(EMProtocol):
    """
    AI Generated:

    This protocol is used to create and manage small molecule structures by drawing
    them interactively using JChemPaint and integrating the resulting molecules
    into a Scipion-compatible workflow.

    It provides a simple graphical entry point for users to manually design or
    sketch chemical compounds, which are then automatically converted into
    standard molecular formats and stored as a SetOfSmallMolecules object.

    Core Concepts
    -------------
    JChemPaint Integration:
        External molecular editor used to draw and edit 2D chemical structures.
        The protocol launches JChemPaint as an external application.

    Molecular Conversion:
        Drawn structures are exported from JChemPaint and converted into
        standard chemical formats (e.g., MOL2) using OpenBabel.

    SetOfSmallMolecules:
        Scipion data structure used to store and manage collections of small
        molecules, each represented as a SmallMolecule object.

    Workflow
    --------
    1. Launch JChemPaint GUI for user-driven molecular drawing.
    2. User creates and saves one or more molecules in the protocol directory.
    3. The protocol scans the output directory for generated files.
    4. Each file is converted into MOL2 format using OpenBabel.
    5. Converted molecules are wrapped into SmallMolecule objects.
    6. A SetOfSmallMolecules is created and populated.
    7. The resulting dataset is exported for downstream use in workflows.

    File Handling
    -------------
    - Input:
        No structured input is required beyond launching the drawing interface.

    - Output Directory:
        All generated files are stored in the protocol working directory
        (including temporary and extra folders).

    - Conversion:
        Non-MOL2 files are converted using OpenBabel before inclusion.

    Output
    ------
    - outputStructROIs:
        A SetOfSmallMolecules containing all molecules drawn in JChemPaint
        and successfully converted.

    Each SmallMolecule includes:
        - Molecular structure file (MOL2)
        - Assigned molecule name derived from filename

    Use Cases
    ---------
    - Manual design of ligands for docking experiments
    - Rapid prototyping of chemical structures
    - Preparation of custom compound libraries
    - Integration of user-designed molecules into computational workflows
    """
    _label = 'Draw small molecules'
    molNumber = 1

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('launchJChem', params.LabelParam,
                      label='Draw new molecules',
                      help='Launches JChemPaint to draw a new molecule')


    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('launchJChemPaint')
        self._insertFunctionStep('createOutput')


    def launchJChemPaint(self):
        pwchem_plugin.runJChemPaint(self, cwd=self._getExtraPath())

    def createOutput(self):
        cFiles = []
        for file in os.listdir(self._getPath()):
            if not file in ['logs', 'tmp', 'extra'] or file.endswith('.sqlite'):
                try:
                    fileName, ext = os.path.splitext(os.path.basename(file))
                    newFile = 'extra/' + fileName + '.mol2'
                    args = ' {} -O {}'.format(file, newFile)
                    pwchem_plugin.runOPENBABEL(self, args=args, cwd=self._getPath())
                    cFiles.append(self._getPath(newFile))
                except:
                    pass

        if len(cFiles) > 0:
            outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
            for molFile in cFiles:
                newMol = SmallMolecule(smallMolFilename=molFile)
                newMol.setMolName(os.path.splitext(os.path.basename(molFile))[0])
                outputSet.append(newMol)

            self._defineOutputs(outputStructROIs=outputSet)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = ['If the created molecules are saved in: \n\n{}\nand in .mol format\n they will be automatically '
                   'incorporated in the Scipion workflow.\nElse, you will need to import the files you generate and '
                   'save somewhere else with the import small molecules '
                   'protocol'.format(os.path.abspath(self._getPath()))]
        return summary

    def _warnings(self):
        ws = ['If the created molecules are saved in: {} and in .mol format, they will be automatically incorporated '
              'in the Scipion workflow.\n'
              'Else, you will need to import the files you generate and save somewhere else with '
              'the import small molecules protocol'.format(os.path.abspath(self._getPath()))]
        return ws

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------






