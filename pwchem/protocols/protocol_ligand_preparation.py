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
Protocol Steps:
0. Input: Get the ligands through the ZINC protocol or import molecules (setofSmallMolecules)
1. Convert to mol2 and pdb format for images (openbabel)
2. Add hydrogens (obabel -h) and remove water
    2.1 Assign charges (Add gasteiger charges (ADT or computeGasteigerCharges de RDKit or babel --partialcharges mmff94 or gasteiger)
3. Generate low energy conformers (openbabel with Confab or RDKIT AllChem.EmbedMolecule)
"""

from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.protocol.params import LEVEL_ADVANCED
from pyworkflow.utils import Message
import pyworkflow.object as pwobj

import os, re, glob, shutil

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, splitConformerFile, appendToConformersFile, relabelAtomsMol2, \
  splitPDBLine, natural_sort


class ProtChemOBabelPrepareLigands(EMProtocol):
    """
    Prepare a set of molecules for use in a docking program (for example, Rosetta DARC).
    Sets the partial atomic charges, generates low-energy conformers. This is done with OpenBabel.
    """

    _label = 'OBabel Ligand preparation'
    _dic_method = {0: "gasteiger",
                  1: "mmff94",
                  2: "qeq",
                  3: "qtpie",
                  4: "eqeq",
                  5: "eem",
                  6: "none"}


    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMols', params.PointerParam, pointerClass="SetOfSmallMolecules",
                      label='Set of small molecules:', allowsNull=False,
                      help='It must be in pdb or mol2 format, you may use the converter')

        H_C = form.addGroup("Charges assignation")
        H_C.addParam('method_charges', params.EnumParam,
                      choices=list(self._dic_method.values()),
                      default=0,
                      label='Assign charges method', allowsNull=False,
                      help='Choose the method to add partial charges to small molecule atoms. OpenBabel is used.')

        H_C.addParam('ph', params.BooleanParam,
                      default=False,
                      label='Do you want to set pH of the ligand?:',  # Next step --> PRopKA
                      help='Set the pH of the ligand environment using OpenBabel program. \n'
                           'In this way, the ligand will present more or less hydrogens depending on the pH. '
                           'Note that openbabel (program that add H) has not a general model of pH-dependent '
                           'protonation. It has a set of rules in a file called phmodel.txt,'
                           ' particularly for amino acids. \n\n For more accurate pKa determination and H addition'
                           ', you will probably want semiempirical quantum calculations or a more complete model '
                           '(It is not trivial).',
                      expertLevel=LEVEL_ADVANCED)

        H_C.addParam('phvalue', params.FloatParam, condition="ph",
                      default=7.4,
                      label='pH value:',
                      help='Set the pH of the ligand environment.',
                      expertLevel=LEVEL_ADVANCED)


        conformers = form.addGroup("Conformers generation")
        conformers.addParam('doConformers', params.BooleanParam, default=False,
                      label='Do you want to generate conformers? ', allowsNull=False,
                      help='You can produce conformers of the ligand in order to do a better rigid docking')
        conformers.addParam('method_conf', params.EnumParam,
                            choices=["OpenBabel Genetic Algorithm", "OpenBabel Confab"],
                            default=0, condition = "doConformers",
                            label='Method of conformers generation',
                            help='Method of conformers generation. If Confab fails due to the impossibility '
                                 'of assigning a force fields (there is a possibility that it may occur), you should'
                                 'use Genetic Algorithm generator ')


        conformers.addParam('number_conf', params.IntParam,
                            default=200, condition = "doConformers",
                            label='Max. number of conformers:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('rmsd_cutoff', params.FloatParam, condition = "method_conf != 0 and doConformers",
                            default=0.5,
                            label='RMSD cutoff:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.',
                            expertLevel=LEVEL_ADVANCED)



    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('addcharges')
        if self.doConformers.get():
            self._insertFunctionStep('conformer_generation')
        self._insertFunctionStep('createOutput')

    def addcharges(self):
        """ Assign the charges using a method available
            in the open-access and free program openbabel
        """
        inMols = self.inputSmallMols.get()
        for mol in inMols:
            fnSmall = mol.getFileName()
            fnMol = os.path.split(fnSmall)[1]      # Name of complete file
            fnRoot, fnFormat = os.path.splitext(fnMol)    # Molecule name: ID, format

            fnSmall = self.reorderAtoms(fnSmall, self._getExtraPath('{}_ordered{}'.format(fnRoot, fnFormat)))

            # 1. Add all hydrogens or add hydrogens depending on the desirable pH with babel (-p)
            # 2. Add and calculate partial charges with different methods
            index_method = self.method_charges.get()
            cmethod = self._dic_method[index_method]

            # With a given pH
            oFile = "{}_prep.mol2".format(fnRoot)
            if self.ph.get():
                args = " -i%s %s -p %s --partialcharge %s -O %s" % (fnFormat[1:], os.path.abspath(fnSmall),
                                                                    str(self.phvalue.get()), cmethod, oFile)

                runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getExtraPath()))

            else:
                args = " -i%s %s -h --partialcharge %s -O %s" % (fnFormat[1:], os.path.abspath(fnSmall), cmethod, oFile)
                runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getExtraPath()))

            oFile = relabelAtomsMol2(os.path.abspath(self._getExtraPath(oFile)))

    def conformer_generation(self):
        """ Generate a number of conformers of the same small molecule in mol2 format with
            openbabel using two different algorithm
        """

        # 3. Generate mol2 conformers file for each molecule with OpenBabel

        for file in glob.glob(self._getExtraPath("*_prep.mol2")):
            fnRoot = os.path.splitext(os.path.split(file)[1])[0]  # ID or filename without -prep.mol2

            if self.method_conf.get() == 0:  # Genetic algorithm
                args = " %s --conformer --nconf %s --score rmsd --writeconformers -O %s_conformers.mol2" %\
                       (os.path.abspath(file), self.number_conf.get()-1, fnRoot)
            else:  # confab
                args = " %s --confab --original --verbose --conf %s --rcutoff %s -O %s_conformers.mol2" % \
                       (os.path.abspath(file), self.number_conf.get()-1, str(self.rmsd_cutoff.get()), fnRoot)

            runOpenBabel(protocol=self, args=args, cwd=os.path.abspath(self._getExtraPath()))


    def createOutput(self):
        """Create a set of Small Molecules as output with the path to:
              - Path to small molecule with H (mol2 format)
              - Path to conformers file (mol2 format)
        """

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='')

        for mol in self.inputSmallMols.get():
            fnSmall = self._getExtraPath("{}_prep.mol2".format(mol.getMolName()))
            if self.doConformers:
                fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                outDir = self._getExtraPath(fnRoot)
                os.mkdir(outDir)
                firstConfFile = self._getTmpPath('{}-{}.mol2'.format(fnRoot, 1))
                shutil.copy(fnSmall, firstConfFile)
                confFile = self._getExtraPath("{}_conformers.mol2".format(fnRoot))
                confFile = appendToConformersFile(confFile, firstConfFile,
                                                  beginning=True)
                confDir = splitConformerFile(confFile, outDir=outDir)
                for molFile in os.listdir(confDir):
                    molFile = os.path.abspath(os.path.join(confDir, molFile))
                    newSmallMol = SmallMolecule(smallMolFilename=molFile)
                    newSmallMol._ConformersFile = pwobj.String(confFile)
                    outputSmallMolecules.append(newSmallMol)
            else:
                newSmallMol = SmallMolecule(smallMolFilename=fnSmall)
                outputSmallMolecules.append(newSmallMol)

        if outputSmallMolecules is not None:
            self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
            self._defineSourceRelation(self.inputSmallMols, outputSmallMolecules)


    # --------------------------- UTILS functions ------------------------------
    def _validate(self):
        """ Validate if the inputs are in mol2 or pdb format
        """
        errors = []
        return errors

    def getLigandCode(self, paramsFile):
      with open(paramsFile) as f:
        code = f.readline().split()[1]
      return code

    def reorderAtoms(self, inFile, outFile):
        '''Atom lines in the file are reordered so the atomNames numbers are in order (C2, C1, C3 -> C1, C2, C3)'''
        atomsDic = {}
        writeFirst, writeLast = '', ''
        with open(inFile) as fIn:
          if inFile.endswith('pdb'):
              for line in fIn:
                  sline = splitPDBLine(line)
                  if sline[0] in ['ATOM', 'HETATM']:
                      atomsDic[sline[2]] = line
                  else:
                      if atomsDic == {}:
                          writeFirst += line
                      else:
                          writeLast += line

          elif inFile.endswith('mol2'):
              atomLines = False
              for line in fIn:
                  if not atomLines and line.startswith('@<TRIPOS>ATOM'):
                      atomLines = True
                      writeFirst += line
                  elif atomLines and line.startswith('@'):
                      atomLines = False
                      writeLast += line
                  elif atomLines:
                      atomsDic[line.split()[1]] = line
                  else:
                      if atomsDic == {}:
                        writeFirst += line
                      else:
                        writeLast += line
          else:
              #Don't touch other kind of files
              return inFile

        with open(outFile, 'w') as f:
            f.write(writeFirst)
            for key in natural_sort(atomsDic.keys()):
                f.write(atomsDic[key])
            f.write(writeLast)
        return outFile






