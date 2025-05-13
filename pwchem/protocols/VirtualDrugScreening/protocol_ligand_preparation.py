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

from pyworkflow.protocol import params
from pyworkflow.protocol.params import LEVEL_ADVANCED
import pyworkflow.object as pwobj

import os, glob

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import runOpenBabel, splitConformerFile, relabelAtomsMol2, natural_sort, getBaseName
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_filter import ProtocolBaseLibraryToSetOfMols


class ProtChemOBabelPrepareLigands(ProtocolBaseLibraryToSetOfMols):
    """
    Prepare a set of molecules for use in a docking program.
    Sets the partial atomic charges, generates low-energy conformers. This is done with OpenBabel.
    """

    _label = 'OBabel Ligand preparation'
    _dic_method = {0: "gasteiger", 1: "mmff94", 2: "qeq", 3: "qtpie", 4: "eqeq", 5: "eem", 6: "none"}
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label='Input')
        form = self.addInputParams(form)

        H_C = form.addGroup("Charges assignation")
        H_C.addParam('method_charges', params.EnumParam,
                     choices=list(self._dic_method.values()),
                     default=0, label='Assign charges method', allowsNull=False,
                     help='Choose the method to add partial charges to small molecule atoms. OpenBabel is used.')

        H_C.addParam('ph', params.BooleanParam, expertLevel=LEVEL_ADVANCED,
                     default=False, label='Do you want to set pH of the ligand?:',  # Next step --> PRopKA
                     help='Set the pH of the ligand environment using OpenBabel program. \n'
                          'In this way, the ligand will present more or less hydrogens depending on the pH. '
                          'Note that openbabel (program that add H) has not a general model of pH-dependent '
                          'protonation. It has a set of rules in a file called phmodel.txt,'
                          ' particularly for amino acids. \n\n For more accurate pKa determination and H addition'
                          ', you will probably want semiempirical quantum calculations or a more complete model '
                          '(It is not trivial).')

        H_C.addParam('phvalue', params.FloatParam, condition="ph",
                     default=7.4, expertLevel=LEVEL_ADVANCED, label='pH value:',
                     help='Set the pH of the ligand environment.')


        conformers = form.addGroup("Conformers generation")
        conformers.addParam('doConformers', params.BooleanParam, default=False,
                            label='Do you want to generate conformers? ', allowsNull=False,
                            help='You can produce conformers of the ligand in order to do a better rigid docking')
        conformers.addParam('method_conf', params.EnumParam,
                            choices=["OpenBabel Genetic Algorithm", "OpenBabel Confab"],
                            default=0, condition = "doConformers", label='Method of conformers generation: ',
                            help='Method of conformers generation. If Confab fails due to the impossibility '
                                 'of assigning a force fields (there is a possibility that it may occur), you should'
                                 'use Genetic Algorithm generator ')


        conformers.addParam('number_conf', params.IntParam,
                            default=10, condition = "doConformers", label='Max. number of conformers: ',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')

        conformers.addParam('rmsd_cutoff', params.FloatParam, condition="method_conf != 0 and doConformers",
                            default=0.5, label='RMSD cutoff: ',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.',
                            expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        aSteps, cSteps = [], []
        nt = self.numberOfThreads.get()
        os.makedirs(self.getPrepDir())
        if nt <= 1: nt = 2

        iStep = self._insertFunctionStep(self.createInputStep, nt-1, prerequisites=[])
        for it in range(nt-1):
            aSteps += [self._insertFunctionStep(self.addChargesStep, it, prerequisites=[iStep])]

        if self.doConformers.get():
            for it in range(nt-1):
                cSteps += [self._insertFunctionStep(self.conformersStep, it, prerequisites=aSteps)]
            aSteps = cSteps
        self._insertFunctionStep(self.createOutput, prerequisites=aSteps)

    def addChargesStep(self, it):
        """ Assign the charges using a method available
            in the open-access and free program openbabel
        """
        failedMols = []
        for fnSmall in self.getInputMolFiles(it):
            if isinstance(fnSmall, SmallMolecule):
              fnSmall = fnSmall.getFileName()
            fnRoot = getBaseName(fnSmall)
            _, fnFormat = os.path.splitext(fnSmall)

            try:
                fnSmall = self.reorderAtoms(fnSmall, self._getExtraPath('{}_ordered{}'.format(fnRoot, fnFormat)))
            except:
                print('Atom reordering could not be performed in {}, this could lead to errors in ahead protocols'
                      .format(fnRoot))

            # 1. Add all hydrogens or add hydrogens depending on the desirable pH with babel (-p)
            # 2. Add and calculate partial charges with different methods
            index_method = self.method_charges.get()
            cmethod = self._dic_method[index_method]

            # With a given pH
            oFile = self.getPrepDir(f"{fnRoot}.mol2")
            if self.ph.get():
                args = " -i%s '%s' -p %s --partialcharge %s -O '%s' " % (fnFormat[1:], os.path.abspath(fnSmall),
                                                                    str(self.phvalue.get()), cmethod, oFile)
            else:
                args = " -i%s '%s' -h --partialcharge %s -O '%s' " % (fnFormat[1:], os.path.abspath(fnSmall), cmethod, oFile)

            if fnFormat == '.smi':
              args += '--gen3D '

            try:
              runOpenBabel(protocol=self, args=args, popen=True)
            except:
              failedMols.append(fnRoot)

            oFile = relabelAtomsMol2(self.getPrepDir(oFile), it)

        if len(failedMols) > 0:
          with open(os.path.abspath(self._getExtraPath(f'failedCharges_{it}.txt')), 'w') as f:
            for molFn in failedMols:
              f.write(molFn + '\n')

        os.remove(self.getInputFile(it))

    def conformersStep(self, it):
        """ Generate a number of conformers of the same small molecule in mol2 format with
            openbabel using two different algorithm
        """
        failedMols = []
        for molFn in os.listdir(self.getPrepDir()):
          fnSmall = self.getPrepDir(molFn)
          fnRoot = getBaseName(fnSmall)

          if self.method_conf.get() == 0:  # Genetic algorithm
              args = " '%s' --conformer --nconf %s --score rmsd --writeconformers -O '%s_conformers.mol2'" %\
                     (os.path.abspath(fnSmall), self.number_conf.get(), fnRoot)
          else:  # confab
              args = " '%s' --confab --original --verbose --conf %s --rcutoff %s -O '%s_conformers.mol2'" % \
                     (os.path.abspath(fnSmall), self.number_conf.get(), str(self.rmsd_cutoff.get()), fnRoot)
          try:
              runOpenBabel(protocol=self, args=args, cwd=self.getPrepDir())
          except:
              failedMols.append(fnRoot)

        if len(failedMols) > 0:
          with open(os.path.abspath(self._getExtraPath(f'failedConfomerGeneration_{it}.txt')), 'w') as f:
            for molFn in failedMols:
              f.write(molFn + '\n')

    def createOutput(self):
        """Create a set of Small Molecules as output with the path to:
              - Path to small molecule with H (mol2 format)
              - Path to conformers file (mol2 format)
        """
        self.mergeErrorFiles()

        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath())
        for fnRoot, mol in self.getOriginalBasenames().items():
            fnSmall = self.getPrepDir(f"{fnRoot}.mol2")
            if os.path.exists(fnSmall) and os.path.getsize(fnSmall) != 0:
                if not self.useLibrary.get():
                  mapFile = mol.writeMapFile(SmallMolecule(smallMolFilename=fnSmall), outDir=self._getExtraPath())

                if self.doConformers:
                    fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]
                    outDir = self._getExtraPath(fnRoot)
                    os.mkdir(outDir)
                    confFile = self.getPrepDir(f"{fnRoot}_conformers.mol2")
                    molFiles = splitConformerFile(confFile, outDir=outDir)
                    for molFile in molFiles:
                        confId = os.path.splitext(molFile)[0].split('-')[-1]

                        newSmallMol = SmallMolecule()
                        if not self.useLibrary.get():
                          newSmallMol.copy(mol, copyId=False)
                          newSmallMol.setMappingFile(pwobj.String(mapFile))

                        newSmallMol.setFileName(molFile)
                        newSmallMol.setConfId(confId)
                        newSmallMol._ConformersFile = pwobj.String(confFile)
                        outputSmallMolecules.append(newSmallMol)
                else:
                    newSmallMol = SmallMolecule(smallMolFilename=fnSmall, molName='guess')
                    if not self.useLibrary.get():
                      newSmallMol.setMappingFile(pwobj.String(mapFile))
                    outputSmallMolecules.append(newSmallMol)

        if outputSmallMolecules is not None:
            outputSmallMolecules.updateMolClass()
            self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

            inObj = self.inputLibrary if self.useLibrary.get() else self.inputSmallMolecules
            self._defineSourceRelation(inObj, outputSmallMolecules)


    # --------------------------- UTILS functions ------------------------------

    def mergeErrorFiles(self):
        lines = []
        with open(self._getExtraPath('failed.txt'), 'w') as f:
            for eFile in glob.glob(self._getExtraPath('failed_*.txt')):
                with open(eFile) as fIn:
                    for line in fIn:
                        lines.append(line)
                os.remove(eFile)
            lines = natural_sort(lines)
            for line in lines:
                f.write(line)

    def getLigandCode(self, paramsFile):
      with open(paramsFile) as f:
        code = f.readline().split()[1]
      return code
    
    def getPrepDir(self, path=''):
      return os.path.join(os.path.abspath(self._getExtraPath('prepared')), path)

    def reorderAtoms(self, inFile, outFile):
        '''Atom lines in the file are reordered so the atomNames numbers are in order (C2, C1, C3 -> C1, C2, C3)'''
        atomsDic = {}
        writeFirst, writeLast = '', ''
        if inFile.endswith('mol2'):
          with open(inFile) as fIn:
          # if inFile.endswith('pdb'):
          #     for line in fIn:
          #         sline = splitPDBLine(line)
          #         if sline[0] in ['ATOM', 'HETATM']:
          #             atomsDic[sline[2]] = line
          #         else:
          #             if atomsDic == {}:
          #                 writeFirst += line
          #             else:
          #                 writeLast += line

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

              with open(outFile, 'w') as f:
                f.write(writeFirst)
                idsDic = {}
                for i, key in enumerate(natural_sort(atomsDic.keys())):
                    idsDic[atomsDic[key].split()[0]] = str(i+1)
                    f.write(str(i+1).rjust(7) + atomsDic[key][7:])

                # Renaming connects
                f.write(writeLast.split('\n')[0] + '\n')
                for line in writeLast.split('\n')[1:]:
                    if line:
                        sline = line.split()
                        f.write(line[:6] + idsDic[sline[1]].rjust(5) + idsDic[sline[2]].rjust(5) + line[16:] + '\n')

        else:
          outFile = inFile

        return outFile

    def getOriginalBasenames(self):
      baseNames = {}
      if self.useLibrary.get():
        inLib = self.inputLibrary.get()
        with open(inLib.getFileName()) as f:
          for line in f:
            baseNames[line.split()[1]] = None
      else:
        for mol in self.inputSmallMolecules.get():
            fnSmall = mol.getFileName()
            baseNames[getBaseName(fnSmall)] = mol.clone()
      return baseNames


