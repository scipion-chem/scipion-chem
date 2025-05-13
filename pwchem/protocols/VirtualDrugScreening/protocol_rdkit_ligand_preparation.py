# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
import os, glob

from pyworkflow.protocol import params
from pyworkflow.protocol.params import LEVEL_ADVANCED
import pyworkflow.object as pwobj

from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import natural_sort, performBatchThreading, getBaseName
from pwchem.constants import RDKIT_DIC
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_filter import ProtocolBaseLibraryToSetOfMols

scriptName = 'ligand_preparation_script.py'

MFF_METHODS = ['MMFF94', 'MMFFp4s']
CONF_METHODS = ['KDG', 'ETDG', 'ETKDG', 'ETKDGv2', 'ETKDGv3', 'srETKDGv3']

class ProtChemRDKitPrepareLigands(ProtocolBaseLibraryToSetOfMols):
    """
    Prepare a set of molecules for use in a docking program (for example, Rosetta DARC).
    Sets the partial atomic charges, generates low-energy conformers. This is done with OpenBabel.
    """

    _label = 'RDKit Ligand preparation'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ Define the input parameters that will be used.
        """
        form.addSection(label='Input')
        form = self.addInputParams(form)

        group = form.addGroup("Molecule management")
        group.addParam('doHydrogens', params.BooleanParam, default=True,
                       label='ReAssign hydrogens: ',
                       help='Hydrogens are reomved (if present) and readded')
        group.addParam('doGasteiger', params.BooleanParam, default=True,
                       label='Recalculate gasteiger charges: ',
                       help='Charges of the molecule are recomputed using the gasteiger method')
        group.addParam('ffMethod', params.EnumParam,
                       choices=MFF_METHODS, default=0, label='Assign force-field to use: ',
                       help='Choose the force filed to perform the ligand optimization')
        group.addParam('splitFrag', params.BooleanParam, default=True,
                       label='Split molecule fragements: ', expertLevel=params.LEVEL_ADVANCED,
                       help='Split non-bonded fragments in the input molecule and prepare them separatedly')
        group.addParam('numAtoms', params.IntParam, expertLevel=params.LEVEL_ADVANCED,
                       default=5, condition="splitFrag", label='Min. number of atoms to keep:',
                       help='Set the minimum number of atoms of a fragment to be kept. This might be used to delete'
                            'water or other molecules stored in the molecule files.')

        conformers = form.addGroup("Conformers generation")
        conformers.addParam('doConformers', params.BooleanParam, default=False,
                            label='Do you want to generate conformers? ',
                            help='You can produce conformers of the ligand in order to do a better rigid docking')

        conformers.addParam('restrainMethod', params.EnumParam,
                            choices=CONF_METHODS,
                            default=4, condition = "doConformers",
                            label='Method of conformers generation',
                            help='Restrains method for conformers generation. '
                                 'https://greglandrum.github.io/rdkit-blog/conformers/exploration/2021/02/22/'
                                 'etkdg-and-distance-constraints.html')
        conformers.addParam('numConf', params.IntParam,
                            default=5, condition = "doConformers",
                            label='Max. number of conformers:',
                            help='Set the number of conformers generated by OpenBabel from the same molecule.')
        conformers.addParam('rmsd_cutoff', params.FloatParam, condition = "doConformers",
                            default=0.5,
                            label='RMSD cutoff:',
                            help='Minimum RMSD between the different output conformers',
                            expertLevel=LEVEL_ADVANCED)

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        # Insert processing steps
        aSteps = []
        nt = self.numberOfThreads.get()
        if nt <= 1: nt = 2

        iStep = self._insertFunctionStep(self.createInputStep, nt - 1, prerequisites=[])
        for it in range(nt-1):
            aSteps += [self._insertFunctionStep(self.preparationStep, it, prerequisites=[iStep])]
        self._insertFunctionStep(self.createOutput, prerequisites=aSteps)

    def preparationStep(self, it):
        """ Preparate the molecules and generate the conformers as specified
        """
        paramsPath = os.path.abspath(self._getExtraPath('inputParams_{}.txt'.format(it)))
        self.writeParamsFile(paramsPath, it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    def createOutput(self):
        """Create a set of Small Molecules as output
        """
        nt = self.numberOfThreads.get()
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='')

        getBasenames = self.getOriginalBasenames()

        outputMols = performBatchThreading(self.generateOutput, getBasenames, nt, cloneItem=False)
        for newSmallMol in outputMols:
          outputSmallMolecules.append(newSmallMol)

        self.mergeFails()

        if len(outputSmallMolecules) > 0:
            outputSmallMolecules.updateMolClass()
            self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
            inObj = self.inputLibrary if self.useLibrary.get() else self.inputSmallMolecules
            self._defineSourceRelation(inObj, outputSmallMolecules)


    # --------------------------- UTILS functions ------------------------------
    def mergeFails(self):
      failedIds = []
      for failFile in glob.glob(self._getExtraPath('failedPreparations_*.txt')):
        with open(failFile) as fr:
          failedIds += fr.read().strip().split('\n')
          os.remove(failFile)

      if len(failedIds) > 0:
        allFails = self.getFailsPath()
        with open(allFails, 'w') as f:
          failedIds.sort()
          f.write('\n'.join(failedIds))

    def getOriginalBasenames(self):
      baseNames = []
      if self.useLibrary.get():
        inLib = self.inputLibrary.get()
        with open(inLib.getFileName()) as f:
          for line in f:
            baseNames.append([line.split()[1], None])
      else:
        for mol in self.inputSmallMolecules.get():
            fnSmall = mol.getFileName()
            baseNames.append([getBaseName(fnSmall), mol.clone()])
      return baseNames

    def generateOutput(self, mols, molLists, it):
      outMols = []
      for basename, mol in mols:
        fPattern = "{}-*.sdf" if self.doConformers else "{}.sdf"
        tempSmall = self._getExtraPath(fPattern.format(basename))
        for molFile in sorted(list(glob.glob(tempSmall))):
          if not self.useLibrary.get():
            mapFile = mol.writeMapFile(SmallMolecule(smallMolFilename=molFile), outDir=self._getExtraPath(),
                                       mapBy='order')
          if os.path.exists(molFile) and os.path.getsize(molFile) != 0:
            confId = os.path.splitext(molFile)[0].split('-')[-1] if self.doConformers else 0

            newSmallMol = SmallMolecule()
            if not self.useLibrary.get():
              newSmallMol.copy(mol, copyId=False)
              newSmallMol.setMappingFile(pwobj.String(mapFile))

            newSmallMol.setFileName(molFile)
            newSmallMol.setConfId(confId)

            outMols.append(newSmallMol)
      molLists[it] = outMols

    def writeParamsFile(self, paramsFile, it):
        molFiles = self.getInputMolFiles(it)
        with open(paramsFile, 'w') as f:
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

            f.write('outputDir: {}\n'.format(os.path.abspath(self._getExtraPath())))
            f.write('doHydrogens: {}\n'.format(self.doHydrogens.get()))
            f.write('doGasteiger: {}\n'.format(self.doGasteiger.get()))
            f.write('ffMethod: {}\n'.format(self.getEnumText('ffMethod')))
            if self.doConformers.get():
                f.write('restrainMethod: {}\n'.format(self.getEnumText('restrainMethod')))
                f.write('numConf: {}\n'.format(self.numConf.get()))
                f.write('rmsThres: {}\n'.format(self.rmsd_cutoff.get()))
            if self.splitFrag.get():
                f.write('numAtoms: {}\n'.format(self.numAtoms.get()))
        return paramsFile

    def getLigandCode(self, paramsFile):
      with open(paramsFile) as f:
        code = f.readline().split()[1]
      return code


    def reorderAtoms(self, inFile, outFile):
        '''Atom lines in the file are reordered so the atomNames numbers are in order (C2, C1, C3 -> C1, C2, C3)'''
        atomsDic = {}
        writeFirst, writeLast = '', ''
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

          if inFile.endswith('mol2'):
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
              #Don't touch other kind of files
              return inFile

        return outFile

    def getFailsPath(self):
      return self._getPath('failedPreparations.txt')

    def _summary(self):
      if os.path.exists(self.getFailsPath()):
        return ['Preparation of some molecules failed. These are listed in {}'.
                  format(self.getFailsPath())]
