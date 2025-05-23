# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# # # *
# # # *
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

from os.path import abspath

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules, SmallMolecule
from pwchem.utils import *
from pwchem import Plugin as pwchemPlugin


scriptName = 'pharmacophore_filtering.py'

class ProtocolPharmacophoreFiltering(EMProtocol):
    """
    Perform the filtering of a set of small molecules that match an input pharmacophore.
    """
    _label = 'Pharmacophore filtering'
    stepsExecutionMode = params.STEPS_PARALLEL

    ##### -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', params.PointerParam, label="Input small molecules: ",
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      help='Select the set of small molecules to be filtered by matching the pharmacophore.')
        form.addParam('maxPerStep', params.IntParam, label="Maximum ligands processed per step: ",
                      expertLevel=params.LEVEL_ADVANCED, default=100,
                      help='Maximum number ligands processed per step')

        form.addParam('inputPharmacophore', params.PointerParam,
                      pointerClass='PharmacophoreChem', allowsNull=False,
                      label="Input pharmacophore: ",
                      help='Select the pharmacophore to use for the filtering.')

        group = form.addGroup('Matching parameters')
        group.addParam('downSample', params.BooleanParam, default=False,
                       label='Use matching downsamplig: ',
                       help="Whether to use downsampling while matching the pharmacophore.\n"
                            "Downsample will produce faster but less robust resuls ")

        group.addParam('nAlignments', params.IntParam, default=5,
                       label='Maximum number of alignments: ',
                       help="Maximum number of molecule-pharmacophore alignments")

        group.addParam('optimize', params.BooleanParam, default=False,
                       label='Optimize molecules to align: ',
                       help="Whether to optimize the molecules that match to pharmacophore."
                            "It may improve the alignment, but it will increase the computational cost")

        group.addParam('maxSSD', params.FloatParam, default=20,
                       label='Maximum deviation (SSD): ',
                       help="Maximum sum squares deviation of the ligand aligned to the pharmacophore to be considered")

        form.addParallelSection(threads=4, mpi=1)


        # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        aSteps = []
        nt = self.numberOfThreads.get()
        inLen = len(self.inputSmallMolecules.get())

        # If the number of input mols is so big, make more subsets than threads
        nSubsets = max(nt - 1, int(inLen / self.maxPerStep.get()))

        # Ensuring there are no more subsets than input molecules
        nSubsets = min(nSubsets, inLen)

        iStep = self._insertFunctionStep(self.createInputStep, nSubsets)
        for it in range(nSubsets):
            cStep = self._insertFunctionStep(self.convertStep, it, prerequisites=[iStep])
            aSteps += [self._insertFunctionStep(self.filterStep, it, prerequisites=cStep)]
        self._insertFunctionStep(self.createOutputStep, prerequisites=aSteps)

    def getInputDir(self):
      return self._getTmpPath()

    def getInputFile(self, it):
      return os.path.join(self.getInputDir(), f'inputLigandFiles_{it}.txt')

    def createInputStep(self, nSubsets):
        ligFiles = []
        for mol in self.inputSmallMolecules.get():
            ligFiles.append(os.path.abspath(mol.getPoseFile()) if mol.getPoseFile()
                            else os.path.abspath(mol.getFileName()))

        inputSubsets = makeSubsets(ligFiles, nSubsets, cloneItem=False)
        for it, fileSet in enumerate(inputSubsets):
            with open(self.getInputFile(it), 'w') as f:
                f.write(' '.join(fileSet))

    def convertStep(self, it):
        tmpDir, outDir = abspath(self._getTmpPath(f'convLigands_{it}')), self.getInputLigandsDir(it)
        os.makedirs(tmpDir), os.makedirs(outDir)

        inFile = self.getInputFile(it)
        with open(inFile) as f:
            molFiles = f.read().strip().split()

        for molFile in molFiles:
            if molFile.split('.')[-1] in ['pdb', 'mol2', 'sdf', 'mol']:
                os.link(molFile, os.path.join(outDir, getBaseFileName(molFile)))
            else:
                os.link(molFile, os.path.join(tmpDir, getBaseFileName(molFile)))

        if len(os.listdir(tmpDir)) > 0:
            # we need the input files in a RDKit readable format (not pdbqt for example)
            args = ' --multiFiles -iD "{}" --pattern "{}" -of pdb --outputDir "{}"'. \
                format(tmpDir, '*', outDir)
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

    def filterStep(self, it):
        paramsPath = self.writeParamsFile(it)

        args = ' {} {}'.format(paramsPath, abspath(self._getPath()))
        pwchemPlugin.runScript(self, scriptName, args, env=RDKIT_DIC, cwd=self._getPath())

    def createOutputStep(self):
        outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath())
        if self.inputPharmacophore.get().getProteinFile():
            outputSmallMolecules.setProteinFile(self.inputPharmacophore.get().getProteinFile())
            outputSmallMolecules.setDocked(True)

        resDic = self.parseResults()
        for mol in self.inputSmallMolecules.get():
            if mol.getPoseFile() and getBaseName(mol.getPoseFile()) in resDic:
                molBase = getBaseName(mol.getPoseFile())
            elif mol.getFileName() and getBaseName(mol.getFileName()) in resDic:
                molBase = getBaseName(mol.getFileName())
            else:
                continue

            newSmallMol = SmallMolecule()
            newSmallMol.copy(mol, copyId=False)
            newSmallMol.setGridId(1)
            newSmallMol.setPoseFile(os.path.relpath(resDic[molBase][0]))
            newSmallMol.setEnergy(resDic[molBase][1])

            poseId = resDic[molBase][0].split('.')[0].split('_')[-1]
            newSmallMol.setPoseId(poseId)
            outputSmallMolecules.append(newSmallMol)

        self._defineOutputs(outputSmallMolecules=outputSmallMolecules)

    # --------------------------- INFO functions -----------------------------------

    def _validate(self):
        vals = []
        return vals

    # --------------------------- UTILS functions -----------------------------------

    def writeParamsFile(self, it):
        paramsPath = abspath(self._getExtraPath(f'inputParams_{it}.txt'))
        convLigandNames = os.listdir(self.getInputLigandsDir(it))

        ligandFiles = []
        for fnLigand in convLigandNames:
            ligandFiles.append(os.path.join(self.getInputLigandsDir(it), fnLigand))

        pharmDic = self.inputPharmacophore.get().pharm2Dic()
        with open(paramsPath, 'w') as f:
            f.write('outputPath:: results.tsv\n')

            f.write('pharmDic:: {}\n'.format(str(pharmDic)))
            f.write('downSample:: {}\n'.format(self.downSample.get()))
            f.write('optimize:: {}\n'.format(self.optimize.get()))
            f.write('nAlignments:: {}\n'.format(self.nAlignments.get()))
            f.write('maxSSD:: {}\n'.format(self.maxSSD.get()))
            
            f.write('ligandFiles:: {}\n'.format(' '.join(ligandFiles)))

        return paramsPath

    def getInputLigandsDir(self, it):
        return abspath(self._getExtraPath(f'inputSmallMolecules_{it}'))

    def parseResults(self):
        resDic = {}
        for file in os.listdir(self._getPath()):
            if file.startswith('deviations'):
                with open(self._getPath(file)) as f:
                    f.readline()
                    for line in f:
                        oriBase, outFile, dev = line.split()
                        resDic[oriBase] = [outFile, dev]
        return resDic