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
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input small molecules: ",
                      help='Select the set of small molecules to be filtered by matching the pharmacophore.')
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
        cStep = self._insertFunctionStep('convertStep')

        aSteps = []
        nt = self.numberOfThreads.get()
        for it in range(nt-1):
            aSteps += [self._insertFunctionStep('filterStep', it, prerequisites=[cStep])]
        self._insertFunctionStep('createOutputStep', prerequisites=aSteps)

    def convertStep(self):
        tmpDir, outDir = abspath(self._getTmpPath('convLigands')), self.getInputLigandsDir()
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        for ligand in self.inputSmallMolecules.get():
            if ligand.getPoseFile():
                inFile = ligand.getPoseFile()
            else:
                inFile = ligand.getFileName()

            if not inFile.split('.')[-1] in ['pdb', 'mol2', 'sdf', 'mol']:
                shutil.copy(inFile, os.path.join(tmpDir, os.path.basename(inFile)))
            else:
                shutil.copy(inFile, os.path.join(outDir, os.path.basename(inFile)))

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
        docked = False
        if self.inputPharmacophore.get().getProteinFile():
            outputSmallMolecules.setProteinFile(self.inputPharmacophore.get().getProteinFile())
            docked = True
            outputSmallMolecules.setDocked(docked)

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
        convLigandNames = os.listdir(self.getInputLigandsDir())
        cLigNames = makeSubsets(convLigandNames, nt=self.numberOfThreads.get()-1, cloneItem=False)[it]

        ligandFiles = []
        for fnLigand in cLigNames:
            ligandFiles.append(os.path.join(self.getInputLigandsDir(), fnLigand))

        pharmDic = self.inputPharmacophore.get().pharm2Dic()

        with open(paramsPath, 'w') as f:
            # Esta linea vale --> es la que lleva al archivo output

            f.write('outputPath:: results.tsv\n')

            f.write('pharmDic:: {}\n'.format(str(pharmDic)))
            f.write('downSample:: {}\n'.format(self.downSample.get()))
            f.write('optimize:: {}\n'.format(self.optimize.get()))
            f.write('nAlignments:: {}\n'.format(self.nAlignments.get()))
            f.write('maxSSD:: {}\n'.format(self.maxSSD.get()))
            
            f.write('ligandFiles:: {}\n'.format(' '.join(ligandFiles)))

        return paramsPath

    def getInputLigandsDir(self):
        return abspath(self._getExtraPath('inputSmallMolecules'))

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