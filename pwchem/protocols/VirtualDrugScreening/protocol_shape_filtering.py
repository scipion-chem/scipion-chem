# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
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

# General imports
import os

# Scipion em imports
from pyworkflow.protocol import params
import pyworkflow.object as pwobj

# Plugin imports
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import RDKIT_DIC
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_filter import ProtocolBaseLibraryToSetOfMols
from pwchem.utils import natural_sort, getBaseName

scriptName = 'shape_distances_script.py'
shapePrograms = ['RDKit', 'Shape-it']


class ProtocolShapeDistances(ProtocolBaseLibraryToSetOfMols):
    """
    Performs shape distance calculation of a set of molecules against another reference molecule
    """

    _label = 'Shape Distance calculation'
    _distanceType = ['Tanimoto Distance', 'Protrude Distance', 'RMSD']
    _shapeitScore = ['TANIMOTO', 'TVERSKY_REF', 'TVERSKY_DB']
    stepsExecutionMode = params.STEPS_PARALLEL

    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form = self.addInputParams(form)

        group = form.addGroup('Reference molecule')
        group.addParam('inputRefSmallMolecules', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Reference Small Molecules: ",
                       help='Select the molecules where the reference molecule is stored')

        group.addParam('inputReferenceMolecule', params.StringParam,
                       label="Reference molecule: ",
                       help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('program', params.EnumParam,
                       choices=shapePrograms, default=1,
                       label='Program to calculate distance: ')
        group.addParam('distanceType', params.EnumParam, default=0, label='Distance type: ',
                       choices=self._distanceType, condition='program==0',
                       help="Chosen distance type to perform the calculation")
        group.addParam('distanceTypeShapeit', params.EnumParam, default=0, label='Shape-it score: ',
                       choices=self._shapeitScore, condition='program==1',
                       help="Shape-it score to take into account"
                            "\nTanimoto: V(overlap) / (V(ref) + V(db) - V(overlap))"
                            "\nTversky_ref: V(overlap) / V(ref)"
                            "\nTversky_db: V(overlap) / V(db)"
                            "\n(https://www.sciencedirect.com/science/article/pii/S109332630800048X?via%3Dihub)")

        group.addParam('prealign', params.BooleanParam, default=True,
                       label='Prealign molecules: ', condition='program==0 and distanceType in [0, 1]',
                       help='Tries to prealign the molecules before the distance calculation by using substructure '
                            'matching. If no substructure matching is found, molecules will be left as they are')
        group.addParam('prealignOrder', params.BooleanParam, default=False,
                       label='Try every atom reordering: ', condition='program==0 and prealign',
                       help='Tries every permutation in the atom order for the molecule alignment. As specified by '
                            'RDKit, for some molecules it will lead to combinatorial explosion, especially if hydrogens'
                            ' are present.')

        form.addParam('hydrogen', params.BooleanParam, default=True, condition='program==0',
                      label='Do you want to ignore the hydrogens in the calculation?')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def mainStep(self, it):
        if self.program.get() == 0:
            self.runRDKitDistance(it)
        else:
            self.runShapeItDistance(it)

    def createOutputStep(self):
        scoreDic = self.getOutputScoreDic()
        if self.useLibrary.get():
          self.defineLibraryOutput(scoreDic)
        else:
          self.defineMolOutput(scoreDic)

    def defineMolOutput(self, scoreDic):
      mols = self.inputSmallMolecules.get()
      newMols = SetOfSmallMolecules.createCopy(mols, self._getPath(), copyInfo=True)
      nSubsets = self.getNSubsets()

      for i, mol in enumerate(mols):
        score = scoreDic[i % nSubsets].pop(0)
        mol._shapeDistance = pwobj.Float(score)
        newMols.append(mol)

      newMols.updateMolClass()
      self._defineOutputs(outputSmallMolecules=newMols)

    def defineLibraryOutput(self, scoreDic):
      inLib = self.inputLibrary.get()
      mapDic = inLib.getLibraryMap(inverted=True, fullLine=True)
      oLibFile = self._getPath('outputLibrary.smi')
      nSubsets = self.getNSubsets()

      with open(oLibFile, 'w') as f:
          for i, (smi, line) in enumerate(mapDic.items()):
              score = scoreDic[i % nSubsets].pop(0)
              f.write(f'{line}\t{score}\n')

      prevHeaders = inLib.getHeaders()
      outputLib = inLib.clone()
      outputLib.setFileName(oLibFile)
      outputLib.setHeaders(prevHeaders + ['Shape_distance'])
      self._defineOutputs(outputLibrary=outputLib)



    # --------------------------- UTILS functions -----------------------------------
    def getOutputScoreDic(self):
        '''Return a dictionary as {it: [scores]} with the scores of each thread
        '''
        scoresDic = {}
        for file in os.listdir(self._getPath()):
          if file.startswith('results_'):
            it = int(getBaseName(file).split('_')[-1])
            if self.program.get() == 1:
              scoresDic[it] = self.parseResultsShapeit(self._getPath(file))
            else:
              scoresDic[it] = self.parseResults(self._getPath(file))
        return scoresDic


    def runRDKitDistance(self, it):
        paramsPath = self.writeParamsFile(it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    def runShapeItDistance(self, i):
        allMolsFile = self.buildMolsFile(i)

        paramsPath = f' -r {self.getRefFile()} -d {allMolsFile} -s results_{i}.tsv --noRef'
        try:
            Plugin.runShapeIt(self, paramsPath, cwd=self._getPath())
        except:
            print('Shape-it could not be executed because it was not installed properly. \nIf you want to use it, '
                  'you will need to install it manually (https://github.com/rdkit/shape-it) and define its home dir '
                  'as "SHAPEIT_HOME=<path_to_shape-it>" in the scipion.conf file')

    def writeParamsFile(self, it):
        molFiles = []
        paramsFile = os.path.abspath(self._getExtraPath(f'inputParams_{it}.txt'))
        f = open(paramsFile, 'w')
        f.write(f'outputPath: results_{it}.tsv\n')
        itMolFiles = self.getInputMolFiles(it)
        for molFile in itMolFiles:
            molFiles.append(os.path.abspath(molFile))

        f.write('distanceType: {}\n'.format(self.getDistance()))
        f.write('prealign: {}\n'.format(self.prealign.get()))
        f.write('prealignOrder: {}\n'.format(self.prealignOrder.get()))
        f.write('ignoreHydrogen: {}\n'.format(self.hydrogen.get()))

        f.write('referenceFile: {}\n'.format(self.getRefFile()))
        f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile
    # --------------- INFO functions -------------------------

    def _citations(self):
        return [
            "@misc{landrum _2021, title={Rdkit.chem.rdmolalign module¶}, url={https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html}, journal={rdkit.Chem.rdMolAlign module - The RDKit 2022.03.1 documentation}, author={Landrum , Greg}, year={2021}}"]

    def _methods(self):
        methods = "This calculation is based on the analysis of the shape of the query molecule/ligand used as a reference in the search for similar structures. A commonly used measure for comparing the structure of two molecules in virtual screening programs is the calculation of the root mean square deviation (RMSD) between the atoms of two molecules. The RMSD is a distance which describes the structural difference between two topologies. The lower the RMSD between two structures, the greater the similarity between them."
        methods2 = "RDKit current tools for calculating distances between molecules of unequal size are Protrude Distance and Tanimoto distance. Protrude Distance focusses on the volume mismatch and is defined as the percentage of the larger molecule which protrudes/exceeds from the smaller molecule. Tanimoto Distance measures similarity between finite sample sets and is defined as the size of the intersection divided by the size of the union of the sample sets."
        methods3 = "RDKit does not have a tool that performs an alignment prior to the calculation of similarity, so the results provided by these tools must be carefully analyzed by users."
        return methods, methods2, methods3

    # --------------------------- UTILS functions -----------------------------------
    def buildMolsFile(self, i):
      molFiles = self.getInputMolFiles(i)
      ext = os.path.splitext(molFiles[0])[1]
      allMolsFile = os.path.abspath(self._getTmpPath(f'allmols_{i}{ext}'))

      with open(allMolsFile, 'w') as f:
        if ext == '.smi':
          for molFile in molFiles:
            with open(molFile) as fIn:
              f.write(fIn.read().split()[0] + '\n')
        else:
          for molFile in molFiles:
            with open(molFile) as fIn:
              f.write(fIn.read())
      return allMolsFile

    def getRefFile(self):
        if not hasattr(self, 'refFile'):
            for mol in self.inputRefSmallMolecules.get():
                if mol.__str__() == self.inputReferenceMolecule.get():
                    self.refFile = os.path.abspath(mol.getFileName())
        return self.refFile

    def getDistance(self):
        return self.getEnumText('distanceType')

    def parseResults(self, outputFile):
        scores = []
        with open(outputFile) as read_tsv:
            for row in read_tsv:
                if row[0] != "#":
                    row1 = row.split("\t")
                    scores.append(row1[1])

        return scores

    def parseResultsShapeit(self, resultsFile):
        scores = []
        with open(resultsFile) as f:
            f.readline()
            for line in f:
                simil = float(line.split()[2 + self.distanceTypeShapeit.get()])
                scores.append(1-simil)
        return scores
