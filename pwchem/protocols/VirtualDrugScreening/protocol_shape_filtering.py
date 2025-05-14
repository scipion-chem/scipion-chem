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
from pwem.protocols import EMProtocol

# Plugin imports
from pwchem import Plugin
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import RDKIT_DIC

scriptName = 'shape_distances_script.py'
shapePrograms = ['RDKit'] #, 'Shape-it']


class ProtocolShapeDistancesFiltering(EMProtocol):
    """
    Performs shape filtering of a set of ligands by the calculation of Tanimoto and Protrude distances between
    a molecule in smi format and a query.
    """

    _label = 'Shape Distance filtering'
    _distanceType = ['Tanimoto Distance', 'Protrude Distance', 'RMSD']
    _shapeitScore = ['TANIMOTO', 'TVERSKY_REF', 'TVERSKY_DB']

    def _defineParams(self, form):
        """ """
        form.addSection(label='Params')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules: ",
                      help='Select the molecules to be filtered')

        form.addParam('inputRefSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Reference Small Molecules: ",
                      help='Select the molecules where the reference molecule is stored')

        form.addParam('inputReferenceMolecule', params.StringParam,
                      label="Reference molecule: ",
                      help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('program', params.EnumParam,
                       choices=shapePrograms, default=0,
                       label='Program to calculate distance: ')
        group.addParam('distanceType', params.EnumParam, default=0, label='Distance type: ',
                       choices=self._distanceType, condition='program==0',
                       help="Chosen distance type to perform the filtering")
        group.addParam('distanceTypeShapeit', params.EnumParam, default=0, label='Shape-it score: ',
                       choices=self._shapeitScore, condition='program==1',
                       help="Shape-it score to take into account"
                            "\nTanimoto: V(overlap) / (V(ref) + V(db) - V(overlap))"
                            "\nTversky_ref: V(overlap) / V(ref)"
                            "\nTversky_db: V(overlap) / V(db)"
                            "\n(https://www.sciencedirect.com/science/article/pii/S109332630800048X?via%3Dihub)")

        group.addParam('prealign', params.BooleanParam, default=True,
                       label='Prealign molecules: ', condition='distanceType in [0, 1]',
                       help='Tries to prealign the molecules before the distance calculation by using substructure '
                            'matching. If no substructure matching is found, molecules will be left as they are')
        group.addParam('prealignOrder', params.BooleanParam, default=False,
                       label='Try every atom reordering: ', condition='prealign',
                       help='Tries every permutation in the atom order for the molecule alignment. As specified by '
                            'RDKit, for some molecules it will lead to combinatorial explosion, especially if hydrogens'
                            ' are present.')

        form.addParam('hydrogen', params.BooleanParam, default=True, condition='program==0',
                      label='Do you want to ignore the hydrogens in the filtering?')

        group.addParam('cut', params.FloatParam, default=0.5, label='Filter cut-off: ',
                       help="Filter cut-off for distance or shape-it score")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def filterStep(self):
        if self.program.get() == 0:
            self.describeFilter()
        else:
            self.shapeitFilter()

    def createOutputStep(self):
        if self.program.get() == 0:
            resultsFile = self._getPath("results.tsv")
            filtered_molecules_dict = self.parseResults(resultsFile)
        else:
            resultsFile = self._getPath("shape-it.tsv")
            filtered_molecules_dict = self.parseResultsShapeit(resultsFile)

        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)
        filtered_molecules = list(filtered_molecules_dict.keys())

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            molFile = os.path.abspath(mol.getFileName())
            if molFile in filtered_molecules:
                mol._shapeDistance = pwobj.String(filtered_molecules_dict[molFile])
                newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)


    # --------------------------- UTILS functions -----------------------------------

    def describeFilter(self):  # , receptorFile):
        mols = self.inputSmallMolecules.get()
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, mols)  # , receptorFile)

        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    def shapeitFilter(self):
        mols = self.inputSmallMolecules.get()
        ext = os.path.splitext(mols.getFirstItem().getFileName())[1]
        allMolsFile = os.path.abspath(self._getTmpPath('allmols' + ext))
        with open(allMolsFile, 'w') as f:
            for mol in mols:
                with open(mol.getFileName()) as fIn:
                    f.write(fIn.read())

        paramsPath = ' -r {} -d {} -s shape-it.tsv'.\
            format(self.getRefFile(), allMolsFile, self.getEnumText('distanceTypeShapeit'), self.cut.get())
        try:
            Plugin.runShapeIt(self, './shape-it', paramsPath, cwd=self._getPath())
        except:
            print('Shape-it could not be executed because it was not installed properly. \nIf you want to use it, '
                  'you will need to install it manually (https://github.com/rdkit/shape-it) and define its home dir '
                  'as "SHAPEIT_HOME=<path_to_shape-it>" in the scipion.conf file')

    def writeParamsFile(self, paramsFile, molsScipion):
        molFiles = []
        f = open(paramsFile, 'w')
        f.write('outputPath: results.tsv\n')
        for mol in molsScipion:
            molFiles.append(os.path.abspath(mol.getFileName()))

        f.write('referenceFile: {}\n'.format(self.getRefFile()))
        f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        f.write('distanceType:{}\n'.format(self.getDistance()))
        f.write('prealign:{}\n'.format(self.prealign.get()))
        f.write('prealignOrder:{}\n'.format(self.prealignOrder.get()))
        f.write('ignoreHydrogen:{}\n'.format(self.hydrogen.get()))
        f.write('cut-off: {}\n'.format(self.cut.get()))

        return paramsFile
    # --------------- INFO functions -------------------------

    def _citations(self):
        return [
            "@misc{landrum _2021, title={Rdkit.chem.rdmolalign module¶}, url={https://www.rdkit.org/docs/source/rdkit.Chem.rdMolAlign.html}, journal={rdkit.Chem.rdMolAlign module - The RDKit 2022.03.1 documentation}, author={Landrum , Greg}, year={2021}}"]

    def _methods(self):
        methods = "This filter is based on the analysis of the shape of the query molecule/ligand used as a reference in the search for similar structures. A commonly used measure for comparing the structure of two molecules in virtual screening programs is the calculation of the root mean square deviation (RMSD) between the atoms of two molecules. The RMSD is a distance which describes the structural difference between two topologies. The lower the RMSD between two structures, the greater the similarity between them."
        methods2 = "RDKit current tools for calculating distances between molecules of unequal size are Protrude Distance and Tanimoto distance. Protrude Distance focusses on the volume mismatch and is defined as the percentage of the larger molecule which protrudes/exceeds from the smaller molecule. Tanimoto Distance measures similarity between finite sample sets and is defined as the size of the intersection divided by the size of the union of the sample sets."
        methods3 = "RDKit does not have a tool that performs an alignment prior to the calculation of similarity, so the results provided by these tools must be carefully analyzed by users."
        return methods, methods2, methods3

    # --------------------------- UTILS functions -----------------------------------
    def getRefFile(self):
        if not hasattr(self, 'refFile'):
            for mol in self.inputRefSmallMolecules.get():
                if mol.__str__() == self.inputReferenceMolecule.get():
                    self.refFile = os.path.abspath(mol.getFileName())
        return self.refFile

    def getDistance(self):
        return self.getEnumText('distanceType')

    def parseResults(self, outputFile):
        molecules = {}
        with open(outputFile) as read_tsv:
            for row in read_tsv:
                if row[0] != "#":
                    row1 = row.split("\t")
                    molecules[row1[0]] = row1[1]

        return molecules

    def parseResultsShapeit(self, resultsFile):
        molecules_dic = {}
        with open(resultsFile) as f:
            lines = f.readlines()
            for i, mol in enumerate(self.inputSmallMolecules.get()):
                score = lines[i+1].split()[2 + self.distanceTypeShapeit.get()]
                if float(score) > self.cut.get():
                    molecules_dic[mol.getFileName()] = score
        return molecules_dic
