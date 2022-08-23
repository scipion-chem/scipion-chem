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

import os, re
from os.path import abspath

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
from pwchem import Plugin as pwchemPlugin
from pwchem.scripts.PBVS_script import *

FEATURE_LABELS_SIMPLE = ["Donor", "Acceptor", "Hydrophobe", "Aromatic"]
FEATURE_LABELS_ADVANCED = ["LumpedHydrophobe", "PosIonizable", "NegIonizable", "ZnBinder"]

scriptName = 'PBVS_script.py'

class ProtocolPharmacophoreFiltering(EMProtocol):
    """
    Perform the construction of a consensus pharmacophore from a set of PDB ligands and
    the PBVS (using the pharmacophore) against a set of molecules in smi, pdb, mol, mol2 or sdf file format.

    """
    #Nombre del protocolo (aparece en grande arriba)
    _label = 'Pharmacophore filtering'

    ##### -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form.addParam('inputLigands', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input reference ligands: ",
                      help='Select the ligands PDB files')

        group = form.addGroup('Features')
        for feat in FEATURE_LABELS_SIMPLE:
            group.addParam(feat, params.BooleanParam, default=True,
                           label='Use {} as feature: '.format(feat),
                           help="Use {} as feature for building the pharmacophore.\nBase features defined in "
                                "https://github.com/rdkit/rdkit/blob/master/Data/BaseFeatures.fdef".format(feat))

        for feat in FEATURE_LABELS_ADVANCED:
            group.addParam(feat, params.BooleanParam, default=False,
                           label='Use {} as feature: '.format(feat), expertLevel=params.LEVEL_ADVANCED,
                           help="Use {} as feature for building the pharmacophore.\nBase features defined in "
                                "https://github.com/rdkit/rdkit/blob/master/Data/BaseFeatures.fdef"
                                "".format(feat))

        form.addSection(label='Clustering')
        group = form.addGroup('Clustering')
        group.addParam('method', params.EnumParam, default=0,
                       label='Clustering method: ', choices=[DBSCAN, KMEANS],
                       help="Which method to use for clustering the molecules features for building the "
                            "pharmacophore")

        group.addParam('eps', params.FloatParam, default=2.0,
                       label='Maximum neighbor distance: ', condition='method==0',
                       help="The maximum distance between two samples for one to be considered as in the neighborhood "
                            "of the other. This is not a maximum bound on the distances of points within a cluster.\n"
                            "https://scikit-learn.org/stable/modules/generated/sklearn.cluster.DBSCAN.html")

        group.addParam('kNumber', params.IntParam, default=3,
                       label='Number of clusters for each feature: ', condition='method==1',
                       help="Number of clusters for each feature using the Kmeans clustering")

        group.addParam('minType', params.EnumParam, default=1,
                       label='Min size meaning: ', choices=[ABSOLUTE, PROP_MOLS, PROP_FEATS],
                       help="How to define the minimum size of a feature cluster. \n"
                            "Absolute: minimum number of features in absolute value\n"
                            "Molecules proportion: minimum proportion of molecules adding a feature to the cluster "
                            "(e.g: 0.5 would mean that at least half of the molecules must have this kind of feature "
                            "in the cluster\n"
                            "Features proportion: minimum proportion of features adding to the cluster "
                            "(compared to the total number of element of each feature)")

        group.addParam('minSize', params.FloatParam, default=0.3,
                       label='Minimum size of a cluster to be considered: ',
                       help="Minimum size of the cluster to be considered (check type of number in 'Min size meaning'"
                            " parameter")

        group.addParam('topClusters', params.IntParam, default=-1,
                       label='Maximum number of clusters to pick for each feature: ',
                       help="Use only the x biggest clusters for each feature for building the pharmacophore.\n"
                            "If -1, all clusters will be used")

        # --------------------------- STEPS functions ------------------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')
        self._insertFunctionStep('filterStep')
        self._insertFunctionStep('createOutputStep')

    def convertStep(self):
        tmpDir, outDir = abspath(self._getTmpPath('convLigands')), self.getInputLigandsDir()
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        for ligand in self.inputLigands.get():
            inFile = ligand.getPoseFile()
            if not inFile.split('.')[-1] in ['pdb', 'mol2', 'sdf', 'mol']:
                shutil.copy(inFile, os.path.join(tmpDir, os.path.basename(inFile)))
            else:
                shutil.copy(inFile, os.path.join(outDir, os.path.basename(inFile)))

        if len(os.listdir(tmpDir)) > 0:
            # we need the input files in a RDKit readable format (not pdbqt for example)
            args = ' --multiFiles -iD "{}" --pattern "{}" -of pdb --outputDir "{}"'. \
                format(tmpDir, '*', outDir)
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

        smiDir = self.getInputSMIDir()
        if not os.path.exists(smiDir):
            os.makedirs(smiDir)
        args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
            format(outDir, '*', smiDir)
        pwchemPlugin.runScript(self, 'obabel_IO.py', args, env='plip', cwd=outDir)

    def filterStep(self):
        paramsPath = self.writeParamsFile()

        args = ' {} {}'.format(paramsPath, abspath(self._getPath()))
        pwchemPlugin.runScript(self, scriptName, args, env='rdkit', cwd=self._getPath())

    # --------------- INFO functions -------------------------

    def _citations(self):
        return ["@article{wójcikowski_zielenkiewicz_siedlecki_2015, title={Open drug discovery toolkit (ODDT): A new open-source player in the Drug Discovery Field}, volume={7}, DOI={10.1186/s13321-015-0078-2}, number={1}, journal={Journal of Cheminformatics}, author={Wójcikowski, Maciej and Zielenkiewicz, Piotr and Siedlecki, Pawel}, year={2015}}"]

    def _methods(self):
        methods ="The first part of the process involves the generation of the consensus pharmacophore from the PDB ligands. The second protocol step is the analysis of starting molecules pharmacophore characteristics using the consensus pharmacophore."
        return methods

    # --------------------------- UTILS functions -----------------------------------

    def createOutputStep(self):
        if False:
            newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

            mols = self.inputSmallMolecules.get()
            for mol in mols:
                file = abspath(mol.getFileName())
                if file in filtered_molecules:
                    newMols.append(mol)

            newMols.updateMolClass()
            self._defineOutputs(outputSmallMolecules=newMols)

    def writeParamsFile(self):
        paramsPath = abspath(self._getExtraPath('inputParams.txt'))
        ligandsFiles = os.listdir(self.getInputLigandsDir())

        outFiles, smiles_ligand = [], []
        smiFileDic = self.getBaseNameDic(self.getInputSMIDir())

        smiDic = {}
        for fnLigand in ligandsFiles:
            ligandBase = getBaseFileName(fnLigand)
            with open(smiFileDic[ligandBase]) as f:
                smile = f.read().split()[0].strip()

            outFiles.append(os.path.join(self.getInputLigandsDir(), fnLigand))
            smiDic[ligandBase] = smile

        with open(paramsPath, 'w') as f:
            # Esta linea vale --> es la que lleva al archivo output

            f.write('outputPath:: results.tsv\n')

            f.write('method:: {}\n'.format(self.getEnumText('method')))

            f.write('kNumber:: {}\n'.format(self.kNumber.get()))
            f.write('eps:: {}\n'.format(self.eps.get()))

            f.write('minType:: {}\n'.format(self.getEnumText('minType')))
            f.write('minSize:: {}\n'.format(self.minSize.get()))
            f.write('topClusters:: {}\n'.format(self.topClusters.get()))

            featStr = 'features:: '
            for feat in FEATURE_LABELS_SIMPLE + FEATURE_LABELS_ADVANCED:
                if getattr(self, feat):
                    featStr += '{} '.format(feat)
            f.write(featStr + "\n")

            f.write('ligandSmiles:: ' + str(smiDic) + "\n")
            f.write('ligandFiles:: {}\n'.format(' '.join(outFiles)))

        return paramsPath

    def getInputLigandsDir(self):
        return abspath(self._getExtraPath('inputLigands'))

    def getInputSMIDir(self):
        return abspath(self._getExtraPath('inputLigandsSMI'))

    def getBaseNameDic(self, inDir):
        bDic = {}
        for file in os.listdir(inDir):
            bDic[getBaseFileName(file)] = abspath(os.path.join(inDir, file))
        return bDic