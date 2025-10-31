# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, pickle
from os.path import abspath

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import PharmFeature, PharmacophoreChem
from pwchem.constants import *
from pwchem.utils import *
from pwchem import Plugin as pwchemPlugin

ABSOLUTE, PROP_MOLS, PROP_FEATS = ['Absolute', 'Molecules proportion', 'Features proportion']
DBSCAN, KMEANS = ['DBSCAN', 'KMeans']

scriptName = 'pharmacophore_generation.py'

class ProtocolPharmacophoreGeneration(EMProtocol):
    """
    Perform the construction of a consensus pharmacophore from a set of docked small molecules
    """
    _label = 'Pharmacophore generation'

    ##### -------------------------- DEFINE param functions ----------------------

    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input reference ligands: ",
                      help='Select the ligands PDB files')

        form.addParam('inputReceptor', params.PointerParam,
                      pointerClass='AtomStruct', allowsNull=True,
                      label="Input receptor: ",  expertLevel=params.LEVEL_ADVANCED,
                      help='Input receptor in case the ligands are not associated with one')

        form.addParam('propRadii', params.FloatParam, default=0.5,
                       label='Pharmacophore radii proportion: ', expertLevel=params.LEVEL_ADVANCED,
                       help="The radii of each feature (sphere) in the pharmacophore is determined by the distance "
                            "from the centriod to one of the points in the cluster. This number determines the "
                            "proportion of points in the cluster that must be included in the resulting sphere."
                            "(e.g: if 0.5, half of the feature points in the cluster will be placed inside the sphere)"
                            "\n Minimum radius is set to 1 (A)")

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
        self._insertFunctionStep('generationStep')
        self._insertFunctionStep('createOutputStep')

    def convertStep(self):
        tmpDir, outDir = abspath(self._getTmpPath('convLigands')), self.getInputLigandsDir()
        if not os.path.exists(tmpDir):
            os.makedirs(tmpDir)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        for ligand in self.inputSmallMolecules.get():
            inFile = ligand.getPoseFile() if ligand.getPoseFile() else ligand.getFileName()
            if not inFile.split('.')[-1] in ['pdb', 'mol2', 'sdf', 'mol']:
                shutil.copy(inFile, os.path.join(tmpDir, os.path.basename(inFile)))
            else:
                shutil.copy(inFile, os.path.join(outDir, os.path.basename(inFile)))

        if len(os.listdir(tmpDir)) > 0:
            # we need the input files in a RDKit readable format (not pdbqt for example)
            args = ' --multiFiles -iD "{}" --pattern "{}" -of pdb --outputDir "{}"'. \
                format(tmpDir, '*', outDir)
            pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

        smiDir = self.getInputSMIDir()
        if not os.path.exists(smiDir):
            os.makedirs(smiDir)
        args = ' --multiFiles -iD "{}" --pattern "{}" -of smi --outputDir "{}"'. \
            format(outDir, '*', smiDir)
        pwchemPlugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)

    def generationStep(self):
        paramsPath = self.writeParamsFile()

        args = ' {} {}'.format(paramsPath, abspath(self._getPath()))
        pwchemPlugin.runScript(self, scriptName, args, env=RDKIT_DIC, cwd=self._getPath())

    def createOutputStep(self):
        cenPath = os.path.abspath(self._getPath('cluster_centers.pkl'))
        with open(cenPath, "rb") as clCenters:
          centers = pickle.load(clCenters)

        radPath = os.path.abspath(self._getPath('cluster_radii.pkl'))
        with open(radPath, "rb") as clRadii:
          radii = pickle.load(clRadii)

        outPharm = PharmacophoreChem().create(outputPath=self._getPath())
        if self.inputSmallMolecules.get().getProteinFile():
            outPharm.setProteinFile(self.inputSmallMolecules.get().getProteinFile())

        elif self.inputReceptor.get().getFileName():
            outPharm.setProteinFile(self.inputReceptor.get().getFileName())

        for feat in radii:
            feat_radii = radii[feat]
            for i, loc in enumerate(centers[feat]):
                pharmFeat = PharmFeature(type=feat, radius=feat_radii[i],
                                         x=loc[0], y=loc[1], z=loc[2])
                outPharm.append(pharmFeat)

        self._defineOutputs(outputPharmacophore=outPharm)

    # --------------------------- UTILS functions -----------------------------------

    def writeParamsFile(self):
        paramsPath = abspath(self._getExtraPath('inputParams.txt'))
        ligandsFiles = os.listdir(self.getInputLigandsDir())

        outFiles, smiles_ligand = [], []
        smiFileDic = self.getBaseNameDic(self.getInputSMIDir())

        smiDic = {}
        for fnLigand in ligandsFiles:
            ligandBase = getBaseName(fnLigand)
            with open(smiFileDic[ligandBase]) as f:
                smile = f.read().split()[0].strip()

            outFiles.append(os.path.join(self.getInputLigandsDir(), fnLigand))
            smiDic[ligandBase] = smile

        with open(paramsPath, 'w') as f:
            # Esta linea vale --> es la que lleva al archivo output

            f.write('method:: {}\n'.format(self.getEnumText('method')))

            f.write('kNumber:: {}\n'.format(self.kNumber.get()))
            f.write('eps:: {}\n'.format(self.eps.get()))

            f.write('minType:: {}\n'.format(self.getEnumText('minType')))
            f.write('minSize:: {}\n'.format(self.minSize.get()))
            f.write('topClusters:: {}\n'.format(self.topClusters.get()))

            f.write('propRadii:: {}\n'.format(self.propRadii.get()))

            featStr = 'features:: '
            for feat in FEATURE_LABELS_SIMPLE + FEATURE_LABELS_ADVANCED:
                if getattr(self, feat):
                    featStr += '{} '.format(feat)
            f.write(featStr + "\n")

            f.write('ligandSmiles:: ' + str(smiDic) + "\n")
            f.write('ligandFiles:: {}\n'.format(' '.join(outFiles)))

        return paramsPath

    def getInputLigandsDir(self):
        return abspath(self._getExtraPath('inputSmallMolecules'))

    def getInputSMIDir(self):
        return abspath(self._getExtraPath('inputSmallMoleculesSMI'))

    def getBaseNameDic(self, inDir):
        bDic = {}
        for file in os.listdir(inDir):
            bDic[getBaseName(file)] = abspath(os.path.join(inDir, file))
        return bDic