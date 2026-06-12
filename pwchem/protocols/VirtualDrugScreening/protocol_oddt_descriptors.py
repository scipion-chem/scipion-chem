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
This protocol is used to generate descriptors for ligands docked to a protein with the Open
Drug Discovery Toolkit (ODDT, https://github.com/oddt/oddt)

"""
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pyworkflow.utils import Message
from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import *
import os

from pwchem import Plugin
from pwchem.utils import fillEmptyAttributes

scriptName = 'descriptors_docking_oddt.py'

class ProtocolODDTDescriptors(EMProtocol):
    """
    AI Generated:

This protocol extracts molecular descriptors for ligands that have been docked to a protein structure using the Open Drug Discovery Toolkit (ODDT).

It is designed to evaluate ligand–protein docking results by computing scoring functions that estimate binding quality based on different empirical or machine-learning-based models.

Core Concepts
-------------
Docked Small Molecules:
    A set of ligands that have been positioned inside a protein binding site using a docking procedure. Each molecule includes its pose information relative to the receptor.

Receptor Protein:
    The protein structure used as the target for docking. It is required for computing structure-dependent descriptors.

Molecular Descriptors:
    Numerical scores that estimate binding affinity or interaction quality between ligand and receptor. These descriptors are computed using ODDT scoring functions.

Supported Scoring Functions
---------------------------
Vina:
    Score derived from the AutoDock Vina scoring function, estimating binding affinity based on empirical energy terms.

RFScore:
    Machine-learning-based scoring function trained on protein–ligand complexes. It evaluates interaction patterns using different model versions.

NNScore:
    Neural network-based scoring function that predicts binding affinity from structural features.

PLECScore:
    Descriptor based on protein–ligand extended connectivity fingerprints capturing interaction patterns in the binding site.

Workflow
--------
1. Input one or more sets of docked small molecules.
2. Extract ligand poses associated with a receptor protein.
3. Select a descriptor type (Vina, RFScore, NNScore, or PLECScore).
4. Optionally select scoring model variants when applicable.
5. Compute interaction descriptors for each ligand–receptor complex.
6. Aggregate results into an annotated molecular dataset.

Output
------
- Set of docked small molecules annotated with computed ODDT descriptors.
- Tabular results containing scoring values per ligand.
- Updated molecular dataset compatible with downstream analysis workflows.

Use Cases
---------
- Evaluation of docking campaigns in virtual screening
- Ranking of ligand poses based on predicted binding affinity
- Post-processing of docking results for lead selection
- Comparative assessment of scoring functions

Notes
-----
- Results depend on docking pose quality and receptor preparation.
- Different scoring functions may produce non-comparable score distributions.
- RFScore and NNScore may vary depending on model version and training set.
"""
    _label = 'ODDT descriptors docking'
    _descChoices = ['Vina', 'RFScore', 'NNScore', 'PLECScore']#, 'Fingerprint']
    _vChoices = ['1', '2', '3']
    #_fpChoices = ['PLEC', 'SimpleInteractionFingerprint', 'ECFP', 'SPLIF']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputSmallMoleculesSets', params.MultiPointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Docked Small Molecules: ",
                       help='Select the docked molecules to be scored')

        group = form.addGroup('Descriptors')
        group.addParam('descChoice', params.EnumParam, default=0, label='Descriptor to calculate: ',
                       choices=self._descChoices,
                       help="Descriptor to calculate from the ODDT")
        group.addParam('vChoice', params.EnumParam, default=0, label='RFScore version: ',
                       choices=self._vChoices, condition='descChoice == 1',
                       help="Fingerprint to calculate from the ODDT")
        # group.addParam('fpChoice', params.EnumParam, default=0, label='Fingerprint to calculate: ',
        #                choices=self._fpChoices, condition='descChoice == 4',
        #                help="Fingerprint to calculate from the ODDT")

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('descriptorStep')
        self._insertFunctionStep('createOutputStep')

    def descriptorStep(self):
        mols = self.getInputMoleculesList()
        receptorFile = self.inputSmallMoleculesSets[0].get().getProteinFile()
        results = self.describeDockings(mols, receptorFile)

    def createOutputStep(self):
        usedIds, i = [], 1
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMoleculesSets[0].get(), self._getPath(), copyInfo=True)
        setattr(newMols, "_oddtDescriptors", String(self._getExtraPath('inputParams.txt')))
        for mol in self.getInputMoleculesList():
            #Id management
            if mol.getObjId() in usedIds:
                newId = i
                i += 1
            else:
                newId = mol.getObjId()
                i = newId + 1
            usedIds.append(newId)
            mol.setObjId(newId)
            newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        for molSet in self.inputSmallMoleculesSets:
            molSet = molSet.get()
            if not molSet.isDocked():
                validations += ['{} is not docked yet'.format(molSet)]

        return validations

    def _warnings(self):
        warnings = []
        if len(self.inputSmallMoleculesSets) > 1:
            warnings.append("The ids of the molecules will be rewritten so they don't collide with those of the "
                            "other input sets")
        protFile = ''
        for molSet in self.inputSmallMoleculesSets:
            molSet = molSet.get()
            if protFile == '':
                protFile = os.path.abspath(molSet.getProteinFile())
            elif not protFile == os.path.abspath(molSet.getProteinFile()):
                validations += ['These sets of small molecules are docked to different proteins.\n'
                                'Be aware that their results might not be comparable']

        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getInputMoleculesList(self):
        mols = []
        inputSets = fillEmptyAttributes(self.inputSmallMoleculesSets)
        for molSet in inputSets:
          for mol in molSet.get():
              mols.append(mol.clone())
        return mols

    def describeDockings(self, molsScipion, receptorFile):
        paramsPath = os.path.abspath(self._getExtraPath('inputParams.txt'))
        self.writeParamsFile(paramsPath, molsScipion, receptorFile)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    def writeParamsFile(self, paramsFile, molsScipion, recFile):
        molFiles = []
        for mol in molsScipion:
          molFiles.append(os.path.abspath(mol.getPoseFile()))

        with open(paramsFile, 'w') as f:
            f.write('outputPath: results.tsv\n')
            f.write('descriptor: {}\n'.format(self.getDescriptor()))
            if self.getDescriptor() == 'RFScore':
                f.write('version: {}\n'.format(self.getEnumText('vChoice')))
            # if self.getDescriptor() == 'Fingerprint':
            #     f.write('version: {}\n'.format(self.getEnumText('fpChoice')))

            f.write('receptorFile: {}\n'.format(os.path.abspath(recFile)))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

        return paramsFile
    
    def getDescriptor(self):
        function = self.getEnumText('descChoice')
        return function




