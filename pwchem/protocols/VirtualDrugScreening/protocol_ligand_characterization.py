# **************************************************************************
# *
# * Authors: ...
# *
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

import numpy as np
import os
from rdkit import Chem
from rdkit.Chem import Descriptors

from pyworkflow.protocol import params
from pyworkflow.object import Float
from pwem.protocols import EMProtocol
from pwchem.objects import SetOfSmallMolecules

class ProtocolLigandCharacterization(EMProtocol):
    _label = 'Ligand characterization'

    def _defineParams(self, form):
        form.addSection(label='Parameters')
        form.addParam('inputSmallMolecules', params.PointerParam,
                      pointerClass='SetOfSmallMolecules', allowsNull=False,
                      label="Input Small Molecules",
                      help='Set of molecules to process.')
        # Nombrar todos los parámetros con booleans
        form.addParam('useConstitutional', params.BooleanParam, default=True,
                    label='Calculate constitutional descriptors')
        form.addParam('useElectronic', params.BooleanParam, default=True,
                    label='Calculate electronic descriptors')
        form.addParam('useTopological', params.BooleanParam, default=True,
                    label='Calculate topological descriptors')
        form.addParam('useGeometrical', params.BooleanParam, default=True,
                    label='Calculate geometrical descriptors')
        form.addParam('useRing', params.BooleanParam, default=True,
                    label='Calculate ring descriptors')
        form.addParam('useFragment', params.BooleanParam, default=True,
                    label='Calculate fragment descriptors')
        form.addParam('useOther', params.BooleanParam, default=True,
                    label='Calculate other descriptors')

    def _insertAllSteps(self):
        self._insertFunctionStep('runDescriptorCalc')
        self._insertFunctionStep('createOutputStep')

    def runDescriptorCalc(self):
        mol_set = self.inputSmallMolecules.get()

        # Get all available RDKit descriptors
        descriptor_list = Descriptors.descList
        self.header = [name for name, _ in descriptor_list]
        self.propertyDict={}

        self. descriptor_categories = {'other': ['MaxAbsEStateIndex', 'MaxEStateIndex', 'MinAbsEStateIndex', 'MinEStateIndex', 'qed', 'SPS', 'MaxPartialCharge', 'MinPartialCharge', 'MaxAbsPartialCharge', 'MinAbsPartialCharge', 'FpDensityMorgan1', 'FpDensityMorgan2', 'FpDensityMorgan3', 'BCUT2D_MWHI', 'BCUT2D_MWLOW', 'BCUT2D_CHGHI', 'BCUT2D_CHGLO', 'BCUT2D_LOGPHI', 'BCUT2D_LOGPLOW', 'BCUT2D_MRHI', 'BCUT2D_MRLOW', 'AvgIpc', 'BalabanJ', 'BertzCT', 'Ipc', 'PEOE_VSA1', 'PEOE_VSA10', 'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5', 'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'SMR_VSA1', 'SMR_VSA10', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5', 'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SlogP_VSA1', 'SlogP_VSA10', 'SlogP_VSA11', 'SlogP_VSA12', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5', 'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'EState_VSA1', 'EState_VSA10', 'EState_VSA11', 'EState_VSA2', 'EState_VSA3', 'EState_VSA4', 'EState_VSA5', 'EState_VSA6', 'EState_VSA7', 'EState_VSA8', 'EState_VSA9', 'VSA_EState1', 'VSA_EState10', 'VSA_EState2', 'VSA_EState3', 'VSA_EState4', 'VSA_EState5', 'VSA_EState6', 'VSA_EState7', 'VSA_EState8', 'VSA_EState9', 'FractionCSP3', 'NHOHCount', 'NOCount', 'NumAliphaticCarbocycles', 'NumAliphaticHeterocycles', 'NumAmideBonds', 'NumAromaticCarbocycles', 'NumAromaticHeterocycles', 'NumAtomStereoCenters', 'NumBridgeheadAtoms', 'NumHeteroatoms', 'NumHeterocycles', 'NumRotatableBonds', 'NumSaturatedCarbocycles', 'NumSaturatedHeterocycles', 'NumSaturatedRings', 'NumSpiroAtoms', 'NumUnspecifiedAtomStereoCenters', 'Phi'], 
                         'constitutional': ['MolWt', 'HeavyAtomMolWt', 'ExactMolWt', 'NumValenceElectrons', 'HeavyAtomCount', 'NumHAcceptors', 'NumHDonors'], 
                         'electronic': ['ExactMolWt', 'NumRadicalElectrons', 'MolLogP', 'MolMR'], 
                         'topological': ['Chi0', 'Chi0n', 'Chi0v', 'Chi1', 'Chi1n', 'Chi1v', 'Chi2n', 'Chi2v', 'Chi3n', 'Chi3v', 'Chi4n', 'Chi4v', 'HallKierAlpha', 'Kappa1', 'Kappa2', 'Kappa3', 'TPSA'], 
                         'geometrical': ['LabuteASA'], 
                         'ring': ['NumAliphaticRings', 'NumAromaticRings', 'RingCount'], 
                         'fragment': ['fr_Al_COO', 'fr_Al_OH', 'fr_Al_OH_noTert', 'fr_ArN', 'fr_Ar_COO', 'fr_Ar_N', 'fr_Ar_NH', 'fr_Ar_OH', 'fr_COO', 'fr_COO2', 'fr_C_O', 'fr_C_O_noCOO', 'fr_C_S', 'fr_HOCCN', 'fr_Imine', 'fr_NH0', 'fr_NH1', 'fr_NH2', 'fr_N_O', 'fr_Ndealkylation1', 'fr_Ndealkylation2', 'fr_Nhpyrrole', 'fr_SH', 'fr_aldehyde', 'fr_alkyl_carbamate', 'fr_alkyl_halide', 'fr_allylic_oxid', 'fr_amide', 'fr_amidine', 'fr_aniline', 'fr_aryl_methyl', 'fr_azide', 'fr_azo', 'fr_barbitur', 'fr_benzene', 'fr_benzodiazepine', 'fr_bicyclic', 'fr_diazo', 'fr_dihydropyridine', 'fr_epoxide', 'fr_ester', 'fr_ether', 'fr_furan', 'fr_guanido', 'fr_halogen', 'fr_hdrzine', 'fr_hdrzone', 'fr_imidazole', 'fr_imide', 'fr_isocyan', 'fr_isothiocyan', 'fr_ketone', 'fr_ketone_Topliss', 'fr_lactam', 'fr_lactone', 'fr_methoxy', 'fr_morpholine', 'fr_nitrile', 'fr_nitro', 'fr_nitro_arom', 'fr_nitro_arom_nonortho', 'fr_nitroso', 'fr_oxazole', 'fr_oxime', 'fr_para_hydroxylation', 'fr_phenol', 'fr_phenol_noOrthoHbond', 'fr_phos_acid', 'fr_phos_ester', 'fr_piperdine', 'fr_piperzine', 'fr_priamide', 'fr_prisulfonamd', 'fr_pyridine', 'fr_quatN', 'fr_sulfide', 'fr_sulfonamd', 'fr_sulfone', 'fr_term_acetylene', 'fr_tetrazole', 'fr_thiazole', 'fr_thiocyan', 'fr_thiophene', 'fr_unbrch_alkane', 'fr_urea']}

        for mol in mol_set:
            path = mol.getFileName()
            rdkit_mol = Chem.MolFromMolFile(path)
            if rdkit_mol is None:
                continue

            row = []
            for name, desc_func in descriptor_list:
                evaluate = False
                for category in self.descriptor_categories:
                    if getattr(self, f"use{category.capitalize()}") and name in self.descriptor_categories[category]:
                        evaluate = True
                        break

                if evaluate:
                    try:
                        value = desc_func(rdkit_mol)
                    except Exception:
                        value = np.nan
                else:
                    value = np.nan

                row.append(value)

            self.propertyDict[mol.molName.get()] = row


    def createOutputStep(self):
        newMols = SetOfSmallMolecules.createCopy(self.inputSmallMolecules.get(), self._getPath(), copyInfo=True)

        mols = self.inputSmallMolecules.get()
        for mol in mols:
            # NUEVO
            if mol.molName.get() not in self.propertyDict:
                continue
            # Hasta aquí
            row = self.propertyDict[mol.molName.get()]
            # for i in range(len(row)):
                # setattr(mol, f"_Property_{self.header[i]}", Float(row[i]))
                # eval(f"mol._Property_{self.header[i]}=pwobj.Float({row[i]})")
                # mol._fingerprintSimilarity = pwobj.String(filtered_molecules_dict[molFile])
            
            for i in range(len(row)):
                prop_name = self.header[i]
                value = row[i]
                # Buscar la categoría del descriptor
                category = 'uncategorized'
                for cat, prop_list in self.descriptor_categories.items():
                    if prop_name in prop_list:
                        category = cat
                        break
                attr_name = f"_Property_{category}_{prop_name}"
                setattr(mol, attr_name, Float(value))

            newMols.append(mol)

        newMols.updateMolClass()
        self._defineOutputs(outputSmallMolecules=newMols)

    def _summary(self):
        pass

    def getPluginCategory(self):
        return "Virtual Screening"
