# **************************************************************************
# *
# * Authors:    Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

# Virtual Drug Screening protocols
from pwchem.protocols.VirtualDrugScreening.protocol_import_smallMolecules import ProtChemImportSmallMolecules
from pwchem.protocols.VirtualDrugScreening.protocol_import_molecules_library import ProtChemImportMoleculesLibrary
from pwchem.protocols.VirtualDrugScreening.protocol_drawSmallMol import ProtDrawMolecules
from pwchem.protocols.VirtualDrugScreening.protocol_extract_ligands import ProtExtractLigands

from pwchem.protocols.VirtualDrugScreening.protocol_receptor_preparation import ProtChemPrepareReceptor
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_preparation import ProtChemOBabelPrepareLigands
from pwchem.protocols.VirtualDrugScreening.protocol_rdkit_ligand_preparation import ProtChemRDKitPrepareLigands

from pwchem.protocols.VirtualDrugScreening.protocol_PAINS_filter import ProtocolPainsRdkitFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_ADME_filter import ProtocolADMEFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_shape_filtering import ProtocolShapeDistancesFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_fingerprints import ProtocolFingerprintFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_filter import ProtocolGeneralLigandFiltering

from pwchem.protocols.VirtualDrugScreening.protocol_define_manual_structROIs import ProtDefineStructROIs
from pwchem.protocols.VirtualDrugScreening.protocol_consensus_structROIs import ProtocolConsensusStructROIs
from pwchem.protocols.VirtualDrugScreening.protocol_extract_interactingMols import ProtExtractInteractingMols

from pwchem.protocols.VirtualDrugScreening.protocol_consensus_docking import ProtocolConsensusDocking
from pwchem.protocols.VirtualDrugScreening.protocol_score_dockings import ProtocolScoreDocking
# from pwchem.protocols.VirtualDrugScreening.protocol_oddt_descriptors import ProtocolODDTDescriptors
from pwchem.protocols.VirtualDrugScreening.protocol_rmsd_dockings import ProtocolRMSDDocking
from pwchem.protocols.VirtualDrugScreening.protocol_define_contact_structROIs import ProtDefineContactStructROIs
from pwchem.protocols.VirtualDrugScreening.protocol_rank_docking_score import ProtocolRankDocking
from pwchem.protocols.VirtualDrugScreening.protocol_scores_correlation import ProtScoreCorrelation

from pwchem.protocols.VirtualDrugScreening.protocol_calculate_SASA import ProtCalculateSASA

from pwchem.protocols.VirtualDrugScreening.protocol_pharmacophore_generation import ProtocolPharmacophoreGeneration
from pwchem.protocols.VirtualDrugScreening.protocol_pharmacophore_modification import ProtocolPharmacophoreModification
from pwchem.protocols.VirtualDrugScreening.protocol_pharmacophore_filtering import ProtocolPharmacophoreFiltering

# Sequences and sequences ROIs protocols
from pwchem.protocols.Sequences.protocol_import_setOfSequences import ProtChemImportSetOfSequences
from pwchem.protocols.Sequences.protocol_define_sequences import ProtDefineSetOfSequences

from pwchem.protocols.Sequences.protocol_import_variants import ProtChemImportVariants
from pwchem.protocols.Sequences.protocol_generate_variants import ProtChemGenerateVariants

from pwchem.protocols.Sequences.protocol_define_sequence_roi import ProtDefineSeqROI
from pwchem.protocols.Sequences.protocol_import_sequence_roi import ProtImportSeqROI
from pwchem.protocols.Sequences.protocol_calculate_conservation import ProtSeqCalculateConservation
from pwchem.protocols.Sequences.protocol_extract_attribute_ROIs import ProtExtractSeqsROI
from pwchem.protocols.Sequences.protocol_operate_sequence_rois import ProtOperateSeqROI
from pwchem.protocols.Sequences.protocol_map_sequence_structure_ROIs import ProtMapSequenceROI
from pwchem.protocols.Sequences.protocol_map_attribute_to_sequence_ROIs import ProtMapAttributeToSeqROIs

from pwchem.protocols.Sequences.protocol_pairWise_alignment import ProtChemPairWiseAlignment
from pwchem.protocols.Sequences.protocol_multipleSequence_alignment import ProtChemMultipleSequenceAlignment

from pwchem.protocols.Sequences.protocol_define_multiepitope import ProtDefineMultiEpitope
from pwchem.protocols.Sequences.protocol_modify_multiepitope import ProtModifyMultiEpitope
from pwchem.protocols.Sequences.protocol_combineScores_sequence_roi import ProtCombineScoresSeqROI
from pwchem.protocols.Sequences.protocol_optimize_multiepitope import ProtOptimizeMultiEpitope


# Databases protocols
from pwchem.protocols.Databases.protocol_import_setOfDatabaseIDs import ProtChemImportSetOfDatabaseIDs
from pwchem.protocols.Databases.protocol_fetch_ligands import ProtocolLigandsFetching
from pwchem.protocols.Databases.protocol_smallMol_identify import ProtChemSmallMolIdentify
from pwchem.protocols.Databases.protocol_uniprot_crossref import ProtChemUniprotCrossRef
from pwchem.protocols.Databases.protocol_ZINC_filter import ProtChemZINCFilter

# General protocols
from pwchem.protocols.General.protocol_converter import ConvertStructures
from pwchem.protocols.General.protocol_add_attribute import ProtAddAttribute
from pwchem.protocols.General.protocol_operate_set import ProtChemOperateSet
from pwchem.protocols.General.protocol_operate_libraries import ProtocolOperateLibrary
from pwchem.protocols.General.protocol_export_csv import ProtChemExportCSV
from pwchem.protocols.General.protocol_pymol import ProtPymolOperate
from pwchem.protocols.General.protocol_ranx_fuse import ProtocolRANXFuse
from pwchem.protocols.General.protocol_operate_libraries import ProtocolOperateLibrary

# Molecular dynamics
from pwchem.protocols.MolecularDynamics.protocol_parametrize_ligand import ProtocolLigandParametrization

