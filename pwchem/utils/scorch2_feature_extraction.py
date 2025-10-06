#!/usr/bin/env python3
"""
SCORCH2 Feature Extraction Tool

This script extracts features from protein-ligand complexes for use in SCORCH2. It processes protein-ligand pairs to extract:

1. Extended Connectivity Interaction Features (ECIF)
2. BINANA interaction features
3. Kier flexibility
4. RDKit descriptors

Usage:
    python scorch2_feature_extraction.py --protein-dir /path/to/proteins --ligand-dir /path/to/ligands --output-dir /path/to/output

Authors: Lin Chen
License: MIT
"""

import os
import pandas as pd
import kier, ecif, binana
from openbabel import pybel
from rdkit import Chem, RDLogger
import json
from functools import partial
from binana import PDB
from multiprocessing import Pool, TimeoutError
from tqdm import tqdm
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Descriptors
import glob
import argparse
import sys


# Mute all RDKit warnings
RDLogger.logger().setLevel(RDLogger.ERROR)

ob_log_handler = pybel.ob.OBMessageHandler()
ob_log_handler.SetOutputLevel(0)
pybel.ob.obErrorLog.StopLogging()


def calculate_ecifs(ligand_pdbqt_block, receptor_content):
    """
    Calculate Extended Connectivity Interaction Features (ECIF) for a ligand-receptor pair.
    
    Args:
        ligand_pdbqt_block: Ligand structure in PDBQT format
        receptor_content: Receptor structure path or content
        
    Returns:
        DataFrame containing ECIF features
    """
    try:
        ECIF_data = ecif.GetECIF(receptor_content, ligand_pdbqt_block, distance_cutoff=6.0)
        ECIFHeaders = [header.replace(';', '') for header in ecif.PossibleECIF]
        ECIF_data = dict(zip(ECIFHeaders, ECIF_data))
        return pd.DataFrame(ECIF_data, index=[0])
    except Exception as e:
        print(f"Error calculating ECIF features: {e}")
        # Return empty DataFrame with expected columns
        return pd.DataFrame({header.replace(';', ''): [0] for header in ecif.PossibleECIF})

def kier_flexibility(ligand_pdbqt_block):
    """
    Calculate Kier flexibility index and other rdkit molecular descriptors.
    
    Args:
        ligand_pdbqt_block: Ligand structure in PDBQT format
        
    Returns:
        Tuple containing (Kier flexibility index, dictionary of required rdkit descriptors)
    """
    # List of rdkit descriptors to calculate
    invariant_rdkit_descriptors = [
        'HeavyAtomMolWt','HeavyAtomCount','NumRotatableBonds', 'RingCount', 'NumAromaticRings', 'NumAliphaticRings',
        'NumSaturatedRings', 'NumAromaticHeterocycles',
        'NumAliphaticHeterocycles', 'NumSaturatedHeterocycles', 'FractionCSP3',
        'Chi0v', 'Chi1v', 'Chi2v', 'Chi3v', 'Chi4v', 'BalabanJ', 'BertzCT',
        'HallKierAlpha', 'Kappa1', 'Kappa2', 'Kappa3',
        'PEOE_VSA1', 'PEOE_VSA2', 'PEOE_VSA3', 'PEOE_VSA4', 'PEOE_VSA5',
        'PEOE_VSA6', 'PEOE_VSA7', 'PEOE_VSA8', 'PEOE_VSA9', 'PEOE_VSA10',
        'PEOE_VSA11', 'PEOE_VSA12', 'PEOE_VSA13', 'PEOE_VSA14',
        'SMR_VSA1', 'SMR_VSA2', 'SMR_VSA3', 'SMR_VSA4', 'SMR_VSA5',
        'SMR_VSA6', 'SMR_VSA7', 'SMR_VSA8', 'SMR_VSA9', 'SMR_VSA10',
        'SlogP_VSA1', 'SlogP_VSA2', 'SlogP_VSA3', 'SlogP_VSA4', 'SlogP_VSA5',
        'SlogP_VSA6', 'SlogP_VSA7', 'SlogP_VSA8', 'SlogP_VSA9', 'SlogP_VSA10',
        'SlogP_VSA11', 'SlogP_VSA12',
        'TPSA'
    ]

    try:
        # Get all available descriptor names
        descriptor_names = [desc[0] for desc in Descriptors._descList]

        # Prepare molecule
        mol = kier.SmilePrep(ligand_pdbqt_block)
        mol.GetRingInfo()
        mol_without_H = Chem.RemoveHs(mol)

        # Calculate all descriptors
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptor_names)
        descriptors = calculator.CalcDescriptors(mol_without_H)

        # Filter to get only the descriptors we need
        features = {}
        for name, value in zip(descriptor_names, descriptors):
            if name in invariant_rdkit_descriptors:
                features[name] = value

        # Return Kier flexibility and other descriptors
        return kier.CalculateFlexibility(mol), features
        
    except Exception as e:
        print(f"Error calculating Kier flexibility: {e}")
        # Return default values
        return 0, {name: 0 for name in invariant_rdkit_descriptors}

def run_binana(ligand_pdbqt_block, receptor_content):
    """
    Calculate BINANA (BINding ANAlyzer) features for a ligand-receptor pair.
    
    Args:
        ligand_pdbqt_block: Ligand structure in PDBQT format
        receptor_content: Receptor structure from binana PDB object
        
    Returns:
        Dictionary containing BINANA interaction features
    """
    binana_features = {}
    main_binana_out = binana.Binana(ligand_pdbqt_block, receptor_content).out

    # define the features we want
    keep_closest_contacts = ["2.5 (HD, OA)",
                            "2.5 (HD, HD)",
                            "2.5 (HD, N)",
                            "2.5 (C, HD)",
                            "2.5 (OA, ZN)",
                            "2.5 (HD, ZN)",
                            "2.5 (A, HD)"]

    keep_close_contacts = ["4.0 (C, C)",
                           "4.0 (HD, OA)",
                           "4.0 (C, HD)",
                           "4.0 (C, N)",
                           "4.0 (A, C)",
                           "4.0 (A, OA)",
                           "4.0 (N, OA)",
                           "4.0 (A, N)",
                           "4.0 (HD, N)",
                           "4.0 (HD, HD)",
                           "4.0 (A, HD)",
                           "4.0 (OA, OA)",
                           "4.0 (C, OA)",
                           "4.0 (N, N)",
                           "4.0 (C, SA)",
                           "4.0 (HD, SA)",
                           "4.0 (OA, SA)",
                           "4.0 (N, SA)",
                           "4.0 (A, A)",
                           "4.0 (HD, S)",
                           "4.0 (S, ZN)",
                           "4.0 (N, ZN)",
                           "4.0 (HD, ZN)",
                           "4.0 (A, SA)",
                           "4.0 (OA, ZN)",
                           "4.0 (C, ZN)",
                           "4.0 (C, NA)",
                           "4.0 (NA, OA)",
                           "4.0 (HD, NA)",
                           "4.0 (N, NA)",
                           "4.0 (A, NA)",
                           "4.0 (BR, C)",
                           "4.0 (HD, P)",
                           "4.0 (F, N)",
                           "4.0 (F, HD)",
                           "4.0 (C, CL)",
                           "4.0 (CL, HD)"]

    keep_ligand_atoms = ["LA N",
                         "LA HD"]

    keep_elsums = ["ElSum (C, C)",
                   "ElSum (HD, OA)",
                   "ElSum (C, HD)",
                   "ElSum (C, N)",
                   "ElSum (A, C)",
                   "ElSum (A, OA)",
                   "ElSum (N, OA)",
                   "ElSum (A, N)",
                   "ElSum (HD, HD)",
                   "ElSum (A, HD)",
                   "ElSum (OA, OA)",
                   "ElSum (C, OA)",
                   "ElSum (N, N)",
                   "ElSum (C, SA)",
                   "ElSum (HD, SA)",
                   "ElSum (OA, SA)",
                   "ElSum (N, SA)",
                   "ElSum (A, A)",
                   "ElSum (N, S)",
                   "ElSum (HD, S)",
                   "ElSum (OA, S)",
                   "ElSum (A, SA)",
                   "ElSum (C, NA)",
                   "ElSum (NA, OA)",
                   "ElSum (HD, NA)",
                   "ElSum (N, NA)",
                   "ElSum (A, NA)",
                   "ElSum (BR, C)",
                   "ElSum (HD, P)",
                   "ElSum (OA, P)",
                   "ElSum (N, P)",
                   "ElSum (C, F)",
                   "ElSum (F, N)",
                   "ElSum (A, F)",
                   "ElSum (CL, OA)",
                   "ElSum (C, CL)",
                   "ElSum (CL, N)",
                   "ElSum (A, CL)"]

    # add closest contacts to binana_features dict
    for contact in keep_closest_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ', '_')
        binana_features[contact] = main_binana_out['closest'].get(binana_name)

    # add close contacts to binana_features dict
    for contact in keep_close_contacts:
        binana_name = contact.split('(')[-1].split(')')[0].replace(', ', '_')
        binana_features[contact] = main_binana_out['close'].get(binana_name)

    # add ligand atoms to binana_features dict as binary tallies
    for atom in keep_ligand_atoms:
        binana_name = atom.split()[-1]
        if main_binana_out['ligand_atoms'].get(binana_name) is None:
            binana_features[atom] = 0
        else:
            binana_features[atom] = 1

    # add electrostatics to binana_features dict
    for elsum in keep_elsums:
        binana_name = elsum.split('(')[-1].split(')')[0].replace(', ', '_')
        binana_features[elsum] = main_binana_out['elsums'].get(binana_name)

    # add active site flexibility features to binana_features
    binana_features["BPF ALPHA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_ALPHA")
    binana_features["BPF ALPHA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_ALPHA")
    binana_features["BPF BETA SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_BETA")
    binana_features["BPF BETA BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_BETA")
    binana_features["BPF OTHER SIDECHAIN"] = main_binana_out['bpfs'].get("SIDECHAIN_OTHER")
    binana_features["BPF OTHER BACKBONE"] = main_binana_out['bpfs'].get("BACKBONE_OTHER")

    # add hydrophobic features to binana_features
    binana_features["HC ALPHA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_ALPHA")
    binana_features["HC ALPHA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_ALPHA")
    binana_features["HC BETA SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_BETA")
    binana_features["HC BETA BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_BETA")
    binana_features["HC OTHER SIDECHAIN"] = main_binana_out['hydrophobics'].get("SIDECHAIN_OTHER")
    binana_features["HC OTHER BACKBONE"] = main_binana_out['hydrophobics'].get("BACKBONE_OTHER")

    # add hydrogen bond features to binana_features
    binana_features["HB ALPHA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_ALPHA")
    binana_features["HB BETA SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE LIGAND"] = main_binana_out['hbonds'].get("HDONOR_LIGAND_BACKBONE_OTHER")
    binana_features["HB ALPHA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_ALPHA")
    binana_features["HB ALPHA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_ALPHA")
    binana_features["HB BETA SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_BETA")
    binana_features["HB BETA BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_BETA")
    binana_features["HB OTHER SIDECHAIN RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_OTHER")
    binana_features["HB OTHER BACKBONE RECEPTOR"] = main_binana_out['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_OTHER")

    # add salt bridge features to binana_features
    binana_features["SB ALPHA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binana_features["SB BETA"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_BETA")
    binana_features["SB OTHER"] = main_binana_out['salt_bridges'].get("SALT-BRIDGE_OTHER")

    # add aromatic stacking features to binana_features
    binana_features["piStack ALPHA"] = main_binana_out['stacking'].get("STACKING ALPHA")
    binana_features["piStack BETA"] = main_binana_out['stacking'].get("STACKING BETA")
    binana_features["piStack OTHER"] = main_binana_out['stacking'].get("STACKING OTHER")
    binana_features["tStack ALPHA"] = main_binana_out['t_stacking'].get("T-SHAPED_ALPHA")
    binana_features["tStack BETA"] = main_binana_out['t_stacking'].get("T-SHAPED_BETA")
    binana_features["tStack OTHER"] = main_binana_out['t_stacking'].get("T-SHAPED_OTHER")

    # add cation pi features to binana_features
    binana_features["catPi BETA LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binana_features["catPi OTHER LIGAND"] = main_binana_out['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    # add rotatable bond count to binana features
    binana_features["nRot"] = main_binana_out['nrot']

    return binana_features

def prune_df_headers(df):
    """
    Filter DataFrame to include only the required feature columns defined in SC1_features.json.
    
    Args:
        df: DataFrame with calculated features
        
    Returns:
        DataFrame with only the required columns
    """
    try:
        # Get path to the features JSON file (in the same directory as this script)
        script_dir = os.path.dirname(os.path.abspath(__file__))
        json_path = os.path.join(script_dir, 'SC1_features.json')
        
        with open(json_path, 'r') as f:
            reference_headers = json.load(f)
            headers_58 = reference_headers.get('492_models_58')
            
            # Check which columns exist in the DataFrame
            missing_cols = [col for col in headers_58 if col not in df.columns]
            if missing_cols:
                print(f"Warning: Missing {len(missing_cols)} columns: {', '.join(missing_cols[:5])}{'...' if len(missing_cols) > 5 else ''}")
                
            # Only include columns that exist in the DataFrame
            available_cols = [col for col in headers_58 if col in df.columns]
            return df[available_cols]
    except Exception as e:
        print(f"Error pruning DataFrame headers: {e}")
        # Return original DataFrame if there's an error
        return df




def process_molecule(molecule, ligand_path, pdbid, protein_path):
    """
    Process a single ligand molecule file and extract features from all poses.

    Args:
        molecule: Filename of the ligand file
        ligand_path: Path to the directory containing ligand files
        pdbid: PDB ID being processed
        protein_path: Path to the protein file

    Returns:
        DataFrame containing features for all poses in the ligand file
    """
    try:
        # Read the ligand file
        ligand_file = os.path.join(ligand_path, molecule)
        with open(ligand_file, 'r') as f:
            lig_text = f.read()

        # Split into individual poses (models)
        lig_poses = lig_text.split('MODEL')
        results = []

        # Process each pose
        for i, pose in enumerate(lig_poses):
            try:
                # Clean up the pose content
                lines = pose.split('\n')
                clean_lines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]

                # Skip if not enough content
                if len(clean_lines) < 3:
                    continue

                # Join cleaned lines back into a string
                pose_text = '\n'.join(clean_lines)

                # Calculate Kier flexibility and RDKit descriptors
                k, rdkit_descriptors = kier_flexibility(pose_text)
                entropy_df = pd.DataFrame([rdkit_descriptors])

                # Calculate BINANA features
                binana_features = run_binana(clean_lines, global_protein_object)
                binana_df = pd.DataFrame([binana_features])

                # Calculate ECIF features
                ecif_df = calculate_ecifs(pose_text, protein_path)

                # Combine all features
                df = pd.concat([ecif_df, binana_df], axis=1)
                df['Kier Flexibility'] = k

                try:
                    # Prune to required columns and add identifier
                    pruned_df = prune_df_headers(df)
                    combined_df = pd.concat([entropy_df, pruned_df], axis=1)
                    combined_df['Id'] = molecule
                    results.append(combined_df)
                except Exception as e:
                    print(f"Error in pruning/combining dataframes for {molecule}: {e}")
                    # Create a basic fallback dataframe to avoid losing computation
                    basic_df = pd.concat([entropy_df, df], axis=1)
                    basic_df['Id'] = molecule
                    results.append(basic_df)

            except Exception as e:
                print(f"Error processing pose {i} in {molecule}: {e}")
                continue

        # Combine results from all poses
        if results:
            try:
                return pd.concat(results, ignore_index=True)
            except Exception as e:
                print(f"Error concatenating results for {molecule}: {e}")
                # If concat fails, return the first result (better than nothing)
                if len(results) > 0:
                    return results[0]
                return pd.DataFrame()
        else:
            return pd.DataFrame()

    except Exception as e:
        print(f"Error processing molecule {molecule}: {e}")
        return pd.DataFrame()



def process_pdbid(pdbid, protein_base_path, molecule_path, des_path, num_cores=None):
    """
    Process a single PDB ID by extracting features from complex.

    Args:
        pdbid: The PDB ID to process
        protein_base_path: Path to directory containing protein PDBQT files
        molecule_path: Path to directory containing molecule PDBQT files
        des_path: Directory to save output files
        num_cores: Number of CPU cores to use (defaults to all available cores minus 1)
    """
    # Find the protein file
    protein_path = glob.glob(f'{protein_base_path}/{pdbid}*.pdbqt')
    if not protein_path:
        print(f'Protein file not found for {pdbid}')
        return
    protein_path = protein_path[0]

    # Check if output file already exists
    output_file = os.path.join(des_path, f'{pdbid}_protein_features.csv')
    if os.path.exists(output_file):
        print(f'PDBID {pdbid} Feature File exists - skipping')
        return

    # Check if molecule directory exists
    molecule_dir = os.path.join(molecule_path, pdbid)
    if not os.path.exists(molecule_dir):
        print(f'Molecules not found for {pdbid}')
        return
    molecules = os.listdir(molecule_dir)

    # Read protein content and start processing
    try:
        def init_worker(protein_path):
            global global_protein_object
            with open(protein_path, 'r') as f:
                protein_content = f.readlines()
            global_protein_object = PDB()
            global_protein_object.load_PDB(protein_path, protein_content)
            global_protein_object.assign_secondary_structure()

        # Determine number of processes to use
        if num_cores is None:
            processes = max(1, os.cpu_count() - 1)
        else:
            processes = min(num_cores, os.cpu_count())

        # Process molecules in parallel
        with Pool(processes=processes, initializer=init_worker, initargs=(protein_path,)) as pool:
            process_func = partial(
                process_molecule,
                ligand_path=molecule_dir,
                pdbid=pdbid,
                protein_path=protein_path  # keep for ECIF
            )
            futures = [pool.apply_async(process_func, (molecule,)) for molecule in molecules]

            results = []
            for i, future in enumerate(tqdm(futures, desc=f"Processing {pdbid} molecules", leave=False)):
                try:
                    result = future.get(timeout=5)  # 5-minute timeout per molecule
                    if not result.empty:
                        results.append(result)
                except TimeoutError:
                    print(f"Processing molecule {i} for {pdbid} timed out")
                except Exception as e:
                    print(f"Error processing {pdbid} molecule {i}: {e}")

            # Combine results and save to file
            if results:
                try:
                    total = pd.concat(results, ignore_index=True)

                    if not total.empty:
                        os.makedirs(des_path, exist_ok=True)
                        total.to_csv(output_file, index=False)
                        print(f"Saved features for {pdbid} ({len(total)} poses)")
                except Exception as e:
                    print(f"Error saving results for {pdbid}: {e}")
            else:
                print(f"No valid results for {pdbid}")

    except Exception as e:
        print(f'Error processing {pdbid}: {e}')

def main(args):
    """
    Main execution function that processes features for protein-ligand pairs.

    Args:
        args: Command line arguments
    """
    des_path = args.output_dir
    protein_base_path = args.protein_dir
    molecule_path = args.ligand_dir
    num_cores = args.num_cores

    # Create output directory if it doesn't exist
    os.makedirs(des_path, exist_ok=True)

    # Get list of PDB IDs
    if args.pdbids:
        pdbids = args.pdbids.split(',')
        print(f"Processing {len(pdbids)} specified PDB IDs")
    else:
        # Extract PDB IDs from filenames in protein directory
        pdbids = [i.split("_")[0] for i in os.listdir(protein_base_path)]
        print(f"Found {len(pdbids)} PDB IDs in {protein_base_path}")

    # Process each PDB ID with progress bar
    for pdbid in tqdm(pdbids, desc="Processing PDB structures"):
        process_pdbid(pdbid, protein_base_path, molecule_path, des_path, num_cores)




if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description="SCORCH2 Feature Extraction Tool",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # Required arguments
    parser.add_argument('--output-dir',
                      dest='output_dir',
                      type=str,
                      required=True,
                      help="Directory to save feature CSV files")

    parser.add_argument('--protein-dir',
                      dest='protein_dir',
                      type=str,
                      required=True,
                      help="Directory containing protein PDBQT files")

    parser.add_argument('--ligand-dir',
                      dest='ligand_dir',
                      type=str,
                      required=True,
                      help="Directory containing ligand PDBQT files")

    # Optional arguments
    parser.add_argument('--pdbids',
                      type=str,
                      default=None,
                      help="Comma-separated list of PDB IDs to process (if not specified, all PDB IDs in protein directory will be processed)")

    parser.add_argument('--num-cores',
                      dest='num_cores',
                      type=int,
                      default=None,
                      help="Number of CPU cores to use (default: all available cores minus 1)")

    parser.add_argument('--verbose',
                      action='store_true',
                      help="Enable verbose output")
    
    # Parse arguments and run main function
    args = parser.parse_args()
    
    # Print banner
    print("======================================")
    print("SCORCH2 Feature Extraction Tool")
    print("======================================")
    print(f"Protein directory: {args.protein_dir}")
    print(f"Ligand directory: {args.ligand_dir}")
    print(f"Output directory: {args.output_dir}")
    
    try:
        main(args)
        print("\nFeature extraction completed successfully")
    except Exception as e:
        print(f"\nError during feature extraction: {e}")
        sys.exit(1)
