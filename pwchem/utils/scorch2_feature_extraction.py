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

obLogHandler = pybel.ob.OBMessageHandler()
obLogHandler.SetOutputLevel(0)
pybel.ob.obErrorLog.StopLogging()


def calculateEcifs(ligandPdbtBlock, receptorContent):
    """
    Calculate Extended Connectivity Interaction Features (ECIF) for a ligand-receptor pair.
    
    Args:
        ligandPdbtBlock: Ligand structure in PDBQT format
        receptorContent: Receptor structure path or content
        
    Returns:
        DataFrame containing ECIF features
    """
    try:
        ECIFData = ecif.GetECIF(receptorContent, ligandPdbtBlock, distance_cutoff=6.0)
        ECIFHeaders = [header.replace(';', '') for header in ecif.PossibleECIF]
        ECIFData = dict(zip(ECIFHeaders, ECIFData))
        return pd.DataFrame(ECIFData, index=[0])
    except Exception as e:
        print(f"Error calculating ECIF features: {e}")
        # Return empty DataFrame with expected columns
        return pd.DataFrame({header.replace(';', ''): [0] for header in ecif.PossibleECIF})

def kierFlexibility(ligandPdbtBlock):
    """
    Calculate Kier flexibility index and other rdkit molecular descriptors.
    
    Args:
        ligandPdbtBlock: Ligand structure in PDBQT format
        
    Returns:
        Tuple containing (Kier flexibility index, dictionary of required rdkit descriptors)
    """
    # List of rdkit descriptors to calculate
    invariantRdkitDescriptors = [
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
        descriptorNames = [desc[0] for desc in Descriptors._descList]

        # Prepare molecule
        mol = kier.SmilePrep(ligandPdbtBlock)
        mol.GetRingInfo()
        molWithoutH = Chem.RemoveHs(mol)

        # Calculate all descriptors
        calculator = MoleculeDescriptors.MolecularDescriptorCalculator(descriptorNames)
        descriptors = calculator.CalcDescriptors(molWithoutH)

        # Filter to get only the descriptors we need
        features = {}
        for name, value in zip(descriptorNames, descriptors):
            if name in invariantRdkitDescriptors:
                features[name] = value

        # Return Kier flexibility and other descriptors
        return kier.CalculateFlexibility(mol), features
        
    except Exception as e:
        print(f"Error calculating Kier flexibility: {e}")
        # Return default values
        return 0, {name: 0 for name in invariantRdkitDescriptors}

def runBinana(ligandPdbtBlock, receptorContent):
    """
    Calculate BINANA (BINding ANAlyzer) features for a ligand-receptor pair.
    
    Args:
        ligandPdbtBlock: Ligand structure in PDBQT format
        receptorContent: Receptor structure from binana PDB object
        
    Returns:
        Dictionary containing BINANA interaction features
    """
    binanaFeatures = {}
    mainBinanaOut = binana.Binana(ligandPdbtBlock, receptorContent).out

    # define the features we want
    keepClosestContacts = ["2.5 (HD, OA)",
                            "2.5 (HD, HD)",
                            "2.5 (HD, N)",
                            "2.5 (C, HD)",
                            "2.5 (OA, ZN)",
                            "2.5 (HD, ZN)",
                            "2.5 (A, HD)"]

    keepCloseContacts = ["4.0 (C, C)",
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

    keepLigandAtoms = ["LA N",
                         "LA HD"]

    keepElsums = ["ElSum (C, C)",
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

    # add closest contacts to binanaFeatures dict
    for contact in keepClosestContacts:
        binanaName = contact.split('(')[-1].split(')')[0].replace(', ', '_')
        binanaFeatures[contact] = mainBinanaOut['closest'].get(binanaName)

    # add close contacts to binanaFeatures dict
    for contact in keepCloseContacts:
        binanaName = contact.split('(')[-1].split(')')[0].replace(', ', '_')
        binanaFeatures[contact] = mainBinanaOut['close'].get(binanaName)

    # add ligand atoms to binanaFeatures dict as binary tallies
    for atom in keepLigandAtoms:
        binanaName = atom.split()[-1]
        if mainBinanaOut['ligand_atoms'].get(binanaName) is None:
            binanaFeatures[atom] = 0
        else:
            binanaFeatures[atom] = 1

    # add electrostatics to binanaFeatures dict
    for elsum in keepElsums:
        binanaName = elsum.split('(')[-1].split(')')[0].replace(', ', '_')
        binanaFeatures[elsum] = mainBinanaOut['elsums'].get(binanaName)

    # add active site flexibility features to binanaFeatures
    binanaFeatures["BPF ALPHA SIDECHAIN"] = mainBinanaOut['bpfs'].get("SIDECHAIN_ALPHA")
    binanaFeatures["BPF ALPHA BACKBONE"] = mainBinanaOut['bpfs'].get("BACKBONE_ALPHA")
    binanaFeatures["BPF BETA SIDECHAIN"] = mainBinanaOut['bpfs'].get("SIDECHAIN_BETA")
    binanaFeatures["BPF BETA BACKBONE"] = mainBinanaOut['bpfs'].get("BACKBONE_BETA")
    binanaFeatures["BPF OTHER SIDECHAIN"] = mainBinanaOut['bpfs'].get("SIDECHAIN_OTHER")
    binanaFeatures["BPF OTHER BACKBONE"] = mainBinanaOut['bpfs'].get("BACKBONE_OTHER")

    # add hydrophobic features to binanaFeatures
    binanaFeatures["HC ALPHA SIDECHAIN"] = mainBinanaOut['hydrophobics'].get("SIDECHAIN_ALPHA")
    binanaFeatures["HC ALPHA BACKBONE"] = mainBinanaOut['hydrophobics'].get("BACKBONE_ALPHA")
    binanaFeatures["HC BETA SIDECHAIN"] = mainBinanaOut['hydrophobics'].get("SIDECHAIN_BETA")
    binanaFeatures["HC BETA BACKBONE"] = mainBinanaOut['hydrophobics'].get("BACKBONE_BETA")
    binanaFeatures["HC OTHER SIDECHAIN"] = mainBinanaOut['hydrophobics'].get("SIDECHAIN_OTHER")
    binanaFeatures["HC OTHER BACKBONE"] = mainBinanaOut['hydrophobics'].get("BACKBONE_OTHER")

    # add hydrogen bond features to binanaFeatures
    binanaFeatures["HB ALPHA SIDECHAIN LIGAND"] = mainBinanaOut['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_ALPHA")
    binanaFeatures["HB BETA SIDECHAIN LIGAND"] = mainBinanaOut['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_BETA")
    binanaFeatures["HB BETA BACKBONE LIGAND"] = mainBinanaOut['hbonds'].get("HDONOR_LIGAND_BACKBONE_BETA")
    binanaFeatures["HB OTHER SIDECHAIN LIGAND"] = mainBinanaOut['hbonds'].get("HDONOR_LIGAND_SIDECHAIN_OTHER")
    binanaFeatures["HB OTHER BACKBONE LIGAND"] = mainBinanaOut['hbonds'].get("HDONOR_LIGAND_BACKBONE_OTHER")
    binanaFeatures["HB ALPHA SIDECHAIN RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_ALPHA")
    binanaFeatures["HB ALPHA BACKBONE RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_ALPHA")
    binanaFeatures["HB BETA SIDECHAIN RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_BETA")
    binanaFeatures["HB BETA BACKBONE RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_BETA")
    binanaFeatures["HB OTHER SIDECHAIN RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_SIDECHAIN_OTHER")
    binanaFeatures["HB OTHER BACKBONE RECEPTOR"] = mainBinanaOut['hbonds'].get("HDONOR_RECEPTOR_BACKBONE_OTHER")

    # add salt bridge features to binanaFeatures
    binanaFeatures["SB ALPHA"] = mainBinanaOut['salt_bridges'].get("SALT-BRIDGE_ALPHA")
    binanaFeatures["SB BETA"] = mainBinanaOut['salt_bridges'].get("SALT-BRIDGE_BETA")
    binanaFeatures["SB OTHER"] = mainBinanaOut['salt_bridges'].get("SALT-BRIDGE_OTHER")

    # add aromatic stacking features to binanaFeatures
    binanaFeatures["piStack ALPHA"] = mainBinanaOut['stacking'].get("STACKING ALPHA")
    binanaFeatures["piStack BETA"] = mainBinanaOut['stacking'].get("STACKING BETA")
    binanaFeatures["piStack OTHER"] = mainBinanaOut['stacking'].get("STACKING OTHER")
    binanaFeatures["tStack ALPHA"] = mainBinanaOut['t_stacking'].get("T-SHAPED_ALPHA")
    binanaFeatures["tStack BETA"] = mainBinanaOut['t_stacking'].get("T-SHAPED_BETA")
    binanaFeatures["tStack OTHER"] = mainBinanaOut['t_stacking'].get("T-SHAPED_OTHER")

    # add cation pi features to binanaFeatures
    binanaFeatures["catPi BETA LIGAND"] = mainBinanaOut['pi_cation'].get("PI-CATION_LIGAND-CHARGED_BETA")
    binanaFeatures["catPi OTHER LIGAND"] = mainBinanaOut['pi_cation'].get("PI-CATION_LIGAND-CHARGED_OTHER")

    # add rotatable bond count to binana features
    binanaFeatures["nRot"] = mainBinanaOut['nrot']

    return binanaFeatures

def pruneDfHeaders(df):
    """
    Filter DataFrame to include only the required feature columns defined in SC1_features.json.
    
    Args:
        df: DataFrame with calculated features
        
    Returns:
        DataFrame with only the required columns
    """
    try:
        # Get path to the features JSON file (in the same directory as this script)
        scriptDir = os.path.dirname(os.path.abspath(__file__))
        jsonPath = os.path.join(scriptDir, 'SC1_features.json')
        
        with open(jsonPath, 'r') as f:
            referenceHeaders = json.load(f)
            headers58 = referenceHeaders.get('492_models_58')
            
            # Check which columns exist in the DataFrame
            missingCols = [col for col in headers58 if col not in df.columns]
            if missingCols:
                print(f"Warning: Missing {len(missingCols)} columns: {', '.join(missingCols[:5])}{'...' if len(missingCols) > 5 else ''}")
                
            # Only include columns that exist in the DataFrame
            availableCols = [col for col in headers58 if col in df.columns]
            return df[availableCols]
    except Exception as e:
        print(f"Error pruning DataFrame headers: {e}")
        # Return original DataFrame if there's an error
        return df




def processMolecule(molecule, ligandPath, pdbid, proteinPath):
    """
    Process a single ligand molecule file and extract features from all poses.

    Args:
        molecule: Filename of the ligand file
        ligandPath: Path to the directory containing ligand files
        pdbid: PDB ID being processed
        proteinPath: Path to the protein file

    Returns:
        DataFrame containing features for all poses in the ligand file
    """
    try:
        # Read the ligand file
        ligandFile = os.path.join(ligandPath, molecule)
        with open(ligandFile, 'r') as f:
            ligText = f.read()

        # Split into individual poses (models)
        ligPoses = ligText.split('MODEL')
        results = []

        # Process each pose
        for i, pose in enumerate(ligPoses):
            try:
                # Clean up the pose content
                lines = pose.split('\n')
                cleanLines = [line for line in lines if not line.strip().lstrip().isnumeric() and 'ENDMDL' not in line]

                # Skip if not enough content
                if len(cleanLines) < 3:
                    continue

                # Join cleaned lines back into a string
                poseText = '\n'.join(cleanLines)

                # Calculate Kier flexibility and RDKit descriptors
                k, rdkitDescriptors = kierFlexibility(poseText)
                entropyDf = pd.DataFrame([rdkitDescriptors])

                # Calculate BINANA features
                binanaFeatures = runBinana(cleanLines, globalProteinObject)
                binanaDf = pd.DataFrame([binanaFeatures])

                # Calculate ECIF features
                ecifDf = calculateEcifs(poseText, proteinPath)

                # Combine all features
                df = pd.concat([ecifDf, binanaDf], axis=1)
                df['Kier Flexibility'] = k

                try:
                    # Prune to required columns and add identifier
                    prunedDf = pruneDfHeaders(df)
                    combinedDf = pd.concat([entropyDf, prunedDf], axis=1)
                    combinedDf['Id'] = molecule
                    results.append(combinedDf)
                except Exception as e:
                    print(f"Error in pruning/combining dataframes for {molecule}: {e}")
                    # Create a basic fallback dataframe to avoid losing computation
                    basicDf = pd.concat([entropyDf, df], axis=1)
                    basicDf['Id'] = molecule
                    results.append(basicDf)

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



def processPdbid(pdbid, proteinBasePath, moleculePath, desPath, numCores=None):
    """
    Process a single PDB ID by extracting features from complex.

    Args:
        pdbid: The PDB ID to process
        proteinBasePath: Path to directory containing protein PDBQT files
        moleculePath: Path to directory containing molecule PDBQT files
        desPath: Directory to save output files
        numCores: Number of CPU cores to use (defaults to all available cores minus 1)
    """
    # Find the protein file
    proteinPath = glob.glob(f'{proteinBasePath}/{pdbid}*.pdbqt')
    if not proteinPath:
        print(f'Protein file not found for {pdbid}')
        return
    proteinPath = proteinPath[0]

    # Check if output file already exists
    outputFile = os.path.join(desPath, f'{pdbid}_protein_features.csv')
    if os.path.exists(outputFile):
        print(f'PDBID {pdbid} Feature File exists - skipping')
        return

    # Check if molecule directory exists
    moleculeDir = os.path.join(moleculePath, pdbid)
    if not os.path.exists(moleculeDir):
        print(f'Molecules not found for {pdbid}')
        return
    molecules = os.listdir(moleculeDir)

    # Read protein content and start processing
    try:
        def initWorker(proteinPath):
            global globalProteinObject
            with open(proteinPath, 'r') as f:
                proteinContent = f.readlines()
            globalProteinObject = PDB()
            globalProteinObject.load_PDB(proteinPath, proteinContent)
            globalProteinObject.assign_secondary_structure()

        # Determine number of processes to use
        if numCores is None:
            processes = max(1, os.cpu_count() - 1)
        else:
            processes = min(numCores, os.cpu_count())

        # Process molecules in parallel
        with Pool(processes=processes, initializer=initWorker, initargs=(proteinPath,)) as pool:
            processFunc = partial(
                processMolecule,
                ligandPath=moleculeDir,
                pdbid=pdbid,
                proteinPath=proteinPath  # keep for ECIF
            )
            futures = [pool.apply_async(processFunc, (molecule,)) for molecule in molecules]

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
                        os.makedirs(desPath, exist_ok=True)
                        total.to_csv(outputFile, index=False)
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
    desPath = args.output_dir
    proteinBasePath = args.protein_dir
    moleculePath = args.ligand_dir
    numCores = args.num_cores

    # Create output directory if it doesn't exist
    os.makedirs(desPath, exist_ok=True)

    # Get list of PDB IDs
    if args.pdbids:
        pdbids = args.pdbids.split(',')
        print(f"Processing {len(pdbids)} specified PDB IDs")
    else:
        # Extract PDB IDs from filenames in protein directory
        pdbids = [i.split("_")[0] for i in os.listdir(proteinBasePath)]
        print(f"Found {len(pdbids)} PDB IDs in {proteinBasePath}")

    # Process each PDB ID with progress bar
    for pdbid in tqdm(pdbids, desc="Processing PDB structures"):
        processPdbid(pdbid, proteinBasePath, moleculePath, desPath, numCores)




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
