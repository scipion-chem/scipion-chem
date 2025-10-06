#!/usr/bin/env python3
"""
SCORCH2 Rescoring Code

This script performs the complete SCORCH2 workflow including feature extraction,
normalization, and rescoring in a single streamlined process.

Usage:
    python scorch2_rescoring.py --protein-dir protein/ --ligand-dir molecule/ \
                               --sc2_ps_model sc2_ps.xgb --sc2_pb_model sc2_pb.xgb \
                               --output results.csv

Authors: Lin Chen
License: MIT
"""

import os
import argparse
import pandas as pd
import numpy as np
import xgboost as xgb
import pickle
import joblib
import subprocess
import sys
import re
from pathlib import Path
from tqdm import tqdm
from datetime import datetime
from typing import Dict, List, Tuple, Optional
import time
import select
import json
import threading
import queue


def run_feature_extraction(protein_dir: str, ligand_dir: str, output_dir: str, num_cores: int = None) -> bool:
    """
    Run feature extraction using the SCORCH2 feature extraction utility.
    
    Args:
        protein_dir: Directory containing protein PDBQT files
        ligand_dir: Directory containing ligand PDBQT files
        output_dir: Directory to save extracted features
        num_cores: Number of CPU cores to use (default: os.cpu_count()-1)
        
    Returns:
        True if successful, False otherwise
    """
    if num_cores is None:
        num_cores = max(1, os.cpu_count() - 1)
    
    print("Running feature extraction...")
    print(f"  Using {num_cores} CPU cores")
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    # Get list of protein targets to estimate total progress
    try:
        protein_files = [f for f in os.listdir(protein_dir) if f.endswith('.pdbqt')]
        pdbids = list(set([f.split("_")[0] for f in protein_files]))
        total_targets = len(pdbids)
        print(f"  Found {total_targets} protein targets")
    except Exception as e:
        print(f"Warning: Could not count targets for progress tracking: {e}")
        total_targets = None
    
    # Run feature extraction with progress monitoring
    cmd = [
        sys.executable, "utils/scorch2_feature_extraction.py",
        "--protein-dir", protein_dir,
        "--ligand-dir", ligand_dir, 
        "--output-dir", output_dir,
        "--num-cores", str(num_cores)
    ]
    
    try:
        # Initialize progress bar
        if total_targets:
            pbar = tqdm(total=total_targets, desc="Processing targets", unit="target")
        else:
            pbar = tqdm(desc="Feature extraction", unit="operation")
        
        print("Starting feature extraction subprocess...")
        
        # Initialize sets
        processed_targets = set()

        # Start the subprocess with unbuffered output
        env = {**os.environ, "PYTHONUNBUFFERED": "1"}
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                  text=True, universal_newlines=True, env=env)

        # Define line processing function
        def process_line(line: str, pbar, processed_targets, is_stderr: bool = False):
            prefix = "⚠️ Subprocess: " if is_stderr else ""
            line = line.strip()
            if not line:
                return

            # Saved features
            match = re.search(r"Saved features for (\w+)", line)
            if match:
                target = match.group(1)
                if target not in processed_targets:
                    processed_targets.add(target)
                    pbar.set_description(f"Completed: {target}")
                    pbar.update(1)
                return

            # Skipped features
            match = re.search(r"PDBID (\w+) Feature File exists", line)
            if match:
                target = match.group(1)
                if target not in processed_targets:
                    processed_targets.add(target)
                    pbar.set_description(f"Skipped: {target}")
                    pbar.update(1)
                return

            # Error messages (ignoring specific molecule errors)
            if any(kw in line.lower() for kw in ["error", "failed", "not found"]):
                if not (line.startswith("Error processing") and "molecule" in line.lower()):
                    pbar.write(f"{prefix}⚠️  {line}")
                return

            # General progress info
            if "processing" in line.lower() and "pdb ids" in line.lower():
                pbar.write(f"{prefix}📊 {line}")
                return

            # Optionally handle other lines if needed
            # pbar.write(f"{prefix}{line}")

        # Define reader functions
        def read_stdout(queue, pipe):
            for line in iter(pipe.readline, ''):
                queue.put(line)
            queue.put(None)

        def read_stderr(queue, pipe):
            for line in iter(pipe.readline, ''):
                queue.put(line)
            queue.put(None)

        # Start reader threads
        stdout_queue = queue.Queue()
        stderr_queue = queue.Queue()
        stdout_thread = threading.Thread(target=read_stdout, args=(stdout_queue, process.stdout), daemon=True)
        stderr_thread = threading.Thread(target=read_stderr, args=(stderr_queue, process.stderr), daemon=True)
        stdout_thread.start()
        stderr_thread.start()

        # Process output in real-time
        # Monitor until process completes
        while process.poll() is None:
            updated = False

            # Process stdout
            try:
                item = stdout_queue.get_nowait()
                if item is not None:
                    process_line(item, pbar, processed_targets, False)
                    updated = True
            except queue.Empty:
                pass

            # Process stderr
            try:
                item = stderr_queue.get_nowait()
                if item is not None:
                    process_line(item, pbar, processed_targets, True)
                    updated = True
            except queue.Empty:
                pass

            if not updated:
                time.sleep(0.1)

        # Process has completed - drain remaining output with timeout
        drain_start = time.time()
        while (time.time() - drain_start < 10) and (not stdout_queue.empty() or not stderr_queue.empty()):
            # Drain stdout
            try:
                item = stdout_queue.get(timeout=1.0)
                if item is not None:
                    process_line(item, pbar, processed_targets, False)
            except queue.Empty:
                pass

            # Drain stderr
            try:
                item = stderr_queue.get(timeout=1.0)
                if item is not None:
                    process_line(item, pbar, processed_targets, True)
            except queue.Empty:
                pass

        # Ensure readers have finished (they should put None when done)
        try:
            while stdout_queue.get(timeout=0.1) is not None:
                pass  # Drain any remaining non-None items (shouldn't happen)
        except queue.Empty:
            pass

        try:
            while stderr_queue.get(timeout=0.1) is not None:
                pass
        except queue.Empty:
            pass

        # Wait for threads to finish
        stdout_thread.join()
        stderr_thread.join()

        pbar.close()
        
        if process.returncode == 0:
            print("Feature extraction completed successfully")
            if total_targets:
                print(f"  Processed {len(processed_targets)}/{total_targets} targets")
            return True
        else:
            print(f"❌ Feature extraction failed with code {process.returncode}")
            return False
                
    except KeyboardInterrupt:
        print("\n⚠️  Feature extraction interrupted by user")
        if 'process' in locals():
            process.terminate()
            try:
                process.wait(timeout=5)
            except:
                process.kill()
        if 'pbar' in locals():
            pbar.close()
        return False
    except Exception as e:
        print(f"❌ Feature extraction failed: {e}")
        if 'process' in locals():
            process.terminate()
            try:
                process.wait(timeout=5)
            except:
                process.kill()
        if 'pbar' in locals():
            pbar.close()
        return False


def normalize_features(feature_file: str, ps_scaler_path: str, pb_scaler_path: str, 
                      output_dir: str) -> Tuple[str, str]:
    """
    Normalize features using the pre-trained scalers.
    
    Args:
        feature_file: Path to the raw feature CSV file
        ps_scaler_path: Path to SC2-PS scaler
        pb_scaler_path: Path to SC2-PB scaler
        output_dir: Directory to save normalized features
        
    Returns:
        Tuple of (ps_normalized_file, pb_normalized_file)
    """
    print("Normalizing features...")
    
    # Load raw features
    df = pd.read_csv(feature_file)
    df.fillna(0, inplace=True,axis=1)
    print(f"Loaded features: {df.shape[0]} compounds, {df.shape[1]} features")
    
    # Check if Id column exists, if not, it's likely the last column
    if 'Id' not in df.columns:
        # Assume the last column is the Id column
        df.columns = list(df.columns[:-1]) + ['Id']
        print("Detected Id column as the last column")
    
    # Create output directories
    ps_output_dir = os.path.join(output_dir, "sc2_ps")
    pb_output_dir = os.path.join(output_dir, "sc2_pb")
    os.makedirs(ps_output_dir, exist_ok=True)
    os.makedirs(pb_output_dir, exist_ok=True)
    
    # Extract base filename
    base_filename = os.path.basename(feature_file).replace('_features.csv', '_normalized.csv')
    
    # Load scalers
    try:
        ps_scaler = joblib.load(ps_scaler_path)
        pb_scaler = joblib.load(pb_scaler_path)
        print("Scalers loaded successfully")
    except Exception as e:
        raise RuntimeError(f"Failed to load scalers: {e}")
    
    # Separate IDs from features first
    ids = df['Id'].copy()
    X_all = df.drop(['Id'], axis=1)
    
    # Get expected columns from scalers
    try:
        ps_expected = list(ps_scaler.feature_names_in_)
        print(f"Aligned features to {len(ps_expected)} expected columns")
    except AttributeError:
        ps_expected = sorted(X_all.columns)  # Fallback: assume sorted order
        print("⚠️ PS scaler feature_names_in_ not available, using sorted columns (may cause issues)")

    try:
        pb_expected = list(pb_scaler.feature_names_in_)
        print(f"Aligned features to {len(pb_expected)} expected columns")
    except AttributeError:
        pb_expected = sorted(X_all.columns)
        print("⚠️ PB scaler feature_names_in_ not available, using sorted columns (may cause issues)")

    # Align features
    X_ps = X_all.reindex(columns=ps_expected, fill_value=0)
    X_pb = X_all.reindex(columns=pb_expected, fill_value=0)
    
    # Normalize features
    try:
        with tqdm(total=len(X_ps), desc="Normalizing PS features", unit="feature") as ps_pbar:
            X_ps_normalized = ps_scaler.transform(X_ps)
            ps_pbar.update(len(X_ps))

        with tqdm(total=len(X_pb), desc="Normalizing PB features", unit="feature") as pb_pbar:
            X_pb_normalized = pb_scaler.transform(X_pb)
            pb_pbar.update(len(X_pb))

        # Create normalized DataFrames
        df_ps_norm = pd.DataFrame(X_ps_normalized, columns=X_ps.columns)
        df_ps_norm.insert(0, 'Id', ids)
        
        df_pb_norm = pd.DataFrame(X_pb_normalized, columns=X_pb.columns)
        df_pb_norm.insert(0, 'Id', ids)
        
        # Save normalized features
        ps_output_file = os.path.join(ps_output_dir, base_filename)
        pb_output_file = os.path.join(pb_output_dir, base_filename)
        
        df_ps_norm.to_csv(ps_output_file, index=False)
        df_pb_norm.to_csv(pb_output_file, index=False)
        
        print(f"PS normalized features saved: {ps_output_file}")
        print(f"PB normalized features saved: {pb_output_file}")
        
        return ps_output_file, pb_output_file
        
    except Exception as e:
        raise RuntimeError(f"Failed to normalize features: {e}")


def load_models(ps_model_path: str, pb_model_path: str, use_gpu: bool = False) -> Tuple[xgb.Booster, xgb.Booster]:
    """
    Load the SC2-PS and SC2-PB XGBoost models.
    
    Args:
        ps_model_path: Path to SC2-PS model file
        pb_model_path: Path to SC2-PB model file  
        use_gpu: Whether to use GPU acceleration
        
    Returns:
        Tuple of (sc2_ps_model, sc2_pb_model)
    """
    try:
        # Load models
        sc2_ps = xgb.Booster()
        sc2_ps.load_model(ps_model_path)
        sc2_pb = xgb.Booster()
        sc2_pb.load_model(pb_model_path)
        
        # Set computation parameters
        params = {'tree_method': 'hist', 'device': 'cuda' if use_gpu else 'cpu'}
        sc2_ps.set_param(params)
        sc2_pb.set_param(params)
        
        print("Successfully loaded models:")
        print(f"  SC2-PS: {os.path.basename(ps_model_path)}")
        print(f"  SC2-PB: {os.path.basename(pb_model_path)}")
        print(f"  Computation device: {'GPU' if use_gpu else 'CPU'}")
        
        return sc2_ps, sc2_pb
        
    except Exception as e:
        raise RuntimeError(f"Failed to load models: {e}")


def load_normalized_features(ps_feature_path: str, pb_feature_path: str) -> Tuple[xgb.DMatrix, xgb.DMatrix, pd.Series]:
    """
    Load normalized feature data for both models.
    
    Args:
        ps_feature_path: Path to SC2-PS normalized feature CSV file
        pb_feature_path: Path to SC2-PB normalized feature CSV file
        
    Returns:
        Tuple of (ps_features, pb_features, ids)
    """
    try:
        # Load PS features
        df_ps = pd.read_csv(ps_feature_path)
        df_ps.fillna(0, inplace=True)
        
        # Load PB features  
        df_pb = pd.read_csv(pb_feature_path)
        df_pb.fillna(0, inplace=True)
        
        # Extract IDs (should be identical in both files)
        ids = df_ps['Id'].copy()
        
        # Prepare feature matrices
        X_ps = df_ps.drop(['Id'], axis=1, errors='ignore')
        X_pb = df_pb.drop(['Id'], axis=1, errors='ignore')
        
        # Convert to XGBoost DMatrix
        ps_features = xgb.DMatrix(X_ps, feature_names=X_ps.columns.tolist())
        pb_features = xgb.DMatrix(X_pb, feature_names=X_pb.columns.tolist())
        
        print("Loaded normalized features:")
        print(f"  PS features: {X_ps.shape[0]} compounds, {X_ps.shape[1]} features")
        print(f"  PB features: {X_pb.shape[0]} compounds, {X_pb.shape[1]} features")
        
        return ps_features, pb_features, ids
        
    except Exception as e:
        raise RuntimeError(f"Failed to load normalized feature data: {e}")


def get_base_compound_name(compound_id: str) -> str:
    """
    Extract base compound name by removing pose information.
    
    Args:
        compound_id: Full compound identifier
        
    Returns:
        Base compound name without pose suffix
    """
    # Split by '_pose' to get the base compound name
    if '_pose' in compound_id:
        return compound_id.split('_pose')[0]
    
    # Fallback: remove common pose suffixes like _out_pose1, etc.
    base_name = re.sub(r'_(out_)?pose\d+$', '', compound_id)
    base_name = re.sub(r'_\d+$', '', base_name)  # Remove trailing numbers
    return base_name


def aggregate_poses(results_df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate multiple poses by selecting the one with highest confidence.
    
    Args:
        results_df: DataFrame with compound results
        
    Returns:
        DataFrame with aggregated results and pose information
    """
    results_df['base_compound'] = results_df['compound_id'].apply(get_base_compound_name)
    
    # Find the pose with maximum confidence for each compound
    idx_max = results_df.groupby('base_compound')['sc2_score'].idxmax()
    aggregated_df = results_df.loc[idx_max].copy()
    
    # Add information about selected pose
    aggregated_df['selected_pose'] = aggregated_df['compound_id'].apply(
        lambda x: x.replace(get_base_compound_name(x), '').lstrip('_') or 'original'
    )
    
    # Count total poses per compound
    pose_counts = results_df.groupby('base_compound').size()
    aggregated_df['total_poses'] = aggregated_df['base_compound'].map(pose_counts)
    
    return aggregated_df


def perform_rescoring(ps_features: xgb.DMatrix, pb_features: xgb.DMatrix, ids: pd.Series,
                     sc2_ps: xgb.Booster, sc2_pb: xgb.Booster, 
                     ps_weight: float = 0.7, pb_weight: float = 0.3) -> pd.DataFrame:
    """
    Perform rescoring using the consensus SCORCH2 model.
    
    Args:
        ps_features: SC2-PS feature matrix
        pb_features: SC2-PB feature matrix
        ids: Compound identifiers
        sc2_ps: SC2-PS model
        sc2_pb: SC2-PB model
        ps_weight: Weight for PS model predictions
        pb_weight: Weight for PB model predictions
        
    Returns:
        DataFrame with rescoring results
    """
    print("Performing rescoring...")
    
    # Make predictions
    preds_ps = sc2_ps.predict(ps_features)
    preds_pb = sc2_pb.predict(pb_features)
    
    # Calculate consensus score
    consensus_scores = preds_ps * ps_weight + preds_pb * pb_weight
    
    # Create results dataframe (without weights in the main data)
    results_df = pd.DataFrame({
        'compound_id': ids,
        'sc2_ps_score': preds_ps,
        'sc2_pb_score': preds_pb, 
        'sc2_score': consensus_scores
    })
    
    # Sort by consensus score (descending)
    results_df = results_df.sort_values('sc2_score', ascending=False).reset_index(drop=True)
    results_df['rank'] = range(1, len(results_df) + 1)
    
    print(f"Rescoring completed for {len(results_df)} compounds")
    return results_df


def save_results(results_df: pd.DataFrame, output_path: str, should_aggregate_poses: bool = False,
                ps_model_path: str = "", pb_model_path: str = "", 
                ps_weight: float = 0.7, pb_weight: float = 0.3) -> None:
    """
    Save rescoring results to CSV file with detailed metadata.
    
    Args:
        results_df: DataFrame with results
        output_path: Path to save results
        should_aggregate_poses: Whether poses were aggregated
        ps_model_path: Path to PS model (for metadata)
        pb_model_path: Path to PB model (for metadata)
        ps_weight: PS model weight
        pb_weight: PB model weight
    """
    output_dir = os.path.dirname(output_path)
    if output_dir and output_dir != '':
        os.makedirs(output_dir, exist_ok=True)
    
    # Prepare final output
    if should_aggregate_poses:
        final_df = aggregate_poses(results_df)
        print(f"Aggregated {len(results_df)} poses into {len(final_df)} unique compounds")
    else:
        final_df = results_df.copy()
        final_df['selected_pose'] = 'N/A'
        final_df['total_poses'] = 1
    
    # Remove PS/PB weights from the output CSV
    output_columns = ['compound_id', 'sc2_ps_score', 'sc2_pb_score', 'sc2_score', 'rank']
    if should_aggregate_poses:
        output_columns.extend(['selected_pose', 'total_poses'])
    
    final_output_df = final_df[output_columns].copy()
    
    # Add metadata header
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    with open(output_path, 'w') as f:
        # Write metadata as comments
        f.write(f"# SCORCH2 Rescoring Results\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write(f"# SC2-PS Model: {os.path.basename(ps_model_path)}\n")
        f.write(f"# SC2-PB Model: {os.path.basename(pb_model_path)}\n")
        f.write(f"# Consensus Weights: PS={ps_weight:.2f}, PB={pb_weight:.2f}\n")
        f.write(f"# Poses Aggregated: {'Yes' if should_aggregate_poses else 'No'}\n")
        f.write(f"# Total Compounds: {len(final_output_df)}\n")
        f.write("#\n")
        
        # Write CSV data
        final_output_df.to_csv(f, index=False)
    
    print(f"Results saved to: {output_path}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"  Mean SC2 Score: {final_df['sc2_score'].mean():.4f}")
    print(f"  Std SC2 Score: {final_df['sc2_score'].std():.4f}")
    print(f"  Score Range: {final_df['sc2_score'].min():.4f} - {final_df['sc2_score'].max():.4f}")
    
    if should_aggregate_poses:
        multi_pose_compounds = final_df[final_df['total_poses'] > 1]
        if len(multi_pose_compounds) > 0:
            print(f"  Compounds with multiple poses: {len(multi_pose_compounds)}")
            print(f"  Average poses per multi-pose compound: {multi_pose_compounds['total_poses'].mean():.1f}")
    
    # Show top results (ensure they are sorted by score descending)
    print("\nTop 5 Compounds:")
    top_compounds = final_df.sort_values('sc2_score', ascending=False).head(5)
    for i, (_, row) in enumerate(top_compounds.iterrows(), 1):
        print(f"  {i}. {row['compound_id']} - Score: {row['sc2_score']:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="SCORCH2 Rescoring Code",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
Examples:
  # Complete workflow with demo data (raw results, no aggregation)
  python scorch2_rescoring.py --protein-dir example_data/protein --ligand-dir example_data/molecule \\
                              --sc2_ps_model /path/to/sc2_ps.xgb --sc2_pb_model /path/to/sc2_pb.xgb \\
                              --ps_scaler /path/to/sc2_ps_scaler --pb_scaler /path/to/sc2_pb_scaler \\
                              --output results.csv --gpu
  
  # Complete workflow with pose aggregation and keep all files
  python scorch2_rescoring.py --protein-dir example_data/protein --ligand-dir example_data/molecule \\
                              --sc2_ps_model /path/to/sc2_ps.xgb --sc2_pb_model /path/to/sc2_pb.xgb \\
                              --ps_scaler /path/to/sc2_ps_scaler --pb_scaler /path/to/sc2_pb_scaler \\
                              --output results.csv --aggregate --keep-temp --res-dir my_results --gpu
  
  # Skip feature extraction if features already exist
  python scorch2_rescoring.py --features existing_features.csv \\
                              --sc2_ps_model /path/to/sc2_ps.xgb --sc2_pb_model /path/to/sc2_pb.xgb \\
                              --ps_scaler /path/to/sc2_ps_scaler --pb_scaler /path/to/sc2_pb_scaler \\
                              --output results.csv
        """
    )
    
    # Model paths (required)
    parser.add_argument('--sc2_ps_model', type=str, required=True,
                       help="Path to the SC2-PS model file")
    parser.add_argument('--sc2_pb_model', type=str, required=True, 
                       help="Path to the SC2-PB model file")
    parser.add_argument('--ps_scaler', type=str, required=True,
                       help="Path to the SC2-PS scaler file")
    parser.add_argument('--pb_scaler', type=str, required=True,
                       help="Path to the SC2-PB scaler file")
    parser.add_argument('--output', type=str,
                       help="Output CSV file path to save rescoring results")
    
    # Input options (mutually exclusive)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('--protein-dir', type=str,
                            help="Directory containing protein PDBQT files (requires --ligand-dir)")
    input_group.add_argument('--features', type=str,
                            help="Path to existing feature CSV file (skip feature extraction)")
    
    # Additional required for protein-dir option
    parser.add_argument('--ligand-dir', type=str,
                       help="Directory containing ligand PDBQT files (required with --protein-dir)")
    
    # Processing options
    parser.add_argument('--num-cores', type=int, default=None,
                       help="Number of CPU cores for feature extraction")
    parser.add_argument('--aggregate', action='store_true',
                       help="Aggregate multiple poses by selecting the highest scoring pose (default: False, keep all poses)")
    parser.add_argument('--ps_weight', type=float, default=0.7,
                       help="Weight for SC2-PS predictions")
    parser.add_argument('--pb_weight', type=float, default=0.3,
                       help="Weight for SC2-PB predictions")
    parser.add_argument('--gpu', action='store_true',
                       help="Use GPU for prediction if available")
    
    # Output options
    parser.add_argument('--res-dir', type=str, default="res",
                       help="Results directory for intermediate and final files")
    parser.add_argument('--keep-temp', action='store_true',
                       help="Keep all intermediate files and save results to res-dir")
    
    args = parser.parse_args()
    
    # Validate arguments
    if args.protein_dir and not args.ligand_dir:
        parser.error("--ligand-dir is required when using --protein-dir")
    
    if abs(args.ps_weight + args.pb_weight - 1.0) > 1e-6:
        parser.error("PS and PB weights must sum to 1.0")
    
    print("=" * 60)
    print("SCORCH2 Rescoring Code")
    print("=" * 60)
    
    try:
        # Load models
        sc2_ps, sc2_pb = load_models(args.sc2_ps_model, args.sc2_pb_model, args.gpu)
        
        # Determine feature files
        if args.features:
            feature_files = [args.features]
            is_multi_target = False
        else:
            print(f"Processing structures from:")
            print(f"  Proteins: {args.protein_dir}")
            print(f"  Ligands: {args.ligand_dir}")
            
            temp_features_dir = os.path.join(args.res_dir, "features")
            
            # Create results directory
            os.makedirs(temp_features_dir, exist_ok=True)
            
            # Use provided num_cores or default
            num_cores = args.num_cores if args.num_cores is not None else max(1, os.cpu_count() - 1)
            
            # Run feature extraction
            if not run_feature_extraction(args.protein_dir, args.ligand_dir, 
                                  temp_features_dir, num_cores):
                raise RuntimeError("Feature extraction failed")
            
            # Collect feature files
            feature_files = [os.path.join(temp_features_dir, f) for f in os.listdir(temp_features_dir) 
                     if f.endswith('_protein_features.csv')]
            is_multi_target = True
            print(f"Found {len(feature_files)} target feature files")

        # Create normalized dir
        temp_normalized_dir = os.path.join(args.res_dir, "normalized_features")
        os.makedirs(temp_normalized_dir, exist_ok=True)

        all_results = [] if is_multi_target else None

        # Process each feature file
        for feature_file in tqdm(feature_files, desc="Processing targets" if is_multi_target else "Processing features"):
            # Extract target name
            base_name = os.path.basename(feature_file)
            target = base_name.replace('_protein_features.csv', '')
            
            # Normalize
            ps_normalized_file, pb_normalized_file = normalize_features(
                feature_file, args.ps_scaler, args.pb_scaler, temp_normalized_dir
            )
            
            # Load
            ps_features, pb_features, ids = load_normalized_features(
                ps_normalized_file, pb_normalized_file
            )
            
            # Rescore
            results_df = perform_rescoring(ps_features, pb_features, ids, sc2_ps, sc2_pb, 
                                   args.ps_weight, args.pb_weight)
            
            # Determine output path
            if is_multi_target:
                if args.keep_temp:
                    # Save target-wise results to res-dir when keeping files
                    output_dir = os.path.join(args.res_dir, 'target_results')
                else:
                    if args.output:
                        output_dir = os.path.join(os.path.dirname(args.output) or '.', 'results')
                    else:
                        output_dir = 'results'
                os.makedirs(output_dir, exist_ok=True)
                if args.output:
                    base_output = os.path.basename(args.output)
                    output_path = os.path.join(output_dir, f"{target}_{base_output}")
                else:
                    output_path = os.path.join(output_dir, f"{target}_results.csv")
            else:
                if args.keep_temp:
                    # Save single target result to res-dir when keeping files
                    output_dir = args.res_dir
                    os.makedirs(output_dir, exist_ok=True)
                    base_output = os.path.basename(args.output) if args.output else "results.csv"
                    output_path = os.path.join(output_dir, base_output)
                else:
                    output_path = args.output if args.output else "scorch2_results.csv"
            
            # Save target-wise results
            save_results(results_df, output_path, args.aggregate, 
                 args.sc2_ps_model, args.sc2_pb_model, args.ps_weight, args.pb_weight)
            if is_multi_target:
                results_df['target'] = target
                all_results.append(results_df.copy())
            print(f"Target-wise results saved to: {output_path}")

        if is_multi_target:
            aggregated_df = pd.concat(all_results, ignore_index=True)
            aggregated_df = aggregated_df.sort_values('sc2_score', ascending=False).reset_index(drop=True)
            aggregated_df['overall_rank'] = range(1, len(aggregated_df) + 1)
            
            # Save aggregated results
            if args.keep_temp:
                # Save aggregated results to res-dir when keeping files
                base_name = os.path.basename(args.output) if args.output else "aggregated_results.csv"
                aggregated_output_path = os.path.join(args.res_dir, f"aggregated_{base_name}")
            else:
                aggregated_output_path = args.output if args.output else "scorch2_aggregated_results.csv"
                
            save_results(aggregated_df, aggregated_output_path, args.aggregate, 
                 args.sc2_ps_model, args.sc2_pb_model, args.ps_weight, args.pb_weight)
            print(f"Aggregated results across {len(feature_files)} targets saved to: {aggregated_output_path}")
            
            # Also save to original output path if different and if output is specified
            if args.keep_temp and args.output and aggregated_output_path != args.output:
                save_results(aggregated_df, args.output, args.aggregate, 
                     args.sc2_ps_model, args.sc2_pb_model, args.ps_weight, args.pb_weight)
                print(f"Aggregated results also saved to original output path: {args.output}")

        # Clean up temporary files if requested
        if not args.keep_temp and os.path.exists(args.res_dir):
            import shutil
            shutil.rmtree(args.res_dir)
            print(f"Temporary files cleaned up")
        else:
            print(f"Keeping all files (features, normalized features, and results) in: {args.res_dir}")
        
        print("\nSCORCH2 rescoring completed successfully!")
        if args.keep_temp:
            print(f"All results and intermediate files saved to: {args.res_dir}")
            if is_multi_target:
                print(f"Target-wise results in: {os.path.join(args.res_dir, 'target_results')}")
                base_name = os.path.basename(args.output) if args.output else "aggregated_results.csv"
                print(f"Aggregated results: {os.path.join(args.res_dir, f'aggregated_{base_name}')}")
        if args.output:
            print(f"Final results saved to: {args.output}")
        else:
            print("Rescoring completed. No specific output file specified.")
        
    except Exception as e:
        print(f"\n❌ Error: {e}")
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
