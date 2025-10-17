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


def runFeatureExtraction(proteinDir: str, ligandDir: str, outputDir: str, numCores: int = None) -> bool:
    """
    Run feature extraction using the SCORCH2 feature extraction utility.
    
    Args:
        proteinDir: Directory containing protein PDBQT files
        ligandDir: Directory containing ligand PDBQT files
        outputDir: Directory to save extracted features
        numCores: Number of CPU cores to use (default: os.cpu_count()-1)
        
    Returns:
        True if successful, False otherwise
    """
    if numCores is None:
        numCores = max(1, os.cpu_count() - 1)
    
    print("Running feature extraction...")
    print(f"  Using {numCores} CPU cores")
    
    # Create output directory
    os.makedirs(outputDir, exist_ok=True)
    
    # Get list of protein targets to estimate total progress
    try:
        proteinFiles = [f for f in os.listdir(proteinDir) if f.endswith('.pdbqt')]
        pdbids = list(set([f.split("_")[0] for f in proteinFiles]))
        totalTargets = len(pdbids)
        print(f"  Found {totalTargets} protein targets")
    except Exception as e:
        print(f"Warning: Could not count targets for progress tracking: {e}")
        totalTargets = None
    
    # Run feature extraction with progress monitoring
    cmd = [
        sys.executable, "utils/scorch2_feature_extraction.py",
        "--protein-dir", proteinDir,
        "--ligand-dir", ligandDir, 
        "--output-dir", outputDir,
        "--num-cores", str(numCores)
    ]
    
    try:
        # Initialize progress bar
        if totalTargets:
            pbar = tqdm(total=totalTargets, desc="Processing targets", unit="target")
        else:
            pbar = tqdm(desc="Feature extraction", unit="operation")
        
        print("Starting feature extraction subprocess...")
        
        # Initialize sets
        processedTargets = set()

        # Start the subprocess with unbuffered output
        env = {**os.environ, "PYTHONUNBUFFERED": "1"}
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, 
                                  text=True, universal_newlines=True, env=env)

        # Define line processing function
        def processLine(line: str, pbar, processedTargets, isStderr: bool = False):
            prefix = "⚠️ Subprocess: " if isStderr else ""
            line = line.strip()
            if not line:
                return

            # Saved features
            match = re.search(r"Saved features for (\w+)", line)
            if match:
                target = match.group(1)
                if target not in processedTargets:
                    processedTargets.add(target)
                    pbar.set_description(f"Completed: {target}")
                    pbar.update(1)
                return

            # Skipped features
            match = re.search(r"PDBID (\w+) Feature File exists", line)
            if match:
                target = match.group(1)
                if target not in processedTargets:
                    processedTargets.add(target)
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
        def readStdout(queue, pipe):
            for line in iter(pipe.readline, ''):
                queue.put(line)
            queue.put(None)

        def readStderr(queue, pipe):
            for line in iter(pipe.readline, ''):
                queue.put(line)
            queue.put(None)

        # Start reader threads
        stdoutQueue = queue.Queue()
        stderrQueue = queue.Queue()
        stdoutThread = threading.Thread(target=readStdout, args=(stdoutQueue, process.stdout), daemon=True)
        stderrThread = threading.Thread(target=readStderr, args=(stderrQueue, process.stderr), daemon=True)
        stdoutThread.start()
        stderrThread.start()

        # Process output in real-time
        # Monitor until process completes
        while process.poll() is None:
            updated = False

            # Process stdout
            try:
                item = stdoutQueue.get_nowait()
                if item is not None:
                    processLine(item, pbar, processedTargets, False)
                    updated = True
            except queue.Empty:
                pass

            # Process stderr
            try:
                item = stderrQueue.get_nowait()
                if item is not None:
                    processLine(item, pbar, processedTargets, True)
                    updated = True
            except queue.Empty:
                pass

            if not updated:
                time.sleep(0.1)

        # Process has completed - drain remaining output with timeout
        drainStart = time.time()
        while (time.time() - drainStart < 10) and (not stdoutQueue.empty() or not stderrQueue.empty()):
            # Drain stdout
            try:
                item = stdoutQueue.get(timeout=1.0)
                if item is not None:
                    processLine(item, pbar, processedTargets, False)
            except queue.Empty:
                pass

            # Drain stderr
            try:
                item = stderrQueue.get(timeout=1.0)
                if item is not None:
                    processLine(item, pbar, processedTargets, True)
            except queue.Empty:
                pass

        # Ensure readers have finished (they should put None when done)
        try:
            while stdoutQueue.get(timeout=0.1) is not None:
                pass  # Drain any remaining non-None items (shouldn't happen)
        except queue.Empty:
            pass

        try:
            while stderrQueue.get(timeout=0.1) is not None:
                pass
        except queue.Empty:
            pass

        # Wait for threads to finish
        stdoutThread.join()
        stderrThread.join()

        pbar.close()
        
        if process.returncode == 0:
            print("Feature extraction completed successfully")
            if totalTargets:
                print(f"  Processed {len(processedTargets)}/{totalTargets} targets")
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


def normalizeFeatures(featureFile: str, psScalerPath: str, pbScalerPath: str, 
                      outputDir: str) -> Tuple[str, str]:
    """
    Normalize features using the pre-trained scalers.
    
    Args:
        featureFile: Path to the raw feature CSV file
        psScalerPath: Path to SC2-PS scaler
        pbScalerPath: Path to SC2-PB scaler
        outputDir: Directory to save normalized features
        
    Returns:
        Tuple of (psNormalizedFile, pbNormalizedFile)
    """
    print("Normalizing features...")
    
    # Load raw features
    df = pd.read_csv(featureFile)
    df.fillna(0, inplace=True,axis=1)
    print(f"Loaded features: {df.shape[0]} compounds, {df.shape[1]} features")
    
    # Check if Id column exists, if not, it's likely the last column
    if 'Id' not in df.columns:
        # Assume the last column is the Id column
        df.columns = list(df.columns[:-1]) + ['Id']
        print("Detected Id column as the last column")
    
    # Create output directories
    psOutputDir = os.path.join(outputDir, "sc2_ps")
    pbOutputDir = os.path.join(outputDir, "sc2_pb")
    os.makedirs(psOutputDir, exist_ok=True)
    os.makedirs(pbOutputDir, exist_ok=True)
    
    # Extract base filename
    baseFilename = os.path.basename(featureFile).replace('_features.csv', '_normalized.csv')
    
    # Load scalers
    try:
        psScaler = joblib.load(psScalerPath)
        pbScaler = joblib.load(pbScalerPath)
        print("Scalers loaded successfully")
    except Exception as e:
        raise RuntimeError(f"Failed to load scalers: {e}")
    
    # Separate IDs from features first
    ids = df['Id'].copy()
    XAll = df.drop(['Id'], axis=1)
    
    # Get expected columns from scalers
    try:
        psExpected = list(psScaler.feature_names_in_)
        print(f"Aligned features to {len(psExpected)} expected columns")
    except AttributeError:
        psExpected = sorted(XAll.columns)  # Fallback: assume sorted order
        print("⚠️ PS scaler feature_names_in_ not available, using sorted columns (may cause issues)")

    try:
        pbExpected = list(pbScaler.feature_names_in_)
        print(f"Aligned features to {len(pbExpected)} expected columns")
    except AttributeError:
        pbExpected = sorted(XAll.columns)
        print("⚠️ PB scaler feature_names_in_ not available, using sorted columns (may cause issues)")

    # Align features
    XPs = XAll.reindex(columns=psExpected, fill_value=0)
    XPb = XAll.reindex(columns=pbExpected, fill_value=0)
    
    # Normalize features
    try:
        with tqdm(total=len(XPs), desc="Normalizing PS features", unit="feature") as psPbar:
            XPsNormalized = psScaler.transform(XPs)
            psPbar.update(len(XPs))

        with tqdm(total=len(XPb), desc="Normalizing PB features", unit="feature") as pbPbar:
            XPbNormalized = pbScaler.transform(XPb)
            pbPbar.update(len(XPb))

        # Create normalized DataFrames
        dfPsNorm = pd.DataFrame(XPsNormalized, columns=XPs.columns)
        dfPsNorm.insert(0, 'Id', ids)
        
        dfPbNorm = pd.DataFrame(XPbNormalized, columns=XPb.columns)
        dfPbNorm.insert(0, 'Id', ids)
        
        # Save normalized features
        psOutputFile = os.path.join(psOutputDir, baseFilename)
        pbOutputFile = os.path.join(pbOutputDir, baseFilename)
        
        dfPsNorm.to_csv(psOutputFile, index=False)
        dfPbNorm.to_csv(pbOutputFile, index=False)
        
        print(f"PS normalized features saved: {psOutputFile}")
        print(f"PB normalized features saved: {pbOutputFile}")
        
        return psOutputFile, pbOutputFile
        
    except Exception as e:
        raise RuntimeError(f"Failed to normalize features: {e}")


def loadModels(psModelPath: str, pbModelPath: str, useGpu: bool = False) -> Tuple[xgb.Booster, xgb.Booster]:
    """
    Load the SC2-PS and SC2-PB XGBoost models.
    
    Args:
        psModelPath: Path to SC2-PS model file
        pbModelPath: Path to SC2-PB model file  
        useGpu: Whether to use GPU acceleration
        
    Returns:
        Tuple of (sc2_ps_model, sc2_pb_model)
    """
    try:
        # Load models
        sc2Ps = xgb.Booster()
        sc2Ps.load_model(psModelPath)
        sc2Pb = xgb.Booster()
        sc2Pb.load_model(pbModelPath)
        
        # Set computation parameters
        params = {'tree_method': 'hist', 'device': 'cuda' if useGpu else 'cpu'}
        sc2Ps.set_param(params)
        sc2Pb.set_param(params)
        
        print("Successfully loaded models:")
        print(f"  SC2-PS: {os.path.basename(psModelPath)}")
        print(f"  SC2-PB: {os.path.basename(pbModelPath)}")
        print(f"  Computation device: {'GPU' if useGpu else 'CPU'}")
        
        return sc2Ps, sc2Pb
        
    except Exception as e:
        raise RuntimeError(f"Failed to load models: {e}")


def loadNormalizedFeatures(psFeaturePath: str, pbFeaturePath: str) -> Tuple[xgb.DMatrix, xgb.DMatrix, pd.Series]:
    """
    Load normalized feature data for both models.
    
    Args:
        psFeaturePath: Path to SC2-PS normalized feature CSV file
        pbFeaturePath: Path to SC2-PB normalized feature CSV file
        
    Returns:
        Tuple of (psFeatures, pbFeatures, ids)
    """
    try:
        # Load PS features
        dfPs = pd.read_csv(psFeaturePath)
        dfPs.fillna(0, inplace=True)
        
        # Load PB features  
        dfPb = pd.read_csv(pbFeaturePath)
        dfPb.fillna(0, inplace=True)
        
        # Extract IDs (should be identical in both files)
        ids = dfPs['Id'].copy()
        
        # Prepare feature matrices
        XPs = dfPs.drop(['Id'], axis=1, errors='ignore')
        XPb = dfPb.drop(['Id'], axis=1, errors='ignore')
        
        # Convert to XGBoost DMatrix
        psFeatures = xgb.DMatrix(XPs, feature_names=XPs.columns.tolist())
        pbFeatures = xgb.DMatrix(XPb, feature_names=XPb.columns.tolist())
        
        print("Loaded normalized features:")
        print(f"  PS features: {XPs.shape[0]} compounds, {XPs.shape[1]} features")
        print(f"  PB features: {XPb.shape[0]} compounds, {XPb.shape[1]} features")
        
        return psFeatures, pbFeatures, ids
        
    except Exception as e:
        raise RuntimeError(f"Failed to load normalized feature data: {e}")


def getBaseCompoundName(compoundId: str) -> str:
    """
    Extract base compound name by removing pose information.
    
    Args:
        compoundId: Full compound identifier
        
    Returns:
        Base compound name without pose suffix
    """
    # Split by '_pose' to get the base compound name
    if '_pose' in compoundId:
        return compoundId.split('_pose')[0]
    
    # Fallback: remove common pose suffixes like _out_pose1, etc.
    baseName = re.sub(r'_(out_)?pose\d+$', '', compoundId)
    baseName = re.sub(r'_\d+$', '', baseName)  # Remove trailing numbers
    return baseName


def aggregatePoses(resultsDf: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate multiple poses by selecting the one with highest confidence.
    
    Args:
        resultsDf: DataFrame with compound results
        
    Returns:
        DataFrame with aggregated results and pose information
    """
    resultsDf['base_compound'] = resultsDf['compound_id'].apply(getBaseCompoundName)
    
    # Find the pose with maximum confidence for each compound
    idxMax = resultsDf.groupby('base_compound')['sc2_score'].idxmax()
    aggregatedDf = resultsDf.loc[idxMax].copy()
    
    # Add information about selected pose
    aggregatedDf['selected_pose'] = aggregatedDf['compound_id'].apply(
        lambda x: x.replace(getBaseCompoundName(x), '').lstrip('_') or 'original'
    )
    
    # Count total poses per compound
    poseCounts = resultsDf.groupby('base_compound').size()
    aggregatedDf['total_poses'] = aggregatedDf['base_compound'].map(poseCounts)
    
    return aggregatedDf


def performRescoring(psFeatures: xgb.DMatrix, pbFeatures: xgb.DMatrix, ids: pd.Series,
                     sc2Ps: xgb.Booster, sc2Pb: xgb.Booster, 
                     psWeight: float = 0.7, pbWeight: float = 0.3) -> pd.DataFrame:
    """
    Perform rescoring using the consensus SCORCH2 model.
    
    Args:
        psFeatures: SC2-PS feature matrix
        pbFeatures: SC2-PB feature matrix
        ids: Compound identifiers
        sc2Ps: SC2-PS model
        sc2Pb: SC2-PB model
        psWeight: Weight for PS model predictions
        pbWeight: Weight for PB model predictions
        
    Returns:
        DataFrame with rescoring results
    """
    print("Performing rescoring...")
    
    # Make predictions
    predsPs = sc2Ps.predict(psFeatures)
    predsPb = sc2Pb.predict(pbFeatures)
    
    # Calculate consensus score
    consensusScores = predsPs * psWeight + predsPb * pbWeight
    
    # Create results dataframe (without weights in the main data)
    resultsDf = pd.DataFrame({
        'compound_id': ids,
        'sc2_ps_score': predsPs,
        'sc2_pb_score': predsPb, 
        'sc2_score': consensusScores
    })
    
    # Sort by consensus score (descending)
    resultsDf = resultsDf.sort_values('sc2_score', ascending=False).reset_index(drop=True)
    resultsDf['rank'] = range(1, len(resultsDf) + 1)
    
    print(f"Rescoring completed for {len(resultsDf)} compounds")
    return resultsDf


def saveResults(resultsDf: pd.DataFrame, outputPath: str, shouldAggregatePoses: bool = False,
                psModelPath: str = "", pbModelPath: str = "", 
                psWeight: float = 0.7, pbWeight: float = 0.3) -> None:
    """
    Save rescoring results to CSV file with detailed metadata.
    
    Args:
        resultsDf: DataFrame with results
        outputPath: Path to save results
        shouldAggregatePoses: Whether poses were aggregated
        psModelPath: Path to PS model (for metadata)
        pbModelPath: Path to PB model (for metadata)
        psWeight: PS model weight
        pbWeight: PB model weight
    """
    outputDir = os.path.dirname(outputPath)
    if outputDir and outputDir != '':
        os.makedirs(outputDir, exist_ok=True)
    
    # Prepare final output
    if shouldAggregatePoses:
        finalDf = aggregatePoses(resultsDf)
        print(f"Aggregated {len(resultsDf)} poses into {len(finalDf)} unique compounds")
    else:
        finalDf = resultsDf.copy()
        finalDf['selected_pose'] = 'N/A'
        finalDf['total_poses'] = 1
    
    # Remove PS/PB weights from the output CSV
    outputColumns = ['compound_id', 'sc2_ps_score', 'sc2_pb_score', 'sc2_score', 'rank']
    if shouldAggregatePoses:
        outputColumns.extend(['selected_pose', 'total_poses'])
    
    finalOutputDf = finalDf[outputColumns].copy()
    
    # Add metadata header
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    
    with open(outputPath, 'w') as f:
        # Write metadata as comments
        f.write(f"# SCORCH2 Rescoring Results\n")
        f.write(f"# Generated: {timestamp}\n")
        f.write(f"# SC2-PS Model: {os.path.basename(psModelPath)}\n")
        f.write(f"# SC2-PB Model: {os.path.basename(pbModelPath)}\n")
        f.write(f"# Consensus Weights: PS={psWeight:.2f}, PB={pbWeight:.2f}\n")
        f.write(f"# Poses Aggregated: {'Yes' if shouldAggregatePoses else 'No'}\n")
        f.write(f"# Total Compounds: {len(finalOutputDf)}\n")
        f.write("#\n")
        
        # Write CSV data
        finalOutputDf.to_csv(f, index=False)
    
    print(f"Results saved to: {outputPath}")
    
    # Print summary statistics
    print("\nSummary Statistics:")
    print(f"  Mean SC2 Score: {finalDf['sc2_score'].mean():.4f}")
    print(f"  Std SC2 Score: {finalDf['sc2_score'].std():.4f}")
    print(f"  Score Range: {finalDf['sc2_score'].min():.4f} - {finalDf['sc2_score'].max():.4f}")
    
    if shouldAggregatePoses:
        multiPoseCompounds = finalDf[finalDf['total_poses'] > 1]
        if len(multiPoseCompounds) > 0:
            print(f"  Compounds with multiple poses: {len(multiPoseCompounds)}")
            print(f"  Average poses per multi-pose compound: {multiPoseCompounds['total_poses'].mean():.1f}")
    
    # Show top results (ensure they are sorted by score descending)
    print("\nTop 5 Compounds:")
    topCompounds = finalDf.sort_values('sc2_score', ascending=False).head(5)
    for i, (_, row) in enumerate(topCompounds.iterrows(), 1):
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
    inputGroup = parser.add_mutually_exclusive_group(required=True)
    inputGroup.add_argument('--protein-dir', type=str,
                            help="Directory containing protein PDBQT files (requires --ligand-dir)")
    inputGroup.add_argument('--features', type=str,
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
        sc2Ps, sc2Pb = loadModels(args.sc2_ps_model, args.sc2_pb_model, args.gpu)
        
        # Determine feature files
        if args.features:
            featureFiles = [args.features]
            isMultiTarget = False
        else:
            print(f"Processing structures from:")
            print(f"  Proteins: {args.protein_dir}")
            print(f"  Ligands: {args.ligand_dir}")
            
            tempFeaturesDir = os.path.join(args.res_dir, "features")
            
            # Create results directory
            os.makedirs(tempFeaturesDir, exist_ok=True)
            
            # Use provided numCores or default
            numCores = args.num_cores if args.num_cores is not None else max(1, os.cpu_count() - 1)
            
            # Run feature extraction
            if not runFeatureExtraction(args.protein_dir, args.ligand_dir, 
                                  tempFeaturesDir, numCores):
                raise RuntimeError("Feature extraction failed")
            
            # Collect feature files
            featureFiles = [os.path.join(tempFeaturesDir, f) for f in os.listdir(tempFeaturesDir) 
                     if f.endswith('_protein_features.csv')]
            isMultiTarget = True
            print(f"Found {len(featureFiles)} target feature files")

        # Create normalized dir
        tempNormalizedDir = os.path.join(args.res_dir, "normalized_features")
        os.makedirs(tempNormalizedDir, exist_ok=True)

        allResults = [] if isMultiTarget else None

        # Process each feature file
        for featureFile in tqdm(featureFiles, desc="Processing targets" if isMultiTarget else "Processing features"):
            # Extract target name
            baseName = os.path.basename(featureFile)
            target = baseName.replace('_protein_features.csv', '')
            
            # Normalize
            psNormalizedFile, pbNormalizedFile = normalizeFeatures(
                featureFile, args.ps_scaler, args.pb_scaler, tempNormalizedDir
            )
            
            # Load
            psFeatures, pbFeatures, ids = loadNormalizedFeatures(
                psNormalizedFile, pbNormalizedFile
            )
            
            # Rescore
            resultsDf = performRescoring(psFeatures, pbFeatures, ids, sc2Ps, sc2Pb, 
                                   args.ps_weight, args.pb_weight)
            
            # Determine output path
            if isMultiTarget:
                if args.keep_temp:
                    # Save target-wise results to res-dir when keeping files
                    outputDir = os.path.join(args.res_dir, 'target_results')
                else:
                    if args.output:
                        outputDir = os.path.join(os.path.dirname(args.output) or '.', 'results')
                    else:
                        outputDir = 'results'
                os.makedirs(outputDir, exist_ok=True)
                if args.output:
                    baseOutput = os.path.basename(args.output)
                    outputPath = os.path.join(outputDir, f"{target}_{baseOutput}")
                else:
                    outputPath = os.path.join(outputDir, f"{target}_results.csv")
            else:
                if args.keep_temp:
                    # Save single target result to res-dir when keeping files
                    outputDir = args.res_dir
                    os.makedirs(outputDir, exist_ok=True)
                    baseOutput = os.path.basename(args.output) if args.output else "results.csv"
                    outputPath = os.path.join(outputDir, baseOutput)
                else:
                    outputPath = args.output if args.output else "scorch2_results.csv"
            
            # Save target-wise results
            saveResults(resultsDf, outputPath, args.aggregate, 
                 args.sc2_ps_model, args.sc2_pb_model, args.ps_weight, args.pb_weight)
            if isMultiTarget:
                resultsDf['target'] = target
                allResults.append(resultsDf.copy())
            print(f"Target-wise results saved to: {outputPath}")

        if isMultiTarget:
            aggregatedDf = pd.concat(allResults, ignore_index=True)
            aggregatedDf = aggregatedDf.sort_values('sc2_score', ascending=False).reset_index(drop=True)
            aggregatedDf['overall_rank'] = range(1, len(aggregatedDf) + 1)
            
            # Save aggregated results
            if args.keep_temp:
                # Save aggregated results to res-dir when keeping files
                baseName = os.path.basename(args.output) if args.output else "aggregated_results.csv"
                aggregatedOutputPath = os.path.join(args.res_dir, f"aggregated_{baseName}")
            else:
                aggregatedOutputPath = args.output if args.output else "scorch2_aggregated_results.csv"
                
            saveResults(aggregatedDf, aggregatedOutputPath, args.aggregate, 
                 args.sc2_ps_model, args.sc2_pb_model, args.ps_weight, args.pb_weight)
            print(f"Aggregated results across {len(featureFiles)} targets saved to: {aggregatedOutputPath}")
            
            # Also save to original output path if different and if output is specified
            if args.keep_temp and args.output and aggregatedOutputPath != args.output:
                saveResults(aggregatedDf, args.output, args.aggregate, 
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
            if isMultiTarget:
                print(f"Target-wise results in: {os.path.join(args.res_dir, 'target_results')}")
                baseName = os.path.basename(args.output) if args.output else "aggregated_results.csv"
                print(f"Aggregated results: {os.path.join(args.res_dir, f'aggregated_{baseName}')}")
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
