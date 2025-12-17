#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generic PoseBusters test runner.

Reads a parameter file describing:
- which test to run
- molecule paths
- test-specific parameters

Loads and transforms molecules internally,
runs the selected PoseBusters test,
and writes results to a JSON file.
"""
import numpy as np
import os
import sys
import pprint
import pandas as pd

from rdkit import Chem


# -------------------------------------------------------------------------
# Utils
# -------------------------------------------------------------------------


def readParams(fname):
    """Read key=value parameter file."""
    params = {}
    with open(fname) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            if '=' not in line:
                continue
            key, val = line.split('=', 1)
            params[key.strip()] = val.strip()
    return params


def str2bool(val):
    return val.lower() in ('true', '1', 'yes', 'y')


def loadMolecule(molfile):
    """
    Load molecule and perform any required transformations.
    This is the ONLY place where molecules are handled.
    """
    if not os.path.isfile(molfile):
        raise FileNotFoundError(f"Molecule file not found: {molfile}")

    ext = os.path.splitext(molfile)[1].lower()

    if ext in ['.sdf', '.mol']:
        suppl = Chem.SDMolSupplier(molfile, removeHs=False)
        mol = suppl[0] if suppl and len(suppl) > 0 else None

    elif ext == '.mol2':
        mol = Chem.MolFromMol2File(molfile, removeHs=False)

    elif ext == '.pdb':
        mol = Chem.MolFromPDBFile(molfile, removeHs=False)

    else:
        raise ValueError(f"Unsupported molecule format: {ext}")

    if mol is None:
        raise RuntimeError(f"Failed to load molecule: {molfile}")

    return mol

def writeResults(results, outputDir):
    """Write PoseBusters results: summary in TXT, details tables as CSV."""
    import os
    import pandas as pd
    import pprint

    os.makedirs(outputDir, exist_ok=True)

    # --- summary ---
    summary = results.get('results', {})
    summary = cleanSummary(summary)

    summaryFile = os.path.join(outputDir, 'summary.txt')
    with open(summaryFile, 'w') as f:
        f.write(pprint.pformat(summary, indent=2))

    # --- details ---
    details = results.get('details', None)
    if details is not None:
        if isinstance(details, dict):
            for key, df in details.items():
                if isinstance(df, pd.DataFrame):
                    csvFile = os.path.join(outputDir, f"{key}.csv")
                    df.to_csv(csvFile, index=False)
                else:
                    with open(summaryFile, 'a') as f:
                        f.write(f"\n=== {key} ===\n")
                        f.write(pprint.pformat(df, indent=2))

        elif isinstance(details, pd.DataFrame):
            csvFile = os.path.join(outputDir, 'details.csv')
            details.to_csv(csvFile, index=False)
        else:
            with open(summaryFile, 'a') as f:
                f.write("\n=== details ===\n")
                f.write(pprint.pformat(details, indent=2))

    print(f"PoseBusters results written to folder: {os.path.basename(outputDir)}")


def cleanSummary(d):
    """Recursively convert numpy types to native Python types."""
    if isinstance(d, dict):
        return {k: cleanSummary(v) for k, v in d.items()}
    elif isinstance(d, list):
        return [cleanSummary(v) for v in d]
    elif isinstance(d, np.generic):
        return d.item()
    else:
        return d
# -------------------------------------------------------------------------
# PoseBusters tests
# -------------------------------------------------------------------------

def ruDdistanceGeometry(params):
    """Run PoseBusters distance geometry test."""
    mol = loadMolecule(params['mol_pred'])
    from posebusters.modules.distance_geometry import check_geometry

    results = check_geometry(
        mol_pred=mol,
        threshold_bad_bond_length=float(
            params.get('threshold_bad_bond_length', 0.2)
        ),
        threshold_clash=float(
            params.get('threshold_clash', 0.2)
        ),
        threshold_bad_angle=float(
            params.get('threshold_bad_angle', 0.2)
        ),
        ignore_hydrogens=str2bool(
            params.get('ignore_hydrogens', 'True')
        ),
        sanitize=str2bool(
            params.get('sanitize', 'True')
        ),
        symmetrize_conjugated_terminal_groups=str2bool(
            params.get('symmetrize_conjugated_terminal_groups', 'True')
        )
    )

    return results

def runEnergyRatio(params):
    """Run PoseBusters energy ratio test."""
    mol = loadMolecule(params['mol_pred'])
    from posebusters.modules.energy_ratio import check_energy_ratio

    results = check_energy_ratio(
        mol_pred=mol,
        threshold_energy_ratio=float(
            params.get('threshold_energy_ratio', 0.7)
        ),
        ensemble_number_conformations=float(
            params.get('ensemble_number_conformations', 100)
        )
    )
    return results

def runFlatness(params):
    """Run PoseBusters flatness test."""
    mol = loadMolecule(params['mol_pred'])
    from posebusters.modules.flatness import check_flatness

    results = check_flatness(
        mol_pred=mol,
        threshold_flatness=float(
            params.get('threshold_flatness', 0.1)
        ),
        check_nonflat=str2bool(
            params.get('check_nonflat', 'False')
        )
    )
    return results

def runIdentity(params):
    """Run PoseBusters flatness test."""
    mol = loadMolecule(params['mol_pred'])
    molTrue = loadMolecule(params['mol_true'])
    from posebusters.modules.identity import check_identity

    results = check_identity(
        mol_pred=mol,
        mol_true=molTrue,
    )
    return results

def runIntermolDistance(params):
    """Run PoseBusters intermolecular distance test."""
    mol = loadMolecule(params['mol_pred'])
    molCond = loadMolecule(params['mol_cond'])
    from posebusters.modules.intermolecular_distance import check_intermolecular_distance

    results = check_intermolecular_distance(
        mol_pred=mol,
        mol_cond=molCond,
        radius_type=str(
            params.get('radius_type', 'vdw')
        ),
        radius_scale=float(
            params.get('radius_scale', 0.8)
        ),
        clash_cutoff=float(
            params.get('clash_cutoff', 0.05)
        ),
        ignore_types = set(params.get('ignore_types', 'hydrogens').split(',')),
        max_distance=float(
            params.get('max_distance', 5.0)
        )
    )
    return results

def runRmsd(params):
    """Run PoseBusters intermolecular distance test."""
    mol = loadMolecule(params['mol_pred'])
    molTrue = loadMolecule(params['mol_true'])
    from posebusters.modules.rmsd import check_rmsd

    results = check_rmsd(
        mol_pred=mol,
        mol_true=molTrue,
        rmsd_threshold=float(
            params.get('rmsd_threshold', 2.0)
        ),
        heavy_only=str2bool(
            params.get('heavy_only', 'True')
        )
    )
    return results

def runVolOverlap(params):
    """Run PoseBusters intermolecular distance test."""
    mol = loadMolecule(params['mol_pred'])
    molCond = loadMolecule(params['mol_cond'])
    from posebusters.modules.volume_overlap import check_volume_overlap

    results = check_volume_overlap(
        mol_pred=mol,
        mol_cond=molCond,
        clash_cutoff=float(
            params.get('clash_cutoff', 0.05)
        ),
        vdw_scale=float(
            params.get('vdw_scale', 0.8)
        ),
        ignore_types = set(params.get('ignore_types', 'hydrogens').split(','))
    )
    return results



# -------------------------------------------------------------------------
# Dispatcher
# -------------------------------------------------------------------------

def runTest(params):
    test = params.get('test')

    if test == 'distance_geometry':
        return ruDdistanceGeometry(params)

    elif test == 'check_energy_ratio':
        return runEnergyRatio(params)
    elif test == 'check_flatness':
        return runFlatness(params)
    elif test == 'check_identity':
        return runIdentity(params)
    elif test == 'check_intermolecular_distance':
        return runIntermolDistance(params)
    elif test == 'check_rmsd':
        return runRmsd(params)
    elif test == 'check_volume_overlap':
        return runVolOverlap(params)

    else:
        raise NotImplementedError(f"Unknown or unsupported test: {test}")


# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

def main():
    paramsFile = sys.argv[1]
    if not os.path.isfile(paramsFile):
        raise FileNotFoundError(paramsFile)

    params = readParams(paramsFile)

    results = runTest(params)

    outDir = params['output']
    writeResults(results, outDir)



if __name__ == '__main__':
    main()