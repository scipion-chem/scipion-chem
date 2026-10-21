# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Joaquin Algorta (ddelhoyo@cnb.csic.es)
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
# *  e-mail address 'joaquin.algorta@cnb.csic.es'
# *
# **************************************************************************

import argparse
import pickle
import MDAnalysis as mda
import prolif as plf
import matplotlib.pyplot as plt
import os
import pandas as pd
import seaborn as sns
from rdkit import DataStructs

def loadUniverse(topoFile, trjFile, gmxTop=None, includeDir=None):
    '''Universe of the system, with the bonds ProLIF needs to perceive the chemistry (H-bonds, aromaticity...):
    read from the GROMACS topology when given and readable, otherwise guessed from the coordinates'''
    trjKwargs = {'format': 'NCDF'} if os.path.splitext(trjFile)[1].lower() == '.netcdf' else {}
    if gmxTop:
        try:
            topKwargs = {'include_dir': includeDir} if includeDir else {}
            u = mda.Universe(gmxTop, trjFile, topology_format='ITP', infer_system=True, **topKwargs, **trjKwargs)
            print(f'Bonds read from the GROMACS topology {gmxTop}')
            return u
        except Exception as e:
            print(f'Could not read the GROMACS topology ({type(e).__name__}: {e})')

    u = mda.Universe(topoFile, trjFile, **trjKwargs)
    if not hasattr(u.atoms, 'bonds') or len(u.atoms.bonds) == 0:
        print(f'{topoFile} has no bonds: guessing them from the atom distances')
        u.guess_TopologyAttrs(to_guess=['elements', 'bonds'])
    return u

def saveResults(fp, ligMol, outPath, prefix):
    df = fp.to_dataframe()
    if df.empty:
        print("!!! WARNING: No interactions were detected.")
        # Create a dummy text file so Scipion knows the job 'finished'
        with open(os.path.join(outPath, 'no_interactions.txt'), 'w') as f:
            f.write("ProLIF analysis finished but found 0 interactions.")
        return

    # Save Fingerprint, and the ligand needed to redraw the interaction network with another threshold
    fp.to_pickle(os.path.join(outPath, f"{prefix}_fingerprint.pkl"))
    with open(os.path.join(outPath, f"{prefix}_ligand.pkl"), 'wb') as f:
        pickle.dump(ligMol, f)
    df.to_csv(os.path.join(outPath, f'{prefix}_interactions.csv'), index=False)

    # Save Barcode
    fp.plot_barcode()
    plt.title("Ligand-Target Interaction Barcode")
    plt.savefig(os.path.join(outPath, f'{prefix}_barcode.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Save LigNetwork
    view = fp.plot_lignetwork(ligMol)
    view.save(os.path.join(outPath, f'{prefix}_network.html'))

    # Save Tanimoto Similarity Matrix
    bitvectors = fp.to_bitvectors()
    similarityMatrixData = []
    for bv in bitvectors:
        similarityMatrixData.append(DataStructs.BulkTanimotoSimilarity(bv, bitvectors))
    similarityMatrix = pd.DataFrame(similarityMatrixData, index=df.index, columns=df.index)

    fig, ax = plt.subplots(figsize=(6, 5), dpi=200)
    colormap = sns.diverging_palette(300, 145, s=90, l=80, sep=30, center="dark", as_cmap=True)

    sns.heatmap(
        similarityMatrix,
        ax=ax,
        square=True,
        cmap=colormap,
        vmin=0,
        vmax=1,
        center=0.5,
        xticklabels=max(1, len(df) // 10),
        yticklabels=max(1, len(df) // 10),
    )
    ax.invert_yaxis()
    plt.yticks(rotation="horizontal")
    plt.title("Tanimoto Similarity Heatmap")
    fig.patch.set_facecolor("white")
    plt.savefig(os.path.join(outPath, f'{prefix}_matrix.png'), dpi=300, bbox_inches='tight')
    plt.close()


if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputTopology <topFile> -t/--inputTraj [<trjFile>]
    -o/--outputName <outputName> -f/--frameNum <frameNum>

    Generates ProLIF interaction fingerprint
    '''
    parser = argparse.ArgumentParser(description='Run ProLIF Ligand-Target interaction analysis')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input system file (PDB or similar), reference')
    parser.add_argument('-o', '--outpuPath', type=str, help='Output path')
    parser.add_argument('-t', '--inputTraj', type=str, help='Input trajectory file')
    parser.add_argument('-l', '--ligandName', type=str, help='ligand name in the topology')
    parser.add_argument('-n', '--outputName', type=str, help='Output name')
    parser.add_argument('-wb', default=False, action='store_true', help='Run a water bridges analysis')
    parser.add_argument('-f', '--frameNum', type=int, help='Frame number', default=10)
    parser.add_argument('--gmxTop', type=str, default=None,
                        help='GROMACS topology (.top) to read the bonds from (the -i file is used if it fails)')
    parser.add_argument('--includeDir', type=str, default=None,
                        help='Directory of the GROMACS force fields #included by the .top')

    args = parser.parse_args()
    topoFile, outPath, trjFile, ligName, outName, frameNum = (
        args.inputFilename, args.outpuPath, args.inputTraj, args.ligandName, args.outputName, args.frameNum)

    u = loadUniverse(topoFile, trjFile, args.gmxTop, args.includeDir)
    u.guess_TopologyAttrs(to_guess=['elements'])

    ligSelection = u.select_atoms(f"resname {ligName}")

    protSelection = plf.select_over_trajectory(u, u.trajectory[::frameNum], "protein and byres around 6.0 group ligand",
                                               ligand=ligSelection)

    if not args.wb:
        fp = plf.Fingerprint()
        fp.run(u.trajectory[::frameNum], ligSelection, protSelection)

    else:
        water_names = "WAT HOH SOL SPC TIP3 T3P TIP4"
        water_selection = u.select_atoms(
            f"resname {water_names} and byres around 8.0 (group ligand or group pocket)",
            ligand=ligSelection,
            pocket=protSelection,
            updating=True,
        )

        fp = plf.Fingerprint(
            ["WaterBridge"], parameters={"WaterBridge": {"water": water_selection, "order": 3}}
        )

        # for MD trajectories
        fp.run(u.trajectory[::frameNum], ligSelection, protSelection)

    ligMol = plf.Molecule.from_mda(ligSelection)
    saveResults(fp, ligMol, outPath, outName)
