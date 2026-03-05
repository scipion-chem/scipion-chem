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
import MDAnalysis as mda
import prolif as plf
import matplotlib.pyplot as plt
import os

def save_results(fp, ligMol, outPath, prefix):
    df = fp.to_dataframe()
    # Check if the dataframe is empty
    if df.empty:
        print(f"!!! WARNING: No interactions were detected. Skipping plots.")
        # Create a dummy text file so Scipion knows the job 'finished'
        with open(os.path.join(outPath, 'no_interactions.txt'), 'w') as f:
            f.write("ProLIF analysis finished but found 0 interactions.")
        return
    fp.to_pickle(os.path.join(outPath,f"{prefix}_fingerprint.pkl"))

    df.to_csv(os.path.join(outPath, f'{prefix}_interactions.csv'), index=False)
    # Save Barcode
    fp.plot_barcode()
    plt.savefig(os.path.join(outPath, f'{prefix}_barcode.png'), dpi=300, bbox_inches='tight')
    plt.close()

    # Save LigNetwork
    view = fp.plot_lignetwork(ligMol)
    view.save(os.path.join(outPath, f'{prefix}_network.html'))

if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputTopology <topFile> -t/--inputTraj [<trjFile>]
    -o/--outputName <outputName> 

    Performs several analysis on a MD trajectory using MDTraj
    '''
    parser = argparse.ArgumentParser(description='Run ProLIF Ligand-Target interaction analysis')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input system file (PDB or similar), reference')
    parser.add_argument('-o', '--outpuPath', type=str, help='Output path')
    parser.add_argument('-t', '--inputTraj', type=str, help='Input trajectory file')
    parser.add_argument('-n', '--outputName', type=str, help='Output name')
    parser.add_argument('-wb', default=False, action='store_true', help='Run a water bridges analysis')

    args = parser.parse_args()
    topoFile, outPath, trjFile, outName = args.inputFilename, args.outpuPath, args.inputTraj, args.outputName

    u = mda.Universe(topoFile, trjFile)
    u.guess_TopologyAttrs(to_guess=['elements'])

    ligSelection = u.select_atoms("resname LIG")

    # protSelection = u.select_atoms("protein and byres around 20.0 group ligand", ligand=ligSelection)
    protSelection = plf.select_over_trajectory(u, u.trajectory[::10], "protein and byres around 6.0 group ligand",
                                                   ligand=ligSelection)

    if not args.wb:
        fp = plf.Fingerprint()
        fp.run(u.trajectory[::10], ligSelection, protSelection)

    else:
        water_names = "WAT HOH SPC TIP3 T3P TIP4"
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
        fp.run(u.trajectory[::10], ligSelection, protSelection)

    ligMol = plf.Molecule.from_mda(ligSelection)
    save_results(fp, ligMol, outPath, outName)


