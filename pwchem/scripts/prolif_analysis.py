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

if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputTopology <topFile> -t/--inputTraj [<trjFile>]
    -o/--outputName <outputName> 

    Performs several analysis on a MD trajectory using MDTraj
    '''
    parser = argparse.ArgumentParser(description='Run ProLIF Ligand-Target interaction analysis')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input system file (PDB or similar), reference')
    parser.add_argument('-o', '--outpuPath', type=str, help='Output name')
    parser.add_argument('-t', '--inputTraj', type=str, help='Input trajectory file')

    args = parser.parse_args()
    topoFile, outPath, trjFile = args.inputFilename, args.outpuPath, args.inputTraj

    u = mda.Universe(topoFile, trjFile)
    ligSelection = u.select_atoms("resname LIG")

    # protSelection = u.select_atoms("protein and byres around 20.0 group ligand", ligand=ligSelection)
    protSelection = plf.select_over_trajectory(u, u.trajectory[::10], "protein and byres around 6.0 resname LIG",
                                                   ligand=ligSelection)
    print(protSelection)

    fp = plf.Fingerprint()
    fp.run(u.trajectory[::10], ligSelection, protSelection)

    fp.to_pickle("fingerprint.pkl")

    fig = fp.plot_barcode()
    plt.savefig(os.path.join(outPath,'interaction_barcode.png'), dpi=300, bbox_inches='tight')

    ligMol = plf.Molecule.from_mda(ligSelection)
    view = fp.plot_lignetwork(ligMol)
    view.save(os.path.join(outPath,'interaction_network.html'))
    # view.save_png("") Solo funciona para notebooks dice

    frame = 0
    # seek specific frame
    u.trajectory[frame]
    ligand_mol = plf.Molecule.from_mda(ligSelection)
    protein_mol = plf.Molecule.from_mda(protSelection)
    # display
    view = fp.plot_3d(ligand_mol, protein_mol, frame=frame, display_all=False)
    view.save(os.path.join(outPath,'interaction_3d.html'))


