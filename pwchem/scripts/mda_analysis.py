# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Scipion-chem team
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

"""Protein-ligand and free-energy-landscape analyses of an MD trajectory, using
MDAnalysis. A single MDAnalysis ``Universe`` is built and reused across every
analysis (no duplicated trajectory loading). Plots are shown with
``plt.show()``, matching the other MDTraj analysis scripts; the Scipion viewer
launches this script with ``popen``, so the plot window does not block the GUI.

Analyses:
  * ``--distance``  Protein-ligand minimum distance per frame.
  * ``--hbonds``    Protein-ligand hydrogen-bond count per frame (both directions).
  * ``--fel``       Free energy landscape F(RMSD, Rg) = -RT ln(P).
"""

import argparse

import numpy as np
import matplotlib.pyplot as plt

import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
from MDAnalysis.analysis.rms import RMSD
from MDAnalysis.analysis.hydrogenbonds.hbond_analysis import HydrogenBondAnalysis

# Boltzmann constant in kcal/(mol K); F = -RT ln(P) is expressed in kcal/mol,
# matching the MD_quick_plot reference.
_KB_KCAL_MOL_K = 0.001987

# ------------------------------------------------------------------
# Trajectory loading (single Universe, reused everywhere)
# ------------------------------------------------------------------

def buildUniverse(topFile, trjFile):
    """Build the one and only MDAnalysis Universe for this run."""
    return mda.Universe(topFile, trjFile)


def resolveLigandSelection(universe, ligandSel=None, ligandID='LIG'):
    """Return a ligand selection string that actually matches atoms."""
    if ligandSel:
        return ligandSel if len(universe.select_atoms(ligandSel)) else None

    for resname in [ligandID, 'LIG', 'UNK', 'UNL', 'DRG']:
        if resname:
            sel = f'resname {resname}'
            if len(universe.select_atoms(sel)):
                return sel
    return None


def getTimesNs(universe):
    """Frame times in ns (falls back to frame index when the trajectory has no
    time information)."""
    times = np.array([ts.time for ts in universe.trajectory])
    if not np.any(times):
        times = np.arange(universe.trajectory.n_frames, dtype=float)
    return times / 1000.0

# ------------------------------------------------------------------
# Core analyses (reusable functions)
# ------------------------------------------------------------------

def compute_protein_ligand_distance(universe, protein_sel, ligand_sel):
    """Minimum protein-ligand distance per frame.

    Returns ``(time_ns, min_distance)`` as numpy arrays, with distances in Å.
    """
    protein = universe.select_atoms(protein_sel)
    ligand = universe.select_atoms(ligand_sel)
    if not len(protein) or not len(ligand):
        raise ValueError('Protein or ligand selection matched no atoms '
                         f'(protein="{protein_sel}": {len(protein)} atoms, '
                         f'ligand="{ligand_sel}": {len(ligand)} atoms).')

    times, minDist = [], []
    for ts in universe.trajectory:
        dmat = distance_array(protein.positions, ligand.positions, box=ts.dimensions)
        minDist.append(dmat.min())
        times.append(ts.time)
    times = np.array(times)
    if not np.any(times):
        times = np.arange(len(minDist), dtype=float)
    return times / 1000.0, np.array(minDist)


def compute_protein_ligand_hbonds(universe, protein_sel, ligand_sel,
                                   d_a_cutoff=3.5, d_h_a_angle_cutoff=150.0):
    """Protein-ligand hydrogen-bond count per frame."""
    if not len(universe.select_atoms(ligand_sel)):
        raise ValueError(f'Ligand selection "{ligand_sel}" matched no atoms.')

    nFrames = universe.trajectory.n_frames
    counts = np.zeros(nFrames, dtype=float)

    def _accumulate(donors_sel, hydrogens_sel, acceptors_sel):
        if not len(universe.select_atoms(hydrogens_sel)) or \
           not len(universe.select_atoms(acceptors_sel)):
            return
        hba = HydrogenBondAnalysis(universe=universe,
                                   donors_sel=donors_sel,
                                   hydrogens_sel=hydrogens_sel,
                                   acceptors_sel=acceptors_sel,
                                   d_a_cutoff=d_a_cutoff,
                                   d_h_a_angle_cutoff=d_h_a_angle_cutoff,
                                   update_selections=False)
        hba.run()
        counts[:] += hba.count_by_time()

    # protein donor -> ligand acceptor
    _accumulate(f'({protein_sel}) and (name N* O* S*)',
                f'({protein_sel}) and (name H*)',
                f'({ligand_sel}) and (name N* O* S* F*)')
    # ligand donor -> protein acceptor
    _accumulate(f'({ligand_sel}) and (name N* O* S*)',
                f'({ligand_sel}) and (name H*)',
                f'({protein_sel}) and (name N* O* S*)')

    return getTimesNs(universe), counts


def compute_rmsd_rg(universe, protein_sel='protein', align_sel='protein and name CA'):
    """RMSD (Å, aligned on ``align_sel``, first frame as reference) and radius of
    gyration (Å, on ``protein_sel``) per frame."""
    rmsdAnal = RMSD(universe, universe, select=align_sel, ref_frame=0)
    rmsdAnal.run()
    rmsd = rmsdAnal.results.rmsd[:, 2]

    protein = universe.select_atoms(protein_sel)
    rg = np.array([protein.radius_of_gyration() for _ in universe.trajectory])
    return rmsd, rg


def compute_free_energy_landscape(rmsd, rg, bins=50, temperature=300):
    """Free energy landscape F(RMSD, Rg) = -RT ln(P)."""
    hist, xedges, yedges = np.histogram2d(rmsd, rg, bins=bins)
    prob = hist / hist.sum()
    # Floor empty bins to a tiny probability so log() is finite everywhere.
    prob[prob == 0] = prob[prob > 0].min() * 0.01

    fel = -_KB_KCAL_MOL_K * temperature * np.log(prob)
    fel -= fel.min()  # normalize global minimum to 0
    return fel, xedges, yedges

# ------------------------------------------------------------------
# Plotting (displayed with plt.show, like the other MDTraj scripts)
# ------------------------------------------------------------------

def plot_distance(timeNs, minDist):
    plt.figure()
    plt.plot(timeNs, minDist, color='#1D3557', lw=1.5, label='Min. distance')
    plt.axhline(minDist.mean(), color='red', ls='--', alpha=0.7,
                label=f'Mean = {minDist.mean():.2f} Å')
    plt.axhline(4.0, color='green', ls=':', alpha=0.5, label='Contact threshold (4 Å)')
    plt.title('Protein-ligand minimum distance')
    plt.xlabel('Time (ns)')
    plt.ylabel('Minimum distance (Å)')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.show()

def plot_hbonds(timeNs, counts):
    import matplotlib.ticker as ticker
    plt.figure()
    plt.fill_between(timeNs, counts, color='#2DC653', alpha=0.35)
    plt.plot(timeNs, counts, color='#1A7A34', lw=1.2, label='H-bonds')
    plt.axhline(counts.mean(), color='red', ls='--', alpha=0.7,
                label=f'Mean = {counts.mean():.2f}')
    plt.gca().yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
    plt.title('Protein-ligand hydrogen bonds')
    plt.xlabel('Time (ns)')
    plt.ylabel('Number of H-bonds')
    plt.grid(alpha=0.3)
    plt.legend()
    plt.show()

def plot_fel(fel, xedges, yedges):
    # histogram2d indexes [x, y]; transpose so the mesh matches (X=RMSD, Y=Rg).
    xmesh, ymesh = np.meshgrid(xedges[:-1], yedges[:-1])

    plt.figure()
    cf = plt.contourf(xmesh, ymesh, fel.T, levels=25, cmap='viridis')
    cbar = plt.colorbar(cf)
    cbar.set_label('Free energy (kcal/mol)')

    # Mark the global free-energy minimum.
    minIdx = np.unravel_index(np.argmin(fel), fel.shape)
    plt.plot(xedges[minIdx[0]], yedges[minIdx[1]], 'r*', markersize=16,
             markeredgecolor='white', markeredgewidth=1, label='Global minimum')

    plt.title('Free energy landscape')
    plt.xlabel('RMSD (Å)')
    plt.ylabel('Radius of gyration (Å)')
    plt.legend()
    plt.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Protein-ligand and FEL analyses of an MD trajectory (MDAnalysis).')
    parser.add_argument('-i', '--inputFilename', type=str, required=True,
                        help='Input topology/system file (PDB, GRO, PRMTOP...).')
    parser.add_argument('-t', '--inputTraj', type=str, required=True,
                        help='Input trajectory file.')

    parser.add_argument('--distance', action='store_true',
                        help='Protein-ligand minimum distance per frame.')
    parser.add_argument('--hbonds', action='store_true',
                        help='Protein-ligand hydrogen-bond count per frame.')
    parser.add_argument('--fel', action='store_true',
                        help='Free energy landscape F(RMSD, Rg).')

    parser.add_argument('--protein-sel', dest='proteinSel', type=str, default='protein',
                        help='MDAnalysis selection for the protein.')
    parser.add_argument('--ligand-id', dest='ligandID', type=str, default='LIG',
                        help='Ligand residue name (from the MD system) used to build the '
                             'ligand selection automatically.')

    parser.add_argument('--bins', type=int, default=50, help='FEL histogram bins.')
    parser.add_argument('--temperature', type=float, default=300,
                        help='Temperature (K) used for the FEL.')

    args = parser.parse_args()

    # Single Universe, reused by every analysis below.
    universe = buildUniverse(args.inputFilename, args.inputTraj)

    if args.distance or args.hbonds:
        ligandSel = resolveLigandSelection(universe, ligandID=args.ligandID)
        if ligandSel is None:
            print(f'No ligand atoms found in the system (residue name "{args.ligandID}"); '
                  'skipping protein-ligand analysis.')
        else:
            if args.distance:
                try:
                    timeNs, minDist = compute_protein_ligand_distance(
                        universe, args.proteinSel, ligandSel)
                    plot_distance(timeNs, minDist)
                except Exception as e:
                    print(f'Protein-ligand distance analysis could not be run: {e}')
            if args.hbonds:
                try:
                    timeNs, counts = compute_protein_ligand_hbonds(
                        universe, args.proteinSel, ligandSel)
                    plot_hbonds(timeNs, counts)
                except Exception as e:
                    print(f'Protein-ligand hydrogen-bond analysis could not be run: {e}')

    if args.fel:
        try:
            rmsd, rg = compute_rmsd_rg(universe, protein_sel=args.proteinSel)
            fel, xedges, yedges = compute_free_energy_landscape(
                rmsd, rg, bins=args.bins, temperature=args.temperature)
            plot_fel(fel, xedges, yedges)
        except Exception as e:
            print(f'Free energy landscape analysis could not be run: {e}')
