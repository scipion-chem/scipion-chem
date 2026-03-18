import argparse
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
from itertools import combinations

def plotPcCoords(traj):
    """PCA of the coordinates"""
    pca = PCA(n_components=2, random_state=42)
    backboneIndices = topo.topology.select('backbone')
    traj.superpose(topo, atom_indices=backboneIndices)

    reducedCartesian = pca.fit_transform(
        traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)
    )

    time_ns = traj.time / 1000

    plt.figure(figsize=(8, 8))
    plt.scatter(reducedCartesian[:, 0], reducedCartesian[:, 1], marker="x", c=time_ns)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Cartesian coordinate PCA")
    cbar = plt.colorbar()
    cbar.set_label("Time [ns]")
    plt.show()


def plotPcaDist(traj):
    """PCA of the distances"""
    pca = PCA(n_components=2, random_state=42)
    caIndices = traj.topology.select('backbone')
    atom_pairs = list(combinations(caIndices, 2))
    pairwiseDistances = mdtraj.geometry.compute_distances(traj, atom_pairs)
    reducedDistances = pca.fit_transform(pairwiseDistances)

    timeNs = traj.time / 1000

    plt.figure(figsize=(8, 8))
    plt.scatter(reducedDistances[:, 0], reducedDistances[:, 1], marker="x", c=timeNs)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Pairwise distance PCA (backbone only)")
    cbar = plt.colorbar()
    cbar.set_label("Time [ns]")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MDTraj DSSP Analysis')
    parser.add_argument('-i', '--inputFilename', type=str, required=True)
    parser.add_argument('-t', '--inputTraj', type=str, required=True)

    parser.add_argument('--pca-coord', action='store_true')
    parser.add_argument('--pca-dist', action='store_true')

    args = parser.parse_args()

    topo = mdtraj.load(args.inputFilename)
    traj = mdtraj.load(args.inputTraj, top=topo)

    if args.pca_coord:
        plotPcCoords(traj)

    elif args.pca_dist:
        plotPcaDist(traj)
