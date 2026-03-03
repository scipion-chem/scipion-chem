import argparse
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from sklearn.decomposition import PCA
from itertools import combinations

def plot_pca_coords(traj):
    """PCA of the coordinates"""
    pca = PCA(n_components=2)
    backbone_indices = topo.topology.select('backbone')
    traj.superpose(topo, atom_indices=backbone_indices)

    reduced_cartesian = pca.fit_transform(
        traj.xyz.reshape(traj.n_frames, traj.n_atoms * 3)
    )

    plt.figure(figsize=(8, 8))
    plt.scatter(reduced_cartesian[:, 0], reduced_cartesian[:, 1], marker="x", c=traj.time)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Cartesian coordinate PCA")
    cbar = plt.colorbar()
    cbar.set_label("Time [ps]")
    plt.show()


def plot_pca_dist(traj):
    """PCA of the distances"""
    pca = PCA(n_components=2)
    ca_indices = traj.topology.select('backbone')
    atom_pairs = list(combinations(ca_indices, 2))
    pairwise_distances = mdtraj.geometry.compute_distances(traj, atom_pairs)
    reduced_distances = pca.fit_transform(pairwise_distances)

    plt.figure(figsize=(8, 8))
    plt.scatter(reduced_distances[:, 0], reduced_distances[:, 1], marker="x", c=traj.time)
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.title("Pairwise distance PCA (backbone only)")
    cbar = plt.colorbar()
    cbar.set_label("Time [ps]")
    plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MDTraj DSSP Analysis')
    parser.add_argument('-i', '--inputFilename', type=str, required=True)
    parser.add_argument('-t', '--inputTraj', type=str, required=True)
    parser.add_argument('-o', '--outputName', type=str, required=True)

    parser.add_argument('--pca-coord', action='store_true')
    parser.add_argument('--pca-dist', action='store_true')

    args = parser.parse_args()

    topo = mdtraj.load(args.inputFilename)
    traj = mdtraj.load(args.inputTraj, top=topo)

    if args.pca_coord:
        plot_pca_coords(traj)

    elif args.pca_dist:
        plot_pca_dist(traj)
