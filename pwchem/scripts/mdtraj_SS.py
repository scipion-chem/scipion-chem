import argparse
import mdtraj
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors


def plot_residue_dist(dssp, n_frames, n_res, outDir):
    """Panel 1: per-residue % SSE stacked bar"""
    helixPct = (dssp == 'H').sum(axis=0) / n_frames * 100
    strandPct = (dssp == 'E').sum(axis=0) / n_frames * 100

    plt.figure(figsize=(10, 5))
    plt.bar(range(n_res), helixPct, color='tomato', width=1.0, label='Helix')
    plt.bar(range(n_res), strandPct, color='steelblue', width=1.0, bottom=helixPct, label='Strand')

    plt.xlim(0, n_res)
    plt.ylim(0, 100)
    plt.xlabel('Residue Index')
    plt.ylabel('Res. % SSE')
    plt.title('Secondary Structure Distribution per Residue')
    plt.legend(loc='upper right')
    plt.show()


def plot_sse_timecourse(dssp, time, n_res, outDir):
    """Panel 2: % SSE per frame"""
    fHelix = (dssp == 'H').sum(axis=1) / n_res * 100
    fStrand = (dssp == 'E').sum(axis=1) / n_res * 100
    fTotal = fHelix + fStrand

    plt.figure(figsize=(14, 4))
    plt.fill_between(time, fTotal, color='gray', alpha=0.3, label='Total SSE')
    plt.plot(time, fHelix, color='tomato', lw=1.2, label='Helix')
    plt.plot(time, fStrand, color='steelblue', lw=1.2, label='Strand')

    plt.xlim(time[0], time[-1])
    plt.xlabel('Time (nsec)')
    plt.ylabel('% SSE')
    plt.title(f'Global SSE Content (Avg: {fTotal.mean():.2f}%)')
    plt.legend(loc='upper right')

    plt.show()


def plot_dssp_heatmap(dssp, time, n_res, outDir):
    """Panel 3: residue × time heatmap"""
    n_frames = len(time)
    dsspNum = np.zeros((n_frames, n_res), dtype=float)
    dsspNum[dssp == 'E'] = 1.0
    dsspNum[dssp == 'H'] = 2.0

    plt.figure(figsize=(14, 10))
    cmap = mcolors.ListedColormap(['#d0d0d0', 'steelblue', 'tomato'])

    im = plt.imshow(dsspNum.T, aspect='auto', origin='upper', cmap=cmap,
                    extent=[time[0], time[-1], n_res, 0],
                    vmin=0, vmax=2, interpolation='none')

    plt.xlabel('Time (nsec)')
    plt.ylabel('Residue Index')

    cbar = plt.colorbar(im, ticks=[0.33, 1.0, 1.67], pad=0.01)
    cbar.ax.set_yticklabels(['Loop', 'Strand', 'Helix'])

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='MDTraj DSSP Analysis')
    parser.add_argument('-i', '--inputFilename', type=str, required=True)
    parser.add_argument('-t', '--inputTraj', type=str, required=True)
    parser.add_argument('-o', '--outputName', type=str, required=True)

    parser.add_argument('--per-residue', action='store_true')
    parser.add_argument('--per-frame', action='store_true')
    parser.add_argument('--heatmap', action='store_true')

    args = parser.parse_args()

    # Load Trajectory
    topo = mdtraj.load(args.inputFilename)
    traj = mdtraj.load(args.inputTraj, top=topo)

    # Compute DSSP Once
    protIdx = traj.topology.select('protein')
    protTraj = traj.atom_slice(protIdx)
    dssp = mdtraj.compute_dssp(protTraj, simplified=True)

    n_frames, n_res = dssp.shape
    time = traj.time / 1000.0

    # 3. Conditional Plotting
    if args.per_residue:
        plot_residue_dist(dssp, n_frames, n_res, args.outputName)

    if args.per_frame:
        plot_sse_timecourse(dssp, time, n_res, args.outputName)

    if args.heatmap:
        plot_dssp_heatmap(dssp, time, n_res, args.outputName)