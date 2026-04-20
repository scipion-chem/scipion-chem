#!/usr/bin/env python
# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Extracted and adapted from TTClust by Thibault TUBIANA
# * Original source: https://github.com/tubiana/TTClust
# * License: GNU GPLv3
# *
# *
# **************************************************************************

import argparse
import datetime
import glob
import os
import sys
import matplotlib as mpl
import matplotlib.pyplot as plt
import mdtraj as md
import numpy as np
import operator
import progressbar as pg
import scipy.cluster.hierarchy as sch
from numba import jit
from hashlib import md5
from scipy.spatial.distance import is_valid_dm
from prettytable import PrettyTable
from sklearn import manifold

__version__ = '4.10.4'

# ==============================================================================
#                     GLOBAL VARIABLES
# ==============================================================================
WIDGETS = [pg.Bar('>'), ' ', pg.ETA(), ' ', pg.ReverseBar('<')]
COORDS = []
COLOR_LIST = ["red", "blue", "lime", "yellow",
              "darkorchid", "deepskyblue",
              "orange", "brown", "gray", "black",
              "darkgreen", "navy"]
DPI = 600


# ==============================================================================
#                          CLASS
# ==============================================================================

class Cluster:
    """
    Simple cluster class object (contains frames numbers, spread, size, ID
    and representative frame)
    """
    def __init__(self, num):
        self.frames = []
        self.spread = -1
        self.size = -1
        self.id = num
        self.representative = -1


# ==============================================================================
#                     TOOL FUNCTIONS
# ==============================================================================

def printScreenLogfile(string):
    """Print string on screen and write it on logfile."""
    print(string)
    LOGFILE.write("{}\n".format(string))
    LOGFILE.flush()


def writeCommandLine():
    """Write command line with quotes for -st, -sr and -sa arguments."""
    LOGFILE.write("command line       : python ")
    i = 0
    while i < len(sys.argv):
        LOGFILE.write("{} ".format(sys.argv[i]))
        if sys.argv[i] in ["-st", "-sr", "-sa"]:
            i += 1
            LOGFILE.write("\"{}\" ".format(sys.argv[i]))
        i += 1
    LOGFILE.write("\n")


def initLog(args, mdTrajectory):
    """Initialise the logfile with run information."""
    topoFile       = args["top"]
    trajFile       = args["traj"]
    selectionStr   = args["select_traj"]
    selectAlign    = args["select_alignement"]
    selectRmsd     = args["select_rmsd"]
    logName        = os.path.splitext(args["logfile"])[0]

    LOGFILE.write("========================================================\n")
    LOGFILE.write("====================  TTCLUST {}  ===================\n".format(__version__))
    LOGFILE.write("========================================================\n\n")
    LOGFILE.write("************ General information ************\n")
    LOGFILE.write("software version   : {}\n".format(__version__))
    LOGFILE.write("Created on         : {}\n".format(datetime.datetime.now()))
    writeCommandLine()
    LOGFILE.write("DESTINATION FOLDER : {}\n".format(os.getcwd() + "/" + logName))
    LOGFILE.write("ARGUMENTS : \n")
    LOGFILE.write("  Selection string :\n")
    LOGFILE.write("      Atoms selected in trajectory = {} \n".format(selectionStr))
    LOGFILE.write("      Atoms selected for alignement = {} \n".format(selectAlign))
    LOGFILE.write("      Atoms selected for RMSD = {} \n".format(selectRmsd))
    LOGFILE.write("  trajectory file  : {} \n".format(','.join(trajFile)))
    LOGFILE.write("   Stride          : {} \n".format(args["stride"]))
    LOGFILE.write("   Number of frames  : {} \n".format(mdTrajectory.n_frames))
    LOGFILE.write("   Number of atoms  : {} \n".format(mdTrajectory.n_atoms))
    LOGFILE.write("  topology file    : {} \n".format(topoFile))
    LOGFILE.write("  method used of clusterring : {}".format(args["method"]))
    LOGFILE.write("\n\n")
    if args["ngroup"]:
        LOGFILE.write("  Number of cluster asked: {}\n".format(args["ngroup"]))
    if args["cutoff"]:
        LOGFILE.write("  cutoff for dendrogram clustering: {}\n".format(args["cutoff"]))


def extractSelectedAtoms(selection, traj, logName, save=False):
    """
    Return a trajectory with only atoms given in arguments (through the
    selection string).
    """
    try:
        subTraj = traj.atom_slice(traj.top.select(selection))
        subTraj.center_coordinates()
        if save:
            printScreenLogfile("NOTE : 'st' argument given. I will save the subtrajectory"
                               " in {0}/{0}.xtc and topology file as {0}/{0}.pdb".format(logName))
            subTraj[0].save_pdb("{0}/{0}.pdb".format(logName))
            subTraj.save_xtc("{0}/{0}.xtc".format(logName))
        return subTraj
    except ValueError:
        print("ERROR : there is an error with your selection string")
        print("        SELECTION STRING : ")
        print("        {}".format(selection))
        print("        > Please check 'http://mdtraj.org/latest/atom_selection.html'")
        exit(1)


def sendErrorMessage(calcType, selectionString, other=""):
    """Print information regarding an invalid selection string."""
    print("ERROR : {} selection string not valid".format(calcType))
    print("        >{}".format(selectionString))
    if other != "":
        print("        >{}".format(other))
    exit(1)


def improveNucleicAcid(selectionString):
    """Improve RNA and DNA selection strings."""
    dna        = "(resname =~ '(5|3)?D([ATGC]){1}(3|5)?$')"
    rna        = "(resname =~ '(3|5)?R?([AUGC]){1}(3|5)?$')"
    backboneNa = "rna or dna and (name =~ \"(P)|(O[35]')|(C[3-5]')\")"
    base       = ("(rna or dna) and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") "
                  "and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H")
    baseRna    = ("(rna) and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") "
                  "and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H")
    baseDna    = ("dna and not (name =~ \"(P)|(O[35]')|(C[3-5]')\") "
                  "and not (name =~ \"(O[24]')|(O[123]P)|(C[12]')\") and not type H")

    if 'base_rna'   in selectionString: selectionString = selectionString.replace('base_rna',   baseRna)
    if 'base_dna'   in selectionString: selectionString = selectionString.replace('base_dna',   baseDna)
    if 'base'       in selectionString: selectionString = selectionString.replace('base',       base)
    if 'backbone_na' in selectionString: selectionString = selectionString.replace('backbone_na', backboneNa)
    if 'dna'        in selectionString: selectionString = selectionString.replace('dna',        dna)
    if 'rna'        in selectionString: selectionString = selectionString.replace('rna',        rna)
    return selectionString


def returnSelectionAtom(useFor, traj, args):
    """Return indices of selected atoms."""
    selectionString = args["select_alignement"]
    try:
        selection = traj.top.select(selectionString)
    except ValueError:
        sendErrorMessage(useFor, selectionString, "Keyword not recognized")

    if len(selection) == 0:
        if selectionString == "backbone":
            selection = traj.top.select(improveNucleicAcid("backbone_na"))
            args["select_alignement"] = "backbone_na"
            args["select_rmsd"]       = "backbone_na"
        if len(selection) == 0 or selection is None:
            sendErrorMessage(useFor, selectionString, "Selection list EMPTY")
        else:
            print("NOTE : Nucleic acids found.")
            print("       Automatic switch to nucleic acid mode")
    return selection


def saveDistMat(distMat, rmsdString, matrixType):
    """Save the numpy distance/linkage matrix for later reuse."""
    name = rmsdString.replace(" ", "_") if rmsdString else "matrix_all"
    try:
        np.save(name, distMat)
        printScreenLogfile("Saving {0} : {1}.npy".format(matrixType, name))
    except OSError:
        print("Selection line too long to be used as filename. TTClust will hash it")
        nameHashed = md5(name.encode()).hexdigest()
        print("{} -> {}".format(name, nameHashed))
        np.save(nameHashed, distMat)
        printScreenLogfile("Saving {0}: {1}.npy".format(matrixType, nameHashed))


def reorderCluster(clusters):
    """
    Reorder the cluster numbers so that the first frame belongs to cluster 1.
    """
    dictOrder = {}
    for i in range(len(clusters)):
        dictOrder[clusters[i].id] = clusters[i].frames[0]

    sortedClusters = sorted(dictOrder.items(), key=operator.itemgetter(1))
    for i in range(len(sortedClusters)):
        dictOrder[sortedClusters[i][0]] = i + 1

    for i in range(len(clusters)):
        clusters[i].id = dictOrder[clusters[i].id]


# ==============================================================================
#                     FUNCTIONS
# ==============================================================================

def parseArg():
    """Parse command-line arguments and return them as a dictionary."""
    arguments = argparse.ArgumentParser(
        description="Clusterize molecular dynamics trajectories "
                    "(Amber, Gromacs, CHARMM, NAMD, PDB)")

    arguments.add_argument('-f', "--traj",    required=True, nargs='+',
                           help="Trajectory file(s).")
    arguments.add_argument('-t', '--top',     default=None,
                           help="Topology file.")
    arguments.add_argument('-s', '--stride',  default=1, type=int,
                           help="Read every Nth frame.")
    arguments.add_argument('-l', '--logfile', default="clustering",
                           help="Logfile basename (default: clustering).")
    arguments.add_argument('-st', '--select_traj',      default="all",
                           help="Selection string for trajectory extraction.")
    arguments.add_argument('-sa', '--select_alignement', default="backbone",
                           help="Selection string for alignment. Use 'none' to skip.")
    arguments.add_argument('-sr', '--select_rmsd',      default="backbone",
                           help="Selection string for RMSD calculation.")
    arguments.add_argument('-m',  '--method', default="ward",
                           help="Linkage method: single/complete/average/weighted/centroid/median/ward.")
    arguments.add_argument('-rs', '--random_seed', type=int,
                           help="Seed for the random number generator used during the KMeans elbow method.")
    arguments.add_argument('-cc', "--cutoff", default=None,
                           help="Distance cutoff for hierarchical clustering.")
    arguments.add_argument('-ng', "--ngroup", default=None,
                           help="Number of clusters. Use 'auto' for elbow-based detection.")
    arguments.add_argument('-aa', "--autoclust", default="Y",
                           help="Autoclustering (Y/n).")
    arguments.add_argument('-i',  '--interactive', default="Y",
                           help="Reuse saved distance matrix if found (Y/n).")
    arguments.add_argument('-axis', '--axis', default="default",
                           help="Use 'frame' to display frame numbers instead of time.")
    arguments.add_argument('-limitmat', '--limitmat', default=100000000, type=int,
                           help="Frame limit for generating the distance matrix image.")

    args = vars(arguments.parse_args())

    args["autoclust"] = args["autoclust"] in ["Y", "y"]

    if args["autoclust"] and args["ngroup"] is None and args["cutoff"] is None:
        args["ngroup"] = "auto"

    return args


def askChoice(args, name):
    """
    Ask the user whether to reuse a found distance matrix file.
    Returns the filename to load, or None to recalculate.
    """
    if args["interactive"].upper() == "Y":
        print("Interactive mode disabled. I will use the saved matrix")
        return name

    print(" I found a distance matrix ({0}) saved. Do you want to use it ?".format(name))
    print("    y/Y - YES")
    print("    n/N - NO")
    print("    o/O - find all other .npy distance matrix")

    choice = input()

    if choice.upper() == "Y":
        printScreenLogfile(" >Distance matrix file detected : {0}".format(name))
        return name
    elif choice.upper() == "N":
        print("Calculation mode activated")
        return None
    elif choice.upper() == "O":
        npyFiles = glob.glob("*.npy")
        for i, f in enumerate(npyFiles):
            print("  {0} - {1}".format(i + 1, f))
        print(" -->Please choose and press Enter")
        choiceFile = input()
        try:
            name = npyFiles[int(choiceFile) - 1]
            return name
        except ValueError:
            print("I didn't understand. Please try again")
            return askChoice(args, name)
    else:
        print("I didn't understand. Please try again")
        return askChoice(args, name)


def searchDistMat(rmsdString, args):
    """Search whether a distance matrix .npy file already exists."""
    name = rmsdString.replace(" ", "_") if rmsdString else "matrix_all"

    nameHashed = md5(name.encode()).hexdigest() + '.npy'
    if not name.endswith(".npy"):
        name += ".npy"

    npyFiles = glob.glob("*.npy")

    if name in npyFiles and args["interactive"].lower() != "n":
        return askChoice(args, name)
    elif nameHashed in npyFiles and args["interactive"].lower() != "n":
        return askChoice(args, nameHashed)
    return None


def calculateRepresentativeFrameSpread(clustersList, distMat):
    """
    Choose the representative frame by calculating the mean RMSD of each
    structure within a cluster against all others.
    """
    print("Searching for representative frames")

    for cluster in clustersList:
        frames = cluster.frames
        meanRmsdPerFrame = {}

        for frameI in frames:
            meanRmsdPerFrame[frameI] = 0
            for frameJ in frames:
                if frameJ != frameI:
                    meanRmsdPerFrame[frameI] += distMat[frameI - 1, frameJ - 1]
            meanRmsdPerFrame[frameI] /= len(frames)

            repre = min(meanRmsdPerFrame, key=meanRmsdPerFrame.get)
            cluster.representative = repre
            cluster.spread = sum(meanRmsdPerFrame.values()) / len(frames) * 10


@jit(nopython=True, parallel=False, cache=True, nogil=True)
def calcRmsd2Frames(ref, frame):
    """RMSD calculation between a reference and a frame (JIT-compiled)."""
    dist = np.zeros(len(frame))
    for atom in range(len(frame)):
        dist[atom] = ((ref[atom][0] - frame[atom][0]) ** 2 +
                      (ref[atom][1] - frame[atom][1]) ** 2 +
                      (ref[atom][2] - frame[atom][2]) ** 2)
    return np.sqrt(dist.mean())


def createDM(traj, args):
    """Calculate or load the pairwise RMSD distance matrix."""
    if args["select_alignement"] != "none":
        alignmentSelection = returnSelectionAtom(useFor="ALIGNEMENT", traj=traj, args=args)
        trajAligned = traj.superpose(traj[0], atom_indices=alignmentSelection, parallel=True)
    else:
        trajAligned = traj

    untouchedRmsdString = args["select_rmsd"]
    rmsdString = improveNucleicAcid(args["select_rmsd"])

    if rmsdString:
        print("NOTE : Extraction of subtrajectory for time optimisation")
        trajAligned = extractSelectedAtoms(rmsdString, trajAligned, args["logname"])

    distances = np.zeros((traj.n_frames, traj.n_frames), dtype=np.float32)

    distanceFile = searchDistMat(untouchedRmsdString, args)

    if distanceFile:
        loaded = np.load(distanceFile)
        if loaded.shape[0] != traj.n_frames:
            printScreenLogfile(
                "WARNING: loaded matrix size ({}) != trajectory frames ({}). "
                "Recalculating.".format(loaded.shape[0], traj.n_frames))
        else:
            printScreenLogfile(" >Distance Matrix File Loaded!")
            return loaded

    pBar = pg.ProgressBar(widgets=WIDGETS, maxval=trajAligned.n_frames).start()
    counter = 0
    for i in range(trajAligned.n_frames):
        for j in range(i + 1, trajAligned.n_frames):
            rmsd = calcRmsd2Frames(trajAligned.xyz[i], trajAligned.xyz[j])
            distances[i][j] = rmsd
            distances[j][i] = rmsd
        pBar.update(counter)
        counter += 1
    pBar.finish()

    print("Calculation ended - saving distance matrix")
    saveDistMat(distances, args["select_rmsd"], "distance matrix")
    return distances


def onClick(event):
    """Get mouse coordinates on the matplotlib dendrogram window."""
    ix, iy = event.xdata, event.ydata
    global COORDS
    COORDS.append((ix, iy))
    if len(COORDS) == 1:
        plt.close(1)


def returnMappingCluster(labels):
    """Assign cluster objects to frames based on label list."""
    clustersList = []
    for clusterNum in set(labels):
        clustersList.append(Cluster(clusterNum))

    for i, clusterNum in enumerate(labels):
        clustersList[clusterNum - 1].frames.append(i)
        if clusterNum != clustersList[clusterNum - 1].id:
            print("{0} - {1}".format(clusterNum, clustersList[clusterNum - 1]))
            sys.exit(1)

    for cluster in clustersList:
        cluster.size = len(cluster.frames)
    return clustersList


def segmentsGain(p1, v, p2):
    """Compute angle gain between three points for elbow detection."""
    vp1  = np.linalg.norm(p1 - v)
    vp2  = np.linalg.norm(p2 - v)
    p1p2 = np.linalg.norm(p1 - p2)
    return np.arccos((vp1 ** 2 + vp2 ** 2 - p1p2 ** 2) / (2 * vp1 * vp2)) / np.pi


def autoClustering(matrix, randomSeed):
    """
    Determine the optimal number of clusters using KMeans and the elbow method.
    Returns the optimal k (int).
    """
    from sklearn.cluster import KMeans
    from scipy.spatial.distance import cdist

    distortions = []
    K = range(2, 15)
    for k in K:
        kMeans = KMeans(n_clusters=k, n_init=10, random_state=randomSeed)
        kMeans.fit(matrix)
        distortions.append(
            sum(np.min(cdist(matrix, kMeans.cluster_centers_, 'euclidean'), axis=1)) / matrix.shape[0])

    criterion = np.array(distortions)
    criterion = (criterion - criterion.min()) / (criterion.max() - criterion.min())

    segGains = np.array(
        [0] +
        [segmentsGain(*[np.array([K[j], criterion[j]]) for j in range(i - 1, i + 2)])
         for i in range(len(K) - 2)] +
        [np.nan])

    segThreshold = 0.99
    kIdx = np.argmax(segGains > segThreshold)

    kMeans = KMeans(n_clusters=kIdx, n_init=10, random_state=randomSeed)
    kMeans.fit(matrix)
    return kIdx


def createClusterTable(traj, args):
    """
    Build the distance matrix and run hierarchical clustering.
    Returns (distances, clusteringResult, linkage, cutoff).
    """
    print("         creating distance matrix")
    distances = createDM(traj, args)

    selectAlign        = args["select_alignement"]
    untouchedSelectRmsd = args["select_rmsd"]
    selectRmsd         = improveNucleicAcid(args["select_rmsd"])
    cutoff             = args["cutoff"]
    nCluster           = args["ngroup"] if args["ngroup"] in ("auto", None) else int(args["ngroup"])

    if selectRmsd is None:
        selectRmsd = "None"

    linkageKey  = selectAlign + '--' + untouchedSelectRmsd + " linkage " + args["method"]
    linkageFile = searchDistMat(linkageKey, args)

    if linkageFile:
        linkage = np.load(linkageFile)
    else:
        print("         Matrix shape: {}".format(distances.shape))
        print("         Scipy linkage in progress. Please wait. It can be long")
        linkageMethods = ['single', 'average', 'complete', 'weighted', 'centroid', 'median', 'ward']
        if args["method"] in linkageMethods:
            if not is_valid_dm(distances):
                print("THE DISTANCE MATRIX IS NOT VALID! Results may be unreliable.")
            else:
                print("The distance matrix is VALID.")
            linkage = sch.linkage(distances, method=args["method"])
        else:
            printScreenLogfile("ERROR : clustering method not recognized")
            printScreenLogfile("      : valid methods: single; complete; average; weighted; centroid; ward.")
            sys.exit(1)
        print("         >Done!")
        print("         ...Saving linkage matrix...")
        saveDistMat(linkage, linkageKey, "linkage matrix")
        print("         >Done!")

    if cutoff:
        cutoff = float(cutoff)
        clusteringResult = sch.fcluster(linkage, cutoff, "distance")
    elif nCluster:
        if nCluster == "auto":
            nCluster = autoClustering(distances, args['random_seed'])
        clusteringResult = sch.fcluster(linkage, t=nCluster, criterion="maxclust")
        nGroup = len(np.unique(clusteringResult))
        cutoff = linkage[-(nGroup - 1), 2]
    else:
        clicked = False
        while not clicked:
            fig = plt.figure()
            fig.canvas.mpl_connect('button_press_event', onClick)
            plt.title("Please click where you want to build clusters")
            sch.dendrogram(linkage)
            plt.show()
            try:
                cutoff = COORDS[0][1]
                clicked = True
                clusteringResult = sch.fcluster(linkage, cutoff, "distance")
            except ValueError:
                print("ERROR : PLEASE CLICK ON THE DENDROGRAM TO CHOOSE YOUR CUTOFF VALUE")
            plt.close()

    printScreenLogfile("  cutoff for clustering : {:.2f}".format(float(cutoff)))
    return distances, clusteringResult, linkage, cutoff


def writeRepresentativeFrame(traj, cluster, logName):
    """Write the representative frame of a cluster as a PDB file."""
    clusterNum = cluster.id
    frame      = cluster.representative
    size       = cluster.size
    traj[frame].save_pdb("{}/C{}-f{}-s{}.pdb".format(logName, clusterNum, frame + 1, size))


def getCmap(numCluster):
    """Return the matplotlib colormap to use for the given number of clusters."""
    global COLOR_LIST
    if numCluster > len(COLOR_LIST):
        cmap = "rainbow_r"
    else:
        COLOR_LIST = COLOR_LIST[:numCluster]
        cmap = mpl.colors.ListedColormap(COLOR_LIST)
    return cmap


def plotBarplot(clustersList, logName, size, traj, args):
    """Plot the linear frame-to-cluster assignment barplot."""
    clustersNumberOrdered = [0] * size
    for cluster in clustersList:
        for frame in cluster.frames:
            clustersNumberOrdered[frame] = cluster.id

    cmap = getCmap(len(clustersList))
    data = np.asmatrix(clustersNumberOrdered)

    _, _ = plt.subplots(figsize=(10, 1.5))

    if traj.time.sum() < 0.0000005 or args["axis"].lower() == "frame":
        timeMin, timeMax = 0, np.shape(data)[1]
        plt.xlabel("Frame")
    else:
        try:
            timeMin, timeMax = traj.time[0] / 1000, traj.time[-1] / 1000
            plt.xlabel("Time (ns)")
            if timeMax >= 1000:
                timeMin /= 1000
                timeMax /= 1000
                plt.xlabel("Time ($\\mu$s)")
        except ValueError:
            timeMin, timeMax = 0, np.shape(data)[1]

    im = plt.imshow(data, aspect='auto', interpolation='none', cmap=cmap,
                    extent=[timeMin, timeMax, 1, 0])
    plt.tight_layout()
    plt.tick_params(axis="y", which='both', left=False, right=False, labelleft=False)
    plt.tick_params(axis="x", direction="out", which='both', top=False)
    colorsList = im.cmap(im.norm(np.unique(clustersNumberOrdered)))

    plt.savefig("{0}/{1}-linear.png".format(logName, logName.split(os.sep)[-1]),
                dpi=DPI, transparent=False)
    plt.close()
    return colorsList


def plotHist(clustersList, logName, colorsList):
    """Plot a histogram of cluster sizes."""
    if mpl.__version__[0] == "2" and "classic" in plt.style.available:
        plt.style.use("classic")

    values = [cl.size for cl in clustersList]
    labels = [cl.id  for cl in clustersList]

    width = 0.7
    index = np.arange(len(values))
    _, ax = plt.subplots()

    bp = ax.bar(index, values, width, color=colorsList, label="Cluster size")
    for rect in bp:
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width() / 2., 1.0 * height,
                '%d' % int(height), ha='center', va='bottom')

    plt.xlabel("Clusters")
    plt.ylabel("Number of members")
    plt.title("Distribution within clusters")
    plt.xticks(index + (width / 2), labels)
    plt.tight_layout()
    plt.savefig("{0}/{1}-hist.png".format(logName, logName.split(os.sep)[-1]),
                dpi=DPI, transparent=False)
    plt.close()


def plotDistmat(distances, logName):
    """Plot the full pairwise RMSD distance matrix."""
    _, _ = plt.subplots()
    plt.imshow(distances, interpolation='none', origin='lower')
    plt.colorbar()
    plt.xlabel("Frame")
    plt.ylabel("Frame")
    plt.title("RMSD between frames (nm)")
    plt.savefig("{0}/{1}-distmat.png".format(logName, logName.split(os.sep)[-1]),
                dpi=DPI, transparent=False)
    plt.close()


def plotDendro(linkage, logName, cutoff, colorList, clustersList):
    """Plot the hierarchical clustering dendrogram."""
    if mpl.__version__[0] == "2" and "classic" in plt.style.available:
            plt.style.use("classic")

    _, ax = plt.subplots()
    colorHex   = [mpl.colors.rgb2hex(x) for x in colorList]
    sch.set_link_color_palette(colorHex)

    colorMember = {}
    for cl in clustersList:
        for frm in cl.frames:
            colorMember[frm] = mpl.colors.rgb2hex(colorList[cl.id - 1])

    linkCols = {}
    for i, i12 in enumerate(linkage[:, :2].astype(int)):
        c1, c2 = (linkCols[x] if x > len(linkage) else colorMember[x] for x in i12)
        linkCols[i + 1 + len(linkage)] = c1 if c1 == c2 else "#808080"

    sch.dendrogram(linkage, color_threshold=float(cutoff),
                   above_threshold_color="#808080",
                   link_color_func=lambda x: linkCols[x])

    plt.title("Clustering Dendrogram")
    ax.set_xticklabels([])
    plt.axhline(y=float(cutoff), color="grey")
    ax.set_ylabel("Distance (AU)")
    ax.set_xlabel("Frames")
    plt.savefig("{0}/{1}-den.png".format(logName, logName.split(os.sep)[-1]),
                format="png", dpi=DPI, transparent=False)
    plt.close()


def symmetrizeMatrix(matrix):
    """Return a symmetric version of a matrix (using the upper triangle)."""
    dim       = matrix.shape[0]
    matrixSym = np.copy(matrix)
    for i in range(dim):
        for j in range(i, dim):
            matrixSym[j, i] = matrix[i, j]
    return matrixSym


def plot2DDistanceProjection(rmsdMat, clustersList, colors, logName):
    """Create a 2D MDS projection of inter-cluster distances."""
    labels   = range(1, len(clustersList) + 1)
    rmsdNorm = rmsdMat / np.max(rmsdMat)
    rmsdNorm = symmetrizeMatrix(rmsdNorm)

    mds     = manifold.MDS(n_components=2, dissimilarity="precomputed", normalized_stress="auto")
    rmsdMds = mds.fit(rmsdNorm)
    coords  = rmsdMds.embedding_

    spreads = np.array([cl.spread for cl in clustersList])
    radii   = np.pi * (25 * spreads ** 2)
    x       = coords[:, 0]
    y       = coords[:, 1]

    plt.figure()
    ax  = plt.subplot(111)
    box = ax.get_position()
    ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

    ax.scatter(x, y, s=radii, c=colors, alpha=0.5)
    for label, xi, yi in zip(labels, x, y):
        plt.annotate(label, xy=(xi, yi), ha='left', va='bottom', fontsize=8)

    lims = list(ax.get_xlim()) + list(ax.get_ylim())
    ax.set_ylim((min(lims), max(lims)))
    ax.set_xlim((min(lims), max(lims)))

    plt.title("Relative distance between clusters")
    plt.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    plt.tick_params(axis='y', which='both', left="off", right="off", labelleft='off')

    maxSize  = max(radii)
    minSize  = min(radii)
    minColor = colors[np.argmin(radii)].copy()
    maxColor = colors[np.argmax(radii)].copy()
    minColor[-1] = 0.5
    maxColor[-1] = 0.5

    legMin = plt.scatter([], [], s=minSize, edgecolor='black', color=minColor)
    legMax = plt.scatter([], [], s=maxSize, edgecolor='black', color=maxColor)
    legend = ax.legend([legMin, legMax],
                       ["{:.2f}".format(min(spreads)), "{:.2f}".format(max(spreads))],
                       ncol=1, frameon=False, fontsize=8,
                       handlelength=2, loc="upper right", borderpad=1.8,
                       handletextpad=1, scatterpoints=1, bbox_to_anchor=(1.3, 0.9))
    legend.set_title('Spread radius', prop={"size": "small"})

    minRmsd      = np.min(rmsdMat[np.nonzero(rmsdMat)])
    maxRmsd      = np.max(rmsdMat[np.nonzero(rmsdMat)])
    textDistance = "RMSD\n   min : {:.2f}$ \\AA$\n   max : {:.2f} $\\AA$".format(minRmsd, maxRmsd)
    ax.annotate(textDistance, xy=(1.05, 0.5), xycoords="axes fraction", fontsize="small")

    plt.savefig("{0}/{1}-dist.png".format(logName, logName.split(os.sep)[-1]),
                format="png", dpi=DPI, transparent=False)
    plt.close()


def generateGraphs(clustersList, output, size, linkage, cutoff, distances, traj, args):
    """Generate all output graphs (barplot, dendrogram, histogram, distance matrix)."""
    colorsList = plotBarplot(clustersList, output, size, traj, args)
    plotDendro(linkage, output, cutoff, colorsList, clustersList)
    plotHist(clustersList, output, colorsList)

    if distances.shape[0] < args["limitmat"]:
        plotDistmat(distances, output)
    else:
        printScreenLogfile("Too many frames! The RMSD distance matrix image will not be generated.")
    return colorsList


def getRmsdCrossCluster(clustersList, distances, logName):
    """
    Compute and print the RMSD between representative frames of all clusters.
    Returns the RMSD matrix (np.array).
    """
    table      = PrettyTable()
    nClusters  = len(clustersList)
    fieldNames = ["Clusters"] + ["C" + str(x) for x in range(1, nClusters + 1)]
    table.field_names = fieldNames
    table.float_format = ".2"

    nonDiagValues = []
    rmsdMatrix    = np.zeros((nClusters, nClusters))

    for i in range(nClusters):
        repr1 = clustersList[i].representative
        for j in range(i + 1, nClusters):
            repr2 = clustersList[j].representative
            rmsd  = distances[repr1][repr2] * 10
            rmsdMatrix[i][j] = rmsdMatrix[j][i] = rmsd
            nonDiagValues.append(rmsd)

    for i in range(nClusters):
        table.add_row(["C" + str(i + 1)] + rmsdMatrix[i].tolist())

    printScreenLogfile("----------------------------")
    printScreenLogfile("RMSD MATRIX BETWEEN CLUSTERS")
    printScreenLogfile(table)

    if nonDiagValues:
        printScreenLogfile("\nAVERAGE RMSD BETWEEN CLUSTERS : {:.2f}".format(np.mean(nonDiagValues)))
    else:
        printScreenLogfile("\nAVERAGE RMSD BETWEEN CLUSTERS : N/A (only 1 cluster)")

    np.savetxt("{0}/RMSD_between_clusters.csv".format(logName), rmsdMatrix, delimiter=";")
    return rmsdMatrix


def clusterAnalysisCall(args):
    """
    Main pipeline: load trajectory, compute distance matrix, cluster, and
    generate all output files and graphs.
    """
    trajFile   = args["traj"]
    topFile    = args["top"]
    selectTraj = improveNucleicAcid(args["select_traj"])
    logName    = os.path.splitext(args["logfile"])[0]

    args["select_traj"]       = improveNucleicAcid(args["select_traj"])
    args["select_alignement"] = improveNucleicAcid(args["select_alignement"])

    print("NOTE : Per default the clustering is made on the BACKBONE of a PROTEIN")
    print("       PLEASE READ THE DOCUMENTATION AT https://www.github.com/tubiana/TTClust FOR PROPER USAGE \n")
    print("======= TRAJECTORY READING =======")

    if len(trajFile) == 1:
        trajFile = trajFile[0]
        if topFile is None and trajFile.endswith(".pdb"):
            traj = md.load_pdb(trajFile)
        else:
            traj = md.load(trajFile, top=topFile, stride=args["stride"])
    elif len(trajFile) > 1:
        print(">Several trajectories given. Will concatenate them.")
        trajList = []
        for t in trajFile:
            if topFile is None and t.endswith(".pdb"):
                trajList.append(md.load_pdb(t))
            else:
                trajList.append(md.load(t, top=topFile, stride=args["stride"]))
        traj = md.join(trajList)
        if traj.timestep > 0:
            traj.time = np.asarray(
                list(range(0, int(len(traj) * traj.timestep), int(traj.timestep))))
    else:
        print("ERROR: no trajectory given. TTClust will stop")
        sys.exit(1)

    initLog(args, traj)

    if selectTraj != "all":
        print("======= EXTRACTION OF SELECTED ATOMS =======")
        traj = extractSelectedAtoms(selectTraj, traj, args["logname"], save=True)

    print("====== Clustering ========")
    distances, clustersLabels, linkage, cutoff = createClusterTable(traj, args)

    printScreenLogfile("\n**** Cluster Results")
    clustersList = returnMappingCluster(clustersLabels)

    print("====== Reordering clusters ======")
    reorderCluster(clustersList)
    clustersList.sort(key=operator.attrgetter("id"))

    print("====== Generating Graph ======")
    colorsList = generateGraphs(clustersList, logName, traj.n_frames, linkage, cutoff, distances, traj, args)

    print("====== Calc. repr. frame  ======")
    calculateRepresentativeFrameSpread(clustersList, distances)

    for cluster in clustersList:
        printScreenLogfile("cluster {}".format(cluster.id))
        printScreenLogfile("    size = {}".format(cluster.size))
        printScreenLogfile("    representative frame={}".format(cluster.representative + 1))
        printScreenLogfile("    spread  : {0:.2f} ".format(cluster.spread))
        printScreenLogfile("    Frames : {} ".format(str([x + 1 for x in cluster.frames])))
        writeRepresentativeFrame(traj, cluster, logName)

    rmsdMatrix = getRmsdCrossCluster(clustersList, distances, logName)

    if len(clustersList) > 1:
        plot2DDistanceProjection(rmsdMatrix, clustersList, colorsList, logName)
    else:
        printScreenLogfile("WARNING: Only 1 cluster found — skipping 2D distance projection graph.")
        printScreenLogfile("         Consider lowering your cutoff value (currently {:.2f}).".format(float(cutoff)))

    return traj


def defineLogfile(log):
    """Define LOGFILE global when called from an external GUI."""
    global LOGFILE
    LOGFILE = log


def main():
    """Entry point: parse arguments, set up output folder, run clustering."""
    print("********************************************************")
    print("******************  TTCLUST {} *********************".format(__version__))
    print("********************************************************\n")

    args = parseArg()
    global LOGFILE

    if os.path.splitext(args["logfile"])[1] == "":
        args["logfile"] = args["logfile"] + ".log"

    logName        = os.path.splitext(args["logfile"])[0]
    args["logname"] = logName

    if not os.path.exists(logName):
        os.makedirs(logName)
    elif not os.path.isdir(logName):
        os.rename(logName, logName + ".bak")
        print("NOTE : A file with the same folder name was found and renamed "
              "to {}.bak".format(logName))
        os.makedirs(logName)

    fileName = args["logfile"].split(os.sep)[-1]
    LOGFILE  = open("{0}/{1}".format(logName, fileName), "w")
    traj     = clusterAnalysisCall(args)
    LOGFILE.close()
    return traj


# ==============================================================================
if __name__ == "__main__":
    main()