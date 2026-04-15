# **************************************************************************
# *
# * Authors: Joaquin Algorta (joaquin.algorta@cnb.csic.es)
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
"""
Clusterize molecular dynamics simulation trajectories to obtain representative
structures and visualize conformational evolution.
Based on TTClust: https://github.com/tubiana/TTClust
             DOI: https://pubs.acs.org/doi/10.1021/acs.jcim.8b00512
"""

import os
import glob

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects.data import AtomStruct, SetOfAtomStructs

from pwchem import Plugin as pwchemPlugin
from pwchem.objects import MDSystem
from pwchem.constants import MDTRAJ_DIC

# ── Clustering ────────────────────────────────────────────────────────────────
CLUSTERING_METHODS = ['ward', 'single', 'complete', 'average', 'weighted', 'centroid', 'median']
CLUSTER_MODE = ['Auto (elbow)', 'Number of groups', 'Cutoff']

# ── Atom selection ────────────────────────────────────────────────────────────
SEL_BACKBONE = 0
SEL_CA = 1
SEL_PROTEIN = 2
SEL_LIGAND = 3
SEL_RESIDUES = 4
SEL_NONE = 5  # alignment only — skips superposition

# Choices shown in the GUI
ALIGN_CHOICES = ['Backbone', 'Alpha carbons (CA)', 'Protein',
                 'Ligand', 'Residue range', 'None (skip alignment)']
RMSD_CHOICES = ['Backbone', 'Alpha carbons (CA)', 'Protein',
                'Ligand', 'Residue range']

# MDTraj syntax for fixed selections
_MDTRAJ_SEL = {
    SEL_BACKBONE: 'backbone',
    SEL_CA: 'name CA',
    SEL_PROTEIN: 'protein',
    # SEL_LIGAND is resolved dynamically from getLigandID()
    SEL_NONE: 'none',
}


def selectionToMdtraj(selIdx, firstRes=None, lastRes=None, ligandId=None):
    """Convert a GUI selection index to an MDTraj-compatible selection string."""
    if selIdx == SEL_RESIDUES:
        # MDTraj resid indices are 0-based
        return 'resid {} to {}'.format(int(firstRes), int(lastRes))
    if selIdx == SEL_LIGAND:
        if ligandId:
            return 'resname {}'.format(ligandId)
        # Fallback if the MDSystem has no ligand ID registered
        return 'not protein and not water'
    return _MDTRAJ_SEL[selIdx]


# ── Protocol ──────────────────────────────────────────────────────────────────

class ProtocolTrajectoryClustering(EMProtocol):
    """Clusterize molecular dynamics simulation trajectories using TTClust."""
    _label = 'Traj Clustering'

    # ── Form ──────────────────────────────────────────────────────────────────
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputMDSystem', params.PointerParam,
                      pointerClass='MDSystem', allowsNull=False,
                      label='MD simulation: ',
                      help='Molecular Dynamics System with an associated trajectory.')

        form.addParam('stride', params.IntParam, default=1,
                      label='Stride',
                      help='Read every Nth frame from the trajectory.')

        # ── Atom selection ────────────────────────────────────────────────────
        group = form.addGroup('Atom selection')

        group.addParam('selectAlign', params.EnumParam,
                       choices=ALIGN_CHOICES, default=SEL_BACKBONE,
                       label='Alignment selection',
                       help='Atoms used for structural superposition before RMSD '
                            'calculation. Choose "None (skip alignment)" if your '
                            'frames are already aligned or for internal-coordinate '
                            'based analyses.')
        line = group.addLine('Align residue range:',
                             condition='selectAlign == {}'.format(SEL_RESIDUES),
                             help='First and last residue index. The wizard shows'
                             'the file in which you can see the residue numbering.')
        line.addParam('firstResAlign', params.IntParam, default=0, label='First')
        line.addParam('lastResAlign', params.IntParam, default=0, label='Last')

        group.addParam('selectRmsd', params.EnumParam,
                       choices=RMSD_CHOICES, default=SEL_BACKBONE,
                       label='RMSD selection',
                       help='Atoms used for the pairwise RMSD distance matrix.\n'
                            'Tip: align on backbone + RMSD on ligand is the classic '
                            'setup for tracking ligand pose changes in a fixed '
                            'protein reference frame.')
        line = group.addLine('RMSD residue range:',
                             condition='selectRmsd == {}'.format(SEL_RESIDUES),
                             help='First and last residue index.'
                                  'The wizard shows the file in which you can see the residue numbering.')
        line.addParam('firstResRmsd', params.IntParam, default=0, label='First')
        line.addParam('lastResRmsd', params.IntParam, default=0, label='Last')

        # ── Clustering ────────────────────────────────────────────────────────
        group = form.addGroup('Clustering')

        group.addParam('clusterMethod', params.EnumParam,
                       choices=CLUSTERING_METHODS, default=0,
                       label='Linkage method',
                       help='Hierarchical clustering linkage method '
                            '(scipy.cluster.hierarchy.linkage).')

        group.addParam('clusterMode', params.EnumParam,
                       choices=CLUSTER_MODE, default=0,
                       label='Cluster determination',
                       help='Auto: KMeans elbow method selects k automatically.\n'
                            'Number of groups: you specify k.\n'
                            'Cutoff: you specify a distance threshold on the dendrogram.')

        group.addParam('nGroup', params.IntParam, default=5,
                       condition='clusterMode == 1',
                       label='Number of clusters',
                       help='Number of clusters (-ng).')

        group.addParam('cutoff', params.FloatParam, default=1.0,
                       condition='clusterMode == 2',
                       label='Distance cutoff',
                       help='Distance cutoff for dendrogram-based clustering (-cc).')

    # ── Steps ─────────────────────────────────────────────────────────────────
    def _insertAllSteps(self):
        self._insertFunctionStep(self.runClusteringStep)
        self._insertFunctionStep(self.createOutputStep)

    def runClusteringStep(self):
        mdsys = self.inputMDSystem.get()
        trajFile = os.path.abspath(mdsys.getTrajectoryFile())

        topPath = mdsys.getTopologyFile()
        topFile = os.path.abspath(
            mdsys.getSystemFile() if os.path.splitext(topPath)[1] == '.top'
            else topPath)

        os.makedirs(self._getExtraPath(), exist_ok=True)
        args = self._buildArgs(trajFile, topFile)
        pwchemPlugin.runScript(self, 'trajectory_clustering.py', args,
                               env=MDTRAJ_DIC, cwd=self._getExtraPath())

    def createOutputStep(self):
        """Collect TTClust representative PDB files into a SetOfAtomStructs."""
        logDir = self._getExtraPath('clustering')

        # TTClust names files: C{clusterID}-f{frame}-s{size}.pdb
        pdbFiles = sorted(
            glob.glob(os.path.join(logDir, 'C*.pdb')),
            key=lambda p: int(os.path.basename(p).split('-')[0][1:]))

        outputSet = SetOfAtomStructs().create(self._getPath())
        for pdbFile in pdbFiles:
            outputSet.append(AtomStruct(filename=os.path.abspath(pdbFile)))

        self._defineOutputs(outputAtomStructs=outputSet)
        self._defineSourceRelation(self.inputMDSystem, outputSet)

    # ── Helpers ───────────────────────────────────────────────────────────────
    def _getLigandId(self):
        """Safely retrieve the ligand residue name from the input MDSystem."""
        try:
            return self.inputMDSystem.get().getLigandID() or None
        except Exception:
            return None

    def _getSelectionStrings(self):
        """Return (alignStr, rmsdStr) as MDTraj selection strings."""
        ligandId = self._getLigandId()

        # int() cast guards against old saved protocol instances that stored
        # the string label (e.g. "backbone") instead of the integer index,
        # which would otherwise cause "Can't set backbone" warnings.
        alignIdx = int(self.selectAlign.get())
        rmsdIdx = int(self.selectRmsd.get())

        alignStr = selectionToMdtraj(alignIdx,
                                     self.firstResAlign.get(),
                                     self.lastResAlign.get(),
                                     ligandId=ligandId)
        rmsdStr = selectionToMdtraj(rmsdIdx,
                                    self.firstResRmsd.get(),
                                    self.lastResRmsd.get(),
                                    ligandId=ligandId)
        return alignStr, rmsdStr

    def _buildArgs(self, trajFile, topFile):
        alignStr, rmsdStr = self._getSelectionStrings()
        mode = int(self.clusterMode.get())

        args = '-f "{}" '.format(trajFile)
        args += '-t "{}" '.format(topFile)
        args += '-s {} '.format(self.stride.get())
        args += '-l clustering.log '
        args += '-sa "{}" '.format(alignStr)
        args += '-sr "{}" '.format(rmsdStr)
        args += '-m {} '.format(CLUSTERING_METHODS[int(self.clusterMethod.get())])
        if self._getLigandId():
            args += f'-st "protein or resname {self._getLigandId()}" '
        else:
            args += '-st "protein" '
        if mode == 0:  # Auto
            args += '-aa Y '
        elif mode == 1:  # N groups
            args += '-ng {} -aa n '.format(self.nGroup.get())
        elif mode == 2:  # Cutoff
            args += '-cc {} -aa n '.format(self.cutoff.get())

        return args

    # ── Info ──────────────────────────────────────────────────────────────────
    def _summary(self):
        summary = []
        mdsys = self.inputMDSystem.get()
        if mdsys:
            summary.append(
                'This protocol has created a set of atom structures as representatives of the clusters.\n '
                'To see the clustering results open the images inside extra/clustering path.\n')
            summary.append('Trajectory: {}'.format(mdsys.getTrajectoryFile()))
        try:
            alignStr, rmsdStr = self._getSelectionStrings()
            summary.append('Alignment selection : {}'.format(alignStr))
            summary.append('RMSD selection      : {}'.format(rmsdStr))
        except Exception:
            pass
        summary.append('Method : {}'.format(CLUSTERING_METHODS[int(self.clusterMethod.get())]))
        summary.append('Mode   : {}'.format(CLUSTER_MODE[int(self.clusterMode.get())]))
        return summary

    def _validate(self):
        errors = []
        mdsys = self.inputMDSystem.get()

        if mdsys and not mdsys.hasTrajectory():
            errors.append('The input MDSystem has no associated trajectory file.')

        if self.clusterMode.get() == 2 and self.cutoff.get() <= 0:
            errors.append('Cutoff must be a positive value.')

        if self.clusterMode.get() == 1 and self.nGroup.get() < 2:
            errors.append('Number of clusters must be at least 2.')

        return errors