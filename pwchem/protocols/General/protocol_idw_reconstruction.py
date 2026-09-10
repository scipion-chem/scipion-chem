# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Carlota Laverón (carlota.laveronvilas@usp.ceu.es)
# *
# * CEU San Pablo University
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
# *  e-mail address 'carlota.laveronvilas@usp.ceu.es'
# *
# **************************************************************************

import os
import copy

import numpy as np
from scipy.spatial import cKDTree

from pyworkflow.constants import BETA
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol
from pwem.objects import SetOfAtomStructs, AtomStruct
from Bio.PDB import PDBIO, PDBParser, MMCIFParser
from pwchem.objects.idw import InvDistTree3D

DEFAULT_SEARCH_RADIUS = 15.0
DEFAULT_IDW_POWER = 2.0
OUTPUT_STRUCTUresName = 'reconstructedStructures'


class ProtocolInverseDistanceWeighting(EMProtocol):
    """
    Reconstructs full-atom protein structures from a C_alpha-only conformational
    ensemble using compact-support Inverse Distance Weighting (IDW) interpolation
    with KD-tree neighbour search.

    The protocol transfers the conformational changes observed in a C_alpha-only
    ensemble onto a full-atom reference structure. For each target conformation,
    C_alpha atoms in the reference and target structures are matched by residue
    identifiers after automatically pairing corresponding chains. The resulting
    C_alpha displacements are then interpolated to all atoms of the reference
    structure using compact-support IDW.

    Workflow
    --------
    1. Receive a full-atom reference structure and a SetOfAtomStructs containing
       C_alpha-only target conformations.
    2. Extract the C_alpha coordinates from the reference structure and use them
       as the source control points for the interpolation.
    3. Determine the number of IDW neighbours either automatically from the
       reference structure or from the user-defined maximum.
    4. For each target conformation:
       - Extract its C_alpha coordinates.
       - Match reference and target chains based on the number of common residue
         identifiers.
       - Match C_alpha atoms using residue number and insertion code.
       - Calculate the displacement of each matched C_alpha atom.
       - Interpolate these displacements to all atoms of the reference structure
         using compact-support IDW.
       - Apply the interpolated displacements to generate a reconstructed
         full-atom structure.
    5. Save each reconstructed structure as a PDB file.
    6. Collect all reconstructed structures into a SetOfAtomStructs.

    Input
    -----
    - inputReference:
        Full-atom reference structure in PDB or CIF format.

        The C_alpha coordinates of this structure are used as the source
        control points, while the complete atomic coordinates are used as the
        reference geometry to be reconstructed.

    - inputEnsemble:
        SetOfAtomStructs containing C_alpha-only protein structures representing
        different conformations of the same protein.

        Each structure provides the target C_alpha coordinates whose
        conformational displacements are transferred to the full-atom reference.

    Parameters
    ----------
    - Search radius R:
        Maximum distance in Angstroms within which reference C_alpha atoms
        contribute to the interpolation of an atom.

        C_alpha atoms located beyond this radius have zero contribution.
        The default value is 15.0 A.

    - Maximum neighbours k:
        Maximum number of nearest C_alpha atoms considered for the IDW
        interpolation.

        When set to 0, k is determined automatically from the reference
        structure based on the mean number of C_alpha atoms within the selected
        search radius.

    - IDW power p:
        Controls how strongly the interpolation favours nearby C_alpha atoms.

        Higher values give greater weight to nearby control points and therefore
        make the deformation more locally influenced. The default value is 2.0.

    Interpolation Method
    --------------------
    The protocol uses compact-support Inverse Distance Weighting (IDW) to
    interpolate the C_alpha displacements throughout the full-atom structure.

    For an atom at position x, only C_alpha control points within the search
    radius R are considered. Their contribution is weighted according to their
    distance from x and the selected IDW power p.

    A KD-tree is used to efficiently identify nearby C_alpha neighbours during
    the interpolation.

    Chain and Residue Matching
    --------------------------
    Reference and target C_alpha atoms are matched independently for each
    target structure.

    Chain IDs do not need to be identical between the reference and target
    structures. Instead, chains are paired according to the number of common
    residue identifiers. Once chains are paired, residues are matched using
    residue number and insertion code.

    Chains without a corresponding match are reported as warnings and are
    excluded from the displacement interpolation.

    Output
    ------
    - reconstructedStructures:
        SetOfAtomStructs containing the reconstructed full-atom structures,
        one for each successfully processed target conformation.

        Each output structure retains the topology and atoms of the full-atom
        reference structure while its atomic coordinates are displaced according
        to the conformational changes observed in the corresponding C_alpha-only
        target structure.

    Summary
    -------
    The protocol reports:

    - Number of reconstructed full-atom structures
    - Search radius R
    - Number of IDW neighbours k, either automatically determined or explicitly
      selected
    - IDW power p

    Use Cases
    ---------
    - Reconstruction of full-atom structures from coarse-grained C_alpha-only
      conformational ensembles
    - Conversion of C_alpha-based conformational sampling into full-atom models
    - Generation of full-atom ensembles from normal-mode or other C_alpha-based
      sampling methods
    - Preparation of conformational ensembles for downstream structural analysis
      or molecular modelling

    Notes
    -----
    The protocol does not perform conventional molecular-mechanics refinement
    after reconstruction. The output coordinates are obtained by transferring
    and interpolating the C_alpha displacements onto the atoms of the reference
    structure.

    The quality of the reconstruction depends on the similarity between the
    reference and target structures, the availability of matching C_alpha atoms,
    and the selected IDW parameters.

    The search radius R and IDW power p control the spatial extent and locality
    of the propagated deformation, respectively.
    """

    _label = 'IDW Backbone Reconstruction'
    _devStatus = BETA

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam(
            'inputReference',
            params.PointerParam,
            pointerClass='AtomStruct', allowsNull=False,
            label='Reference full-atom structure: ',
            help='PDB/CIF file containing ALL atoms. Its C_alfa positions '
                 'serve as source control points for the IDW interpolation.',
        )

        form.addParam(
            'inputEnsemble',
            params.PointerParam,
            pointerClass='SetOfAtomStructs', allowsNull=False,
            label='Target C_alfa ensemble: ',
            help='Set of C_alfa-only models (one per conformer), e.g. from '
                 'prody2 - ANM MC walks. Each model provides destination '
                 'C_alfa positions.',
        )

        form.addParam(
            'searchRadius',
            params.FloatParam,
            default=DEFAULT_SEARCH_RADIUS,
            label='Search radius R (A): ',
            help='Maximum distance in Angstroms to search for C_alfa neighbours. '
                 'Ca beyond this radius contribute exactly 0. '
                 'Typical values: 10-20 A.\n\n'
                 'Formula: w_k = ((R - d) / (R * d)) ^ p',
        )

        form.addParam(
            'maxNeighbours',
            params.IntParam,
            default=0,
            label='Maximum neighbours k: ',
            help='Maximum number of nearest C_alfa atoms used for the IDW interpolation. '
                 'Set to 0 to determine k automatically from the reference structure. '
                 'If fewer C_alfa atoms are found within the search radius, all available '
                 'neighbours are used.',
        )

        form.addParam(
            'idwPower',
            params.FloatParam,
            default=DEFAULT_IDW_POWER,
            label='IDW power p: ',
            help='Power controlling the influence of distance in the IDW interpolation. '
                 'Higher values give more weight to nearby C_alfa atoms. '
                 'Default: 2.0.',
        )

    def _insertAllSteps(self):
        self._insertFunctionStep('reconstructStep')
        self._insertFunctionStep('createOutputStep')

    def reconstructStep(self):

        R = self.searchRadius.get()
        p = self.idwPower.get()

        refPath = self.inputReference.get().getFileName()
        refStruct = self._parseStructure(refPath, 'reference')

        srcCaCoords, srcCaKeys = self._extractCaCoords(refStruct)
        allAtomCoords, allAtoms = self._extractAllAtomCoords(refStruct)

        if srcCaCoords.shape[0] == 0:
            raise RuntimeError(
                'No C_alfa atoms found in reference structure: %s' % refPath
            )

        # Automatic k selection based on mean number of C_alfa within R
        treRref = cKDTree(srcCaCoords)
        counts = [len(treRref.query_ball_point(ca, r=R)) for ca in srcCaCoords]

        if self.maxNeighbours.get() == 0:
            k = max(int(np.mean(counts) * 2), 1)
            self.info('Auto k = %d (mean C_alfa within R=%.1f A: %.1f)'
                      % (k, R, np.mean(counts)))
        else:
            k = self.maxNeighbours.get()

        self.info('Reference: %d C_alfa, %d total atoms, R=%.1f A, k=%d, p=%.2f'
                  % (len(srcCaKeys), len(allAtoms), R, k, p))

        outDir = self._getExtraPath('reconstructed')
        os.makedirs(outDir, exist_ok=True)

        ensemble = self.inputEnsemble.get()
        nModels = len(ensemble)
        self.info('Ensemble: %d model(s).' % nModels)

        io = PDBIO()

        for i, targetAs in enumerate(ensemble):
            targetPath = targetAs.getFileName()
            targetStruct = self._parseStructure(targetPath, 'target_%d' % i)
            dstCaCoords, dstCaKeys = self._extractCaCoords(targetStruct)

            srcAligned, dstAligned = self._alignCaByKey(
                srcCaCoords, srcCaKeys,
                dstCaCoords, dstCaKeys,
            )

            if srcAligned.shape[0] == 0:
                self.warning('No matching C_alfa for model %d - skipping.' % i)
                continue

            self.info('[%d/%d] %d matched C_alfa pairs'
                      % (i + 1, nModels, srcAligned.shape[0]))

            idw = InvDistTree3D(srcAligned, leafsize=10)

            newCoords = self._reconstructAtoms(
                srcCa=srcAligned,
                dstCa=dstAligned,
                allAtoms=allAtomCoords,
                idw=idw,
                R=R,
                k=k,
                power=p
            )

            outStruct = copy.deepcopy(refStruct)
            _, outAtoms = self._extractAllAtomCoords(outStruct)
            self._applyCoordsToStructure(outAtoms, newCoords)

            baseName = os.path.splitext(os.path.basename(targetPath))[0]
            extension = os.path.splitext(targetPath)[1]

            outPath = os.path.join(
                outDir,
                baseName + '_reconstruct' + extension
            )
            
            io.set_structure(outStruct)
            io.save(outPath)
            self._fixPdbFormat(outPath)

        self.info('Reconstruction complete.')

    def createOutputStep(self):
        outDir = self._getExtraPath('reconstructed')
        pdbFiles = sorted(f for f in os.listdir(outDir) if f.endswith('.pdb'))

        if not pdbFiles:
            raise RuntimeError('No reconstructed PDB files in %s' % outDir)

        outputSet = SetOfAtomStructs.create(self._getPath())

        for pdbFile in pdbFiles:
            asObj = AtomStruct()
            asObj.setFileName(os.path.join(outDir, pdbFile))
            outputSet.append(asObj)

        self._defineOutputs(**{OUTPUT_STRUCTUresName: outputSet})
        self._defineSourceRelation(self.inputReference, outputSet)
        self._defineSourceRelation(self.inputEnsemble, outputSet)
        self.info('Output: %d structure(s).' % len(outputSet))

    def _parseStructure(self, path, structId):
        ext = os.path.splitext(path)[1].lower()
        parser = MMCIFParser(QUIET=True) if ext in ('.cif', '.mmcif') \
            else PDBParser(QUIET=True)
        return parser.get_structure(structId, path)

    def _alignCaByKey(self, srcCoords, srcKeys, dstCoords, dstKeys):
        """
        Align C-alpha atoms between reference and target structures.

        Chain IDs are ignored. Reference and target chains are matched based
        on the number of common residue identifiers. Once chains are paired,
        residues are matched by residue number and insertion code.

        This allows a chain to be renamed (e.g. I -> A) without affecting
        the reconstruction.
        """
        srcChains = {}
        dstChains = {}

        for i, key in enumerate(srcKeys):
            chainId, resNum, insCode = key
            srcChains.setdefault(chainId, []).append(i)

        for i, key in enumerate(dstKeys):
            chainId, resNum, insCode = key
            dstChains.setdefault(chainId, []).append(i)

        srcResidues = {}

        for chainId, indices in srcChains.items():
            srcResidues[chainId] = {
                (srcKeys[i][1], srcKeys[i][2])
                for i in indices
            }

        dstResidues = {}

        for chainId, indices in dstChains.items():
            dstResidues[chainId] = {
                (dstKeys[i][1], dstKeys[i][2])
                for i in indices
            }

        unusedDstChains = set(dstChains.keys())
        chainPairs = []

        for srcChainId in srcChains:

            if not unusedDstChains:
                break

            bestDstChain = None
            bestCommon = 0

            for dstChainId in unusedDstChains:

                common = len(
                    srcResidues[srcChainId] &
                    dstResidues[dstChainId]
                )

                if common > bestCommon:
                    bestCommon = common
                    bestDstChain = dstChainId

            if bestDstChain is not None and bestCommon > 0:
                chainPairs.append(
                    (srcChainId, bestDstChain, bestCommon)
                )

                unusedDstChains.remove(bestDstChain)

        srcAligned = []
        dstAligned = []

        matchedSrcChains = set()
        matchedDstChains = set()

        for srcChainId, dstChainId, common in chainPairs:

            matchedSrcChains.add(srcChainId)
            matchedDstChains.add(dstChainId)

            # Target residue -> C-alpha index
            dstMap = {}

            for idx in dstChains[dstChainId]:
                residueKey = (
                    dstKeys[idx][1],
                    dstKeys[idx][2]
                )
                dstMap[residueKey] = idx

            chainMatches = 0

            for srcIdx in srcChains[srcChainId]:

                residueKey = (
                    srcKeys[srcIdx][1],
                    srcKeys[srcIdx][2]
                )

                if residueKey in dstMap:
                    dstIdx = dstMap[residueKey]

                    srcAligned.append(srcCoords[srcIdx])
                    dstAligned.append(dstCoords[dstIdx])

                    chainMatches += 1

            self.info(
                'Matched chain %s -> %s: %d C_alfa'
                % (srcChainId, dstChainId, chainMatches)
            )

        unmatchedSrc = [
            chainId
            for chainId in srcChains
            if chainId not in matchedSrcChains
        ]

        unmatchedDst = [
            chainId
            for chainId in dstChains
            if chainId not in matchedDstChains
        ]

        if unmatchedSrc:
            self.warning(
                'Reference chains with no matching target chain: %s'
                % unmatchedSrc
            )

        if unmatchedDst:
            self.warning(
                'Target chains with no matching reference chain: %s'
                % unmatchedDst
            )

        if not srcAligned:
            return (
                np.empty((0, 3)),
                np.empty((0, 3))
            )

        return (
            np.array(srcAligned, dtype=np.float64),
            np.array(dstAligned, dtype=np.float64)
        )

    def _reconstructAtoms(self, srcCa, dstCa, allAtoms, idw,
                          R=15.0, k=8, power=2.0, leafsize=10):
        srcCa = np.asarray(srcCa, dtype=np.float64)
        dstCa = np.asarray(dstCa, dtype=np.float64)
        allAtoms = np.asarray(allAtoms, dtype=np.float64)

        if srcCa.shape != dstCa.shape:
            raise ValueError(
                "srcCa and dstCa must have the same number of ca, "
                "got %s vs %s" % (srcCa.shape, dstCa.shape)
            )

        displacements = np.asarray(dstCa, dtype=np.float64) - np.asarray(srcCa, dtype=np.float64)
        interpolated = idw(allAtoms, displacements, R=R, k=k, p=power)
        return allAtoms + interpolated

    def _extractCaCoords(self, structure):
        coords, keys = [], []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    if 'CA' not in residue:
                        continue
                    coords.append(residue['CA'].get_vector().get_array())
                    keys.append((chain.id, residue.id[1], residue.id[2].strip()))
            break
        return np.array(coords, dtype=np.float64), keys

    def _extractAllAtomCoords(self, structure):
        coords, atoms = [], []
        for model in structure:
            for chain in model:
                for residue in chain:
                    if residue.id[0] != ' ':
                        continue
                    for atom in residue:
                        coords.append(atom.get_vector().get_array())
                        atoms.append(atom)
            break
        return np.array(coords, dtype=np.float64), atoms

    def _applyCoordsToStructure(self, atoms, newCoords):
        for atom, coord in zip(atoms, newCoords):
            atom.set_coord(coord)

    def _summary(self):
        out = getattr(self, OUTPUT_STRUCTUresName, None)

        if self.isFinished() and out is not None:
            return [
                'Reconstructed %d full-atom structure(s).' % len(out),
                'R=%.1f A, k=%s, p=%.2f.'
                % (
                    self.searchRadius.get(),
                    'auto' if self.maxNeighbours.get() == 0
                    else str(self.maxNeighbours.get()),
                    self.idwPower.get()
                )
            ]

        return ['Protocol not finished yet.']

    def _methods(self):
        return [
            'Full-atom structures were reconstructed from a Ca-only ensemble '
            'using compact-support IDW with KD-tree neighbour search '
            '(R=%.1f A, p=%.2f). '
            'Formula: w_k = ((R - d(x,x_k)) / (R * d(x,x_k)))^p.'
            % (
                self.searchRadius.get(),
                self.idwPower.get()
            )
        ]

    def _validate(self):
        errors = []
        if self.inputReference.get() is None:
            errors.append('A reference full-atom structure must be provided.')
        ens = self.inputEnsemble.get()
        if ens is None:
            errors.append('A target Ca ensemble must be provided.')
        elif len(ens) == 0:
            errors.append('The target ensemble is empty.')
        if self.searchRadius.get() <= 0:
            errors.append('Search radius R must be positive.')
        if self.idwPower.get() <= 0:
            errors.append('IDW power p must be positive.')
        return errors

    @staticmethod
    def _fixPdbFormat(pdbPath):
        """Strips original TER/END lines and regenerates them. a TER is
        inserted right before each chain's first HETATM"""

        with open(pdbPath, 'r') as f:
            lines = f.readlines()

        fixed = []
        atomCount = 0
        lastAtom = None
        inHetatm = False

        for line in lines:
            record = line[:6].strip()

            if record == 'ATOM':
                inHetatm = False
                atomCount += 1
                lastAtom = line
                fixed.append(line)

            elif record == 'HETATM':
                if not inHetatm and lastAtom is not None:
                    atomCount += 1
                    resName = lastAtom[17:20]
                    chain = lastAtom[21]
                    resseq = lastAtom[22:26]
                    fixed.append('TER   %5d      %3s %s%4s\n' % (
                        atomCount, resName, chain, resseq
                    ))
                    inHetatm = True
                atomCount += 1
                fixed.append(line)

            elif record in ('TER', 'END'):
                continue

            else:
                fixed.append(line)

        if not inHetatm and lastAtom is not None:
            atomCount += 1
            resName = lastAtom[17:20]
            chain = lastAtom[21]
            resseq = lastAtom[22:26]
            fixed.append('TER   %5d      %3s %s%4s\n' % (
                atomCount, resName, chain, resseq
            ))

        fixed.append('END\n')

        with open(pdbPath, 'w') as f:
            f.writelines(fixed)