# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
# *          Judith Maestro (judith.maestro@cnb.csic.es)
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pwem.convert.atom_struct import toCIF, AtomicStructHandler, addScipionAttribute
from pyworkflow.protocol import params
from pwchem.objects import SetOfStructROIs, StructROI, SequenceChem
from pwchem.utils import getBaseName, createPocketFile, parseResidueCoords


class ProtROIVoting(EMProtocol):
    """
    Protocol to combine several SetOfStructROIs by residue voting.

    AI Generated:

        ProtROIVoting - User Manual

        Overview
        --------
        The ProtROIVoting protocol integrates the results of several independent
        structural ROI (region of interest) detections performed on the same
        protein structure -e.g. repeated pocket-detection runs, docking-derived
        contact residues, or ROIs defined on different structural models- into
        a single consensus result based on residue voting.

        For each input SetOfStructROIs ("model"), the protocol collects the set
        of protein residues in contact with at least one of its ROIs (one vote
        per model, regardless of how many ROIs/pockets in that model contain the
        residue). Residues are then ranked by how many models voted for them,
        and this frequency is embedded in a new SetOfStructROIs, in per-chain
        sequence attributes, and in the original structure so it can be
        inspected in ChimeraX.

        Inputs
        ------
        - **roisList**:
            A list of SetOfStructROIs objects (pointers), each one typically
            coming from a different ROI-detection run (e.g. several
            ProtDefineStructROIs / ProtocolConsensusStructROIs executions, or
            pockets computed on different conformations/models) over the *same*
            protein structure. All the input sets must describe the same
            structure (same chains and residue numbering), since residues are
            compared directly by chain and residue number across sets.

        Workflow
        --------
        1. **Voting**:
           - For every input SetOfStructROIs, the contact residues of all its
             ROIs are decoded and merged into a single set (so a residue counts
             once per model, even if it appears in several of that model's
             pockets).
           - A global counter tallies, for every residue, in how many of the
             input models it appeared.

        2. **Structural ROI reconstruction**:
           - For every voted residue, its atom coordinates are looked up in the
             (shared) protein structure and written out as a single-residue
             StructROI, with its contacts recalculated on the actual structure.
           - ROIs are numbered by decreasing vote count.
           - Each ROI stores its raw vote count (`_frequency`) and its
             percentage relative to the most-voted residue (`_percentage`).

        3. **Per-chain sequence output**:
           - For every protein chain that has at least one voted residue, a
             SequenceChem is created with a `frequency` per-residue attribute
             (0 for residues that were never voted).

        4. **Structure attribute output**:
           - The `frequency` value of every residue in the structure (0 if it
             was never voted) is embedded as a Scipion attribute in a copy of
             the protein structure, so it can be visualized in 3D.

        Outputs
        -------
        - **outputStructROIs**:
            SetOfStructROIs containing one ROI per voted residue, ordered by
            descending vote count, each with `_frequency` and `_percentage`
            attributes.

        - **outputSequence_<chainId>** (one per chain with voted residues):
            SequenceChem of that chain with a `frequency` per-residue attribute.

        - **outputAtomStruct**:
            AtomStruct (cif) of the input protein with the per-residue
            `frequency` embedded as a Scipion attribute, meant to be explored
            with the ROI Voting viewer (3D coloring in ChimeraX, histogram and
            sequence plot of the voting frequency).

        Practical Recommendations
        -------------------------
        - Use this protocol to find the contact residues/pockets that are
          consistently detected across several independent runs (different
          pocket-detection algorithms, different conformations/models, or
          repeated stochastic predictions), rather than trusting a single run.
        - All inputs should come from the same protein (same chains/numbering);
          mixing SetOfStructROIs from different structures will produce
          meaningless votes.
        - Open the "Viewer ROI Voting" on this protocol to inspect the results:
          color the structure by voting frequency in ChimeraX, or plot/inspect
          the frequency along the sequence.

        Interpretation
        --------------
        The `_frequency`/`frequency` value of a residue is the number of input
        models in which that residue was in contact with at least one ROI.
        `_percentage` normalizes this count so that the most-voted residue is
        100%. Residues with high frequency/percentage are the ones most
        consistently identified as relevant across the different inputs.

        Warnings
        --------
        - Input SetOfStructROIs describing different structures (or the same
          structure with different chain/residue numbering) will not be
          compared correctly, since voting matches residues by chain and
          residue number.
        - If none of the inputs report any contact residue, no outputs are
          generated.

        Final Perspective
        -----------------
        ProtROIVoting provides a simple consensus mechanism to combine several
        structural ROI detections into a single, ranked, and visualizable set of
        the most consistently identified residues, complementing more complex
        consensus approaches (e.g. ProtocolConsensusStructROIs) with a
        lightweight, residue-level voting scheme.
    """

    _label = 'ROI voting (multi-model integration)'
    _ATTRNAME = 'frequency'
    _OUTNAME = 'outputAtomStruct'

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('roisList', params.MultiPointerParam,
                      label='Input Sets of Structural ROIs',
                      pointerClass='SetOfStructROIs',
                      help='Select all SetOfStructROIs to combine.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.computeVotingStep)

    def computeVotingStep(self):
        residueCounts = {}
        proteinFile = None

        for p in self.roisList:
            roiSet = p.get()
            if proteinFile is None:
                proteinFile = roiSet.getProteinFile()

            modelResidues = set()
            for roi in roiSet:
                modelResidues.update(roi.getDecodedCResidues())

            for r in modelResidues:
                residueCounts[r] = residueCounts.get(r, 0) + 1

        if not residueCounts:
            self.warning("No contact residues found in the input ROI sets.")
            return

        maxCount = max(residueCounts.values())
        sortedResidues = sorted(residueCounts.items(), key=lambda x: x[1], reverse=True)
        resCoordsDic = parseResidueCoords(proteinFile)

        outSet = SetOfStructROIs(filename=self._getPath('ROIVoting.sqlite'))
        for i, (residue, count) in enumerate(sortedResidues):
            coords = resCoordsDic.get(residue)
            if not coords:
                self.warning(f"Residue {residue} not found in {proteinFile}, skipping")
                continue

            pocketFile = self._getExtraPath(f'roi_{residue}.cif')
            createPocketFile(coords, i + 1, pocketFile)

            roi = StructROI(pocketFile, proteinFile)
            roi.setObjId(i + 1)
            roi.calculateContacts()
            roi._frequency = pwobj.Integer(count)
            roi._percentage = pwobj.Float(round((count / maxCount) * 100.0, 2))
            outSet.append(roi)

        if len(outSet) > 0:
            outSet.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outSet)
            self._defineSourceRelation(self.roisList, self.outputStructROIs)

            self.createSequenceOutputs(outSet, residueCounts, maxCount)
            self.createAtomStructOutput(proteinFile, residueCounts, resCoordsDic)

    def createAtomStructOutput(self, proteinFile, residueCounts, resCoordsDic):
        '''Embeds the per-residue voting frequency as a Scipion attribute in the protein cif file,
        so it can be inspected with the generic ChimeraAttributeViewer (3D coloring, histogram...).
        Every residue in the structure gets a value (0 if not voted): pwem's replaceOcuppancyWithAttribute
        crashes on residues missing from the attribute dict, so it must not be left sparse.'''
        attrDic = {residue.replace('_', ':', 1): 0 for residue in resCoordsDic}
        attrDic.update({residue.replace('_', ':', 1): count for residue, count in residueCounts.items()})

        ash = AtomicStructHandler()
        inpCif = toCIF(proteinFile, self._getTmpPath('inputStruct.cif'))
        cifDic = addScipionAttribute(ash.readLowLevel(inpCif), attrDic, self._ATTRNAME)

        outStructFileName = self._getPath('outputStructureVoting.cif')
        ash._writeLowLevel(outStructFileName, cifDic)
        self._defineOutputs(outputAtomStruct=AtomStruct(filename=outStructFileName))

    def createSequenceOutputs(self, roiSet, residueCounts, maxCount):
        '''Defines one SequenceChem output (with a per-residue frequency attribute) for each
        protein chain that has at least one voted contact residue'''
        seqDic = roiSet.getProteinSequencesDic()
        resIdxDic = roiSet.getProteinSequencesResIdsDic()
        base = getBaseName(roiSet.getProteinFile())

        chainsWithROIs = sorted({residue.rsplit('_', 1)[0] for residue in residueCounts})
        for chainId in chainsWithROIs:
            if chainId not in seqDic:
                continue

            seqStr = seqDic[chainId]
            freqValues = [0] * len(seqStr)
            for resId, seqIdx in resIdxDic[chainId].items():
                # Trailing non-standard residues (waters, heteroatoms) are counted here but
                # trimmed from seqStr by getSequenceFromChain, so they fall outside its range
                if seqIdx < len(seqStr):
                    freqValues[seqIdx] = residueCounts.get(f'{chainId}_{resId}', 0)

            outSeq = SequenceChem(
                name=f'{base}_{chainId}',
                sequence=seqStr,
                id=f'{base}_{chainId}',
                attributesFile=self._getExtraPath(f'sequenceAttributes_{chainId}.txt')
            )
            outSeq.addAttributes({'frequency': freqValues})
            self._defineOutputs(**{f'outputSequence_{chainId}': outSeq})

    def _summary(self):
        summary = []
        if hasattr(self, 'outputStructROIs'):
            summary.append(f"Voting results: {self.outputStructROIs.getSize()} residues.")
        else:
            summary.append("No residues detected or protocol not executed yet.")
        return summary

    def _methods(self):
        return ["Residues are scored by frequency across ROI runs (one vote per model), "
                "and normalized so that the highest count equals 100%."]