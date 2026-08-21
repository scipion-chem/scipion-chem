# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
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

"""
This protocol combines the results from several ROI sets and counts how many
times each residue appears across models (one count per model, not per pocket).
It outputs a SetOfStructROIs and a SequenceChem with per-residue frequency attributes.
"""

import pyworkflow.object as pwobj
from pwem.protocols import EMProtocol
from pwem.objects import AtomStruct
from pwem.convert.atom_struct import toCIF, AtomicStructHandler, addScipionAttribute
from pyworkflow.protocol import params
from pwchem.objects import SetOfStructROIs, StructROI, SequenceChem
from pwchem.utils import getBaseName, createPocketFile, parseResidueCoords


class ProtROIVoting(EMProtocol):
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