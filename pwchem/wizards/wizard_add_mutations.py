# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Judith Maestro Ciria
# *
# * Biocomputing Unit, CNB-CSIC
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

import json

from pwem.wizards import EmWizard
import pwem.convert as emconv

from pwchem.utils import RESIDUES3TO1, MODIFIED_RESIDUES3TO1

AA_THREE_TO_ONE = {**RESIDUES3TO1, **MODIFIED_RESIDUES3TO1}


class AddMutationsWizard(EmWizard):
    def getPositions(self, form):
        protocol = form.protocol
        allRanPos = protocol.RangPositions.get().split(", ")
        return allRanPos

    def getaaTo(self, form):
        protocol = form.protocol
        return "X" if protocol.mutSaturation else str(protocol.mutResidue.get())

    def getchainResidues(self, form):
        protocol = form.protocol
        structureHandler = emconv.AtomicStructHandler()
        structureHandler.read(protocol.inputAtomStruct.get().getFileName())
        structureHandler.getStructure()
        _, modelsFirstResidue = structureHandler.getModelsChains()
        chainResidues = {}

        for modelID, chains in modelsFirstResidue.items():
            for chainID, residues in chains.items():
                if chainID not in chainResidues:
                    chainResidues[chainID] = {}
                for residue in residues:
                    resId = residue[0]
                    resType = residue[1]
                    chainResidues[chainID][resId] = resType

        return chainResidues

    def getchain(self, form):
        protocol = form.protocol
        chain = json.loads(protocol.mutChain.get())['chain']
        return chain

    def getROIOrigen(self, form):
        protocol = form.protocol
        roiOrigin = protocol.ROIOrigin.get()
        return roiOrigin

    def getSructROI(self, form):
        protocol = form.protocol
        inputStructROI = protocol.inputStructROI.get()
        return inputStructROI

    def getROIChain(self, form):
        protocol = form.protocol
        chainStr = protocol.ROIChain.get()
        if chainStr and chainStr.strip():
            return json.loads(chainStr)['chain']
        return None

    def _getMutationsFromRange(self, form, chainResidues, aaTo):
        chain = self.getchain(form)
        residuesDict = chainResidues.get(chain, {})
        mutations = []
        seen = set()
        for ranPos in self.getPositions(form):
            first, last = ranPos.split("-")
            for pos in range(int(first), int(last) + 1):
                if pos in residuesDict:
                    aaFrom = AA_THREE_TO_ONE[residuesDict[pos]]
                    mutation = '{}{}{}{}'.format(aaFrom, chain, pos, aaTo)
                    if mutation not in seen:
                        seen.add(mutation)
                        mutations.append(mutation)
        return mutations

    def _getMutationsFromROI(self, form, chainResidues, aaTo):
        roiChain = self.getROIChain(form)
        mutations = []
        seen = set()
        for item in self.getSructROI(form):
            for roi in item.getDecodedCResidues():
                chain, pos = roi.split("_")
                pos = int(pos)
                if roiChain and chain != roiChain:
                    continue
                residuesDict = chainResidues.get(chain, {})
                if pos in residuesDict:
                    aaFrom = AA_THREE_TO_ONE[residuesDict[pos]]
                    mutation = '{}{}{}{}'.format(aaFrom, chain, pos, aaTo)
                    if mutation not in seen:
                        seen.add(mutation)
                        mutations.append(mutation)
        return mutations

    def getMutations(self, form):
        aaTo = self.getaaTo(form)
        chainResidues = self.getchainResidues(form)
        roiOrigin = self.getROIOrigen(form)

        if roiOrigin == 0:
            return self._getMutationsFromRange(form, chainResidues, aaTo)
        return self._getMutationsFromROI(form, chainResidues, aaTo)

    def show(self, form, *params):
        protocol = form.protocol
        mutations = self.getMutations(form)

        toMutateList = protocol.toMutateList.get()
        existing = {line.strip() for line in toMutateList.strip().split("\n") if line.strip()}
        newMutations = [m for m in mutations if m not in existing]
        toMutateList += "\n" + "\n".join(newMutations)
        form.setVar('toMutateList', toMutateList.strip())


class ClearMutationsWizard(EmWizard):
  def show(self, form, *params):
    form.setVar('toMutateList', '')
