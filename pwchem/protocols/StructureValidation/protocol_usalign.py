# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


from pathlib import Path
from pwem.convert import cifToPdb

from pyworkflow.object import Float
from pyworkflow.protocol import params, STEPS_PARALLEL
from pwem.protocols import EMProtocol

from pwem.objects import AtomStruct
from pwchem.utils import os
import re
from pwchem import Plugin, USALIGN_DIC


class ProtocolUSalign(EMProtocol):
    """
    AI Generated:
    """
    _label = 'US-align structures'
    stepsExecutionMode = params.STEPS_PARALLEL
    program = os.path.join(Plugin.getVar(USALIGN_DIC['home']), 'USalign/USalign')

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure1', params.PointerParam, pointerClass='AtomStruct',
                        label='Input structure: ',
                        help='Input structure to align.')
        form.addParam('chain1',
            params.StringParam, default='', label='Structure 1 chains',
            expertLevel=params.LEVEL_ADVANCED,
            help='Chains to include form the first structure'
        )
        form.addParam('inputStructure2', params.PointerParam, pointerClass='AtomStruct',
                        label='Input reference structure: ',
                        help='Input reference structure to align to.')
        form.addParam('chain2',
                        params.StringParam, default='', label='Structure 2 chains',
                        expertLevel=params.LEVEL_ADVANCED,
                        help='Chains to include form the second structure'
        )
        group = form.addGroup('Parameters')
        group.addParam('moleculeType',params.EnumParam,choices=['Protein and DNA/RNA', 'Protein', 'RNA/DNA'],
            default=0, label='Molecule type: ',
            help='Type of molecules to align.'
        )
        group.addParam('alignmentMode', params.EnumParam, choices=['Monomer', 'Complex (multimer), Biological assembly'],
            default=0, label='Alignment mode: ', help='Multimeric alignment options.'
        )
        group.addParam('fastMode', params.BooleanParam, default=False,
            label='Fast alignment. May be worse quality alignment if selected.'
        )
        group.addParam('includeHetatm', params.BooleanParam, default=False,
            label='Include HETATM residues'
        )

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.runUSalignStep)
        self._insertFunctionStep(self.createOutputStep)

    def runUSalignStep(self):
        struct1 = self.inputStructure1.get().getFileName()
        struct2 = self.inputStructure2.get().getFileName()

        if self.moleculeType.get() == 0:
            molType = 'auto'
        elif self.moleculeType.get() == 1:
            molType = 'prot'
        else:
            molType = 'rna'

        ter = self.getTerParam()

        args = [os.path.abspath(struct1),
                os.path.abspath(struct2),
                '-mol', molType,
                '-mm', self.alignmentMode.get(),
                '-o', 'superposed'
                ]
        if ter is not None:
            args += ['-ter', ter]

        if self.fastMode .get(): args.append('-fast')
        if self.includeHetatm .get(): args += ['-het', 1]

        summaryFile = "summary.txt"

        cmd = (
            f'{self.program} '
            f'{os.path.abspath(struct1)} '
            f'{os.path.abspath(struct2)} '
            f'-mol {molType} '
            f'-mm {self.alignmentMode.get()} '
            f'-o superposed '
            f'2>&1 | tee {summaryFile}'
        )

        self.runJob(
            "/bin/bash",
            f'-c "{cmd}"',
            cwd=self._getExtraPath(),
            env=Plugin.getEnviron()
        )

    def createOutputStep(self):
        outStruct = AtomStruct()
        outStruct.setFileName(self._getExtraPath("superposed.cif"))

        self._defineOutputs(outputStructure=outStruct)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summaryFile = self._getExtraPath("summary.txt")

        if not os.path.exists(summaryFile):
            return ["US-align has not been run yet."]

        rmsd = None
        alnLength = None
        seqId = None
        tmReference = None
        tmInput = None

        with open(summaryFile) as f:
            for line in f:

                if line.startswith("Aligned length="):
                    m = re.search(
                        r'Aligned length=\s*(\d+), RMSD=\s*([\d.]+), .*Seq_ID=.*=\s*([\d.]+)',
                        line
                    )
                    if m:
                        alnLength = m.group(1)
                        rmsd = m.group(2)
                        seqId = m.group(3)

                elif line.startswith("TM-score="):
                    m = re.search(
                        r'TM-score=\s*([\d.]+).*Structure_(\d)',
                        line
                    )

                    if m:
                        score = m.group(1)

                        if m.group(2) == "2":
                            tmReference = score
                        else:
                            tmInput = score

        summary = [
            f"Reference structure: {os.path.basename(self.inputStructure2.get().getFileName())}"
        ]

        if tmReference:
            summary.append(f"TM-score (normalized by reference): {tmReference}")

        if tmInput:
            summary.append(f"TM-score (normalized by input): {tmInput}")

        if rmsd:
            summary.append(f"RMSD: {rmsd} Å")

        if alnLength:
            summary.append(f"Aligned length: {alnLength} residues")

        if seqId:
            summary.append(f"Sequence identity: {float(seqId) * 100:.1f}%")

        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        validations = []
        return validations

    def _warnings(self):
        warnings = []
        return warnings

    # --------------------------- UTILS functions -----------------------------------
    def getTerParam(self):
        if self.alignmentMode.get() == 0:
            ter = None
        elif self.alignmentMode.get() == 1:
            ter = 1
        else:
            ter = 0
        return ter


