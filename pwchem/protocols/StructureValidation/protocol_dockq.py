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
from pwchem import Plugin, DOCKQ_DIC

class ProtocolDockQ(EMProtocol):
    """
    AI Generated:
    """
    _label = 'DockQ quality measure'
    stepsExecutionMode = params.STEPS_PARALLEL

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure1', params.PointerParam, pointerClass='AtomStruct',
                        label='Input structure: ',
                        help='Input structure to align.')
        form.addParam('inputStructure2', params.PointerParam, pointerClass='AtomStruct',
                        label='Input reference structure: ',
                        help='Input reference structure to align to.')

        group = form.addGroup('Parameters')
        #group.addParam('moleculeType',params.EnumParam,choices=['Protein or DNA/RNA', 'Small molecule docking'],
        #    default=0, label='Evaluation type: ',
        #    help='Type of molecules to measure.'
        #)
        group.addParam('mapping',params.StringParam,default='',expertLevel=params.LEVEL_ADVANCED,
            label='Chain mapping',
            help='Specify a chain mapping between model and native structure. If the native contains two chains "H" and "L" while the model contains two chains "A" and "B", and chain A is a model of native chain H and chain B is a model\n '
                'of native chain L, the flag can be set as: "--mapping AB:HL". This can also help limit the search to specific native interfaces. For example, if the native is a tetramer (ABCD) but the user is only interested in the \n'
                'interface between chains B and C, the flag can be set as: "--mapping :BC" or the equivalent "--mapping *:BC".'
        )
        group.addParam('align', params.BooleanParam, default=True, expertLevel=params.LEVEL_ADVANCED,
            label='Align: ',
            help='Do not align, use residue numbering.'
        )
        group.addParam('allowedMismatches', params.IntParam, default=5,
            label='Allowed sequence mismatches: ',
            help='Number of allowed mismatches when mapping model sequence to native sequence.'
        )
        group.addParam( 'optDockQF1', params.BooleanParam, default=False,
            label='Optimize DockQ_F1: ',
            help=' Optimize on DockQ_F1 instead of DockQ.'
        )

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.runDockQStep)
        #self._insertFunctionStep(self.createOutputStep)

    def runDockQStep(self):
        struct1 = os.path.abspath(self.inputStructure1.get().getFileName())
        struct2 = os.path.abspath(self.inputStructure2.get().getFileName())

        cmd = f'DockQ "{struct1}" "{struct2}"'

        if self.mapping.get():
            cmd += f' --mapping "{self.mapping.get()}"'

        if not self.align.get():
            cmd += ' --no_align'

        cmd += f' --allowed_mismatches {self.allowedMismatches.get()}'

        if self.optDockQF1.get():
            cmd += ' --optDockQF1'

        summaryFile = os.path.abspath(self._getExtraPath("summary.txt"))

        cmd += f' 2>&1 | tee {summaryFile}'

        self.runJob(
            Plugin.getEnvActivationCommand(DOCKQ_DIC),
            f'&& {cmd}',
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
            return ["DockQ has not been run yet."]

        totalDockQ = None
        mapping = None
        dockQ = None
        iRMSD = None
        lRMSD = None
        fnat = None
        f1 = None
        clashes = None
        with open(summaryFile) as f:
            for line in f:
                line = line.strip()
                if line.startswith("Total DockQ"):
                    m = re.search(
                        r"Total DockQ over (\d+) native interfaces: ([\d.]+) with (.+) model:native mapping",
                        line
                    )
                    if m:
                        totalDockQ = m.group(2)
                        mapping = m.group(3)
                elif line.startswith("DockQ:"):
                    dockQ = line.split(":")[1].strip()
                elif line.startswith("iRMSD:"):
                    iRMSD = line.split(":")[1].strip()
                elif line.startswith("LRMSD:"):
                    lRMSD = line.split(":")[1].strip()
                elif line.startswith("fnat:"):
                    fnat = line.split(":")[1].strip()
                elif line.startswith("F1:"):
                    f1 = line.split(":")[1].strip()
                elif line.startswith("clashes:"):
                    clashes = line.split(":")[1].strip()
        summary = []
        if totalDockQ:
            summary.append(f"Total DockQ score: {totalDockQ}")
        if mapping:
            summary.append(f"Chain mapping: {mapping}")
        if dockQ:
            summary.append(f"Interface DockQ score: {dockQ}")
        if iRMSD:
            summary.append(f"Interface RMSD (iRMSD): {iRMSD} Å")
        if lRMSD:
            summary.append(f"Ligand RMSD (LRMSD): {lRMSD} Å")
        if fnat:
            summary.append(f"Fraction of native contacts (Fnat): {fnat}")
        if f1:
            summary.append(f"Interface contact F1 score: {f1}")
        if clashes:
            summary.append(f"Interface clashes: {clashes}")

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


