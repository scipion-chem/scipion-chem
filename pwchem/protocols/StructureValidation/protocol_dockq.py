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

    This protocol evaluates the quality of a predicted protein-protein or
    protein-nucleic acid complex by comparing it against a reference complex
    using the DockQ scoring framework.

    DockQ combines several interface quality metrics into a single continuous
    score that correlates with CAPRI quality classifications, making it a
    standard metric for assessing docking and complex structure predictions.

    Core Concepts
    -------------
    DockQ:
        Composite score ranging from 0 to 1 that combines interface RMSD,
        ligand RMSD, and the fraction of native contacts into a single measure
        of docking quality. Higher values indicate better agreement with the
        reference complex.

    iRMSD (Interface RMSD):
        Root-mean-square deviation computed using interface residues after
        optimal superposition. Measures the accuracy of the binding interface.

    LRMSD (Ligand RMSD):
        RMSD of the ligand chains after superposing the receptor. Measures how
        accurately the ligand is positioned relative to the receptor.

    Fnat:
        Fraction of native residue-residue contacts reproduced by the predicted
        complex.

    F1 Score:
        Harmonic mean of interface precision and recall, evaluating how well
        interface contacts are recovered.

    Chain Mapping:
        Defines the correspondence between chains in the predicted and
        reference complexes. DockQ can automatically determine the mapping or
        use a user-provided mapping.

    Workflow
    --------
    1. Load the predicted complex.
    2. Load the reference (native) complex.
    3. Optionally define chain mappings between both structures.
    4. Optionally disable structural alignment and compare structures using
       residue numbering.
    5. Execute DockQ.
    6. Parse the reported quality metrics.
    7. Store the DockQ score and interface RMSD as attributes of the output
       structure.

    Input
    -----
    - inputStructure1:
        Predicted complex to evaluate.

    - inputStructure2:
        Reference (native) complex.

    Parameters
    ----------
    - Chain mapping:
        Optional mapping between chains of the predicted and reference
        complexes. This can also restrict the evaluation to specific
        interfaces.

    - Align:
        Whether DockQ should structurally align both complexes before
        comparison. If disabled, residue numbering is used directly.

    - Allowed sequence mismatches:
        Maximum number of mismatches permitted while mapping sequences between
        model and reference.

    - Optimize DockQ_F1:
        Optimizes the interface F1 score instead of the standard DockQ score.

    Output
    ------
    - outputStructure:
        Structure associated with the docking evaluation.

        The output contains the following attributes:

        - DockQ:
            Overall docking quality score.

        - iRMSD:
            Interface RMSD between predicted and reference complexes.

    Reported Metrics
    ----------------
    The protocol extracts and summarizes:

    - Overall DockQ score
    - Chain mapping used during evaluation
    - Interface DockQ score
    - Interface RMSD (iRMSD)
    - Ligand RMSD (LRMSD)
    - Fraction of native contacts (Fnat)
    - Interface contact F1 score
    - Number of interface clashes

    Use Cases
    ---------
    - Evaluating protein-protein docking predictions
    - Benchmarking AlphaFold-Multimer or similar complex prediction methods
    - Comparing predicted complexes against experimental structures
    - Assessing docking quality in CAPRI-style evaluations
    - Measuring interface accuracy after docking refinement

    Notes
    -----
    - DockQ provides a more comprehensive assessment than RMSD alone by
      combining geometric and contact-based measures into a single score.

    - When multiple interfaces are present, DockQ reports an overall score as
      well as interface-specific metrics, depending on the evaluated complex.
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
        self._insertFunctionStep(self.createOutputStep)

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
        dockQ, rmsd = self.getScores()

        outStruct.DockQ = Float()
        outStruct.setAttributeValue('DockQ', dockQ)

        outStruct.iRMSD = Float()
        outStruct.setAttributeValue('iRMSD', rmsd)

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
    def getScores(self):
        summaryFile = self._getExtraPath("summary.txt")
        outDockQ =None
        outiRMSD =None
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
                elif line.startswith("DockQ:"):
                    dockQ = line.split(":")[1].strip()
                elif line.startswith("iRMSD:"):
                    iRMSD = line.split(":")[1].strip()
        if totalDockQ: outDockQ = totalDockQ
        else: outDockQ = dockQ

        return outDockQ, outiRMSD


