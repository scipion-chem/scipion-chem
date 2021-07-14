# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam
import pyworkflow.object as pwobj
from pwchem import Plugin as pwchem_plugin
from pyworkflow.utils.path import makePath, createLink, cleanPattern
from pwchem.objects import SetOfSmallMolecules, SmallMolecule

class ProtBioinformaticsAutodock(EMProtocol):
    """Perform a docking experiment with autodock. See the help at
       http://autodock.scripps.edu/faqs-help/manual/autodock-4-2-user-guide/AutoDock4.2_UserGuide.pdf"""
    _label = 'autodock'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputGrid', PointerParam, pointerClass="AutodockGrid",
                       label='Input grid:', allowsNull=False,
                       help="The grid must be prepared for autodock")
        form.addParam('inputLibrary', PointerParam, pointerClass="SetOfSmallMolecules",
                       label='Ligand library:', allowsNull=False,
                       help="The library must be prepared for autodock")
        form.addParam('rmsTol', FloatParam, label='Cluster tolerance (A)', default=2.0)

        form.addSection(label="Genetic algorithm")
        form.addParam('gaPop', IntParam, label='Population size', default=150)
        form.addParam('gaNumEvals', IntParam, label='Number of evaluations', default=2500000)
        form.addParam('gaNumGens', IntParam, label='Number of generations', default=27000)
        form.addParam('gaElitism', IntParam, label='Elitism', default=1)
        form.addParam('gaMutationRate', FloatParam, label='Mutation rate', default=0.02)
        form.addParam('gaCrossOverRate', FloatParam, label='Crossover rate', default=0.8)
        form.addParam('gaWindowSize', IntParam, label='Window size', default=10)
        form.addParam('lsFreq', FloatParam, label='Local search frequency', default=0.06)
        form.addParam('gaRun', IntParam, label='Number of runs', default=10)

        form.addSection(label="Solis & Wets algorithm")
        form.addParam('swMaxIts', IntParam, label='Max. Number of iterations', default=300)
        form.addParam('swMaxSucc', IntParam, label='Successes in a raw', default=4,
                      help='Number of successes before changing rho')
        form.addParam('swMaxFail', IntParam, label='Failures in a raw', default=4,
                      help='Number of failures before changing rho')
        form.addParam('swRho', FloatParam, label='Initial variance', default=1.0,
                      help='It defines the size of the local search')
        form.addParam('swLbRho', FloatParam, label='Variance lower bound', default=0.01)
        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        fnGridDir = self.inputGrid.get().getFileName()
        dockSteps = []
        for smallMol in self.inputLibrary.get():
            fnSmall = smallMol.getFileName()

            stepId = self._insertFunctionStep('dockStep', fnGridDir, fnSmall, prerequisites=[])
            dockSteps.append(stepId)
        self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

    def dockStep(self, fnGridDir, fnSmall):
        fnReceptor = os.path.join(fnGridDir,"atomStruct.pdbqt")
        fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
        fnSmallDir = self._getExtraPath(fnBase)
        makePath(fnSmallDir)
        fnDPF = os.path.join(fnSmallDir,fnBase+".dpf")
        args = " -l %s -r %s -o %s"%(fnSmall, fnReceptor, fnDPF)

        args += " -p ga_pop_size=%d"%self.gaPop.get()
        args += " -p ga_num_evals=%d"%self.gaNumEvals.get()
        args += " -p ga_num_generations=%d"%self.gaNumGens.get()
        args += " -p ga_elitism=%d"%self.gaElitism.get()
        args += " -p ga_mutation_rate=%f"%self.gaMutationRate.get()
        args += " -p ga_crossover_rate=%f"%self.gaCrossOverRate.get()
        args += " -p ga_window_size=%d"%self.gaWindowSize.get()
        args += " -p sw_max_its=%d"%self.swMaxIts.get()
        args += " -p sw_max_succ=%d"%self.swMaxSucc.get()
        args += " -p sw_max_fail=%d"%self.swMaxFail.get()
        args += " -p sw_rho=%f"%self.swRho.get()
        args += " -p sw_lb_rho=%d"%self.swLbRho.get()
        args += " -p ls_search_freq=%f"%self.lsFreq.get()
        args += " -p ga_run=%d"%self.gaRun.get()
        args += " -p rmstol=%f"%self.rmsTol.get()

        self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                    pwchem_plugin.getADTPath('Utilities24/prepare_dpf42.py')+args)

        fnSmallLocal = os.path.split(fnSmall)[1]
        createLink(fnSmall,os.path.join(fnSmallDir,fnSmallLocal))
        createLink(fnReceptor,os.path.join(fnSmallDir,"atomStruct.pdbqt"))

        args = " -r atomStruct.pdbqt -l %s -o library.gpf"%fnSmallLocal
        self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                    pwchem_plugin.getADTPath('Utilities24/prepare_gpf4.py') + args,
                    cwd=fnSmallDir)

        args = "-p library.gpf -l library.glg"
        self.runJob(pwchem_plugin.getAutodockPath("autogrid4"), args, cwd=fnSmallDir)

        args = "-p %s.dpf -l %s.dlg"%(fnBase,fnBase)
        self.runJob(pwchem_plugin.getAutodockPath("autodock4"), args, cwd=fnSmallDir)

        # Clean a bit
        cleanPattern(os.path.join(fnSmallDir,"atomStruct.*.map"))

    def createOutputStep(self):
        outputSetBest = SetOfSmallMolecules().create(path=self._getPath(),suffix='Best')
        outputSet = SetOfSmallMolecules().create(path=self._getPath())
        for smallMol in self.inputLibrary.get():
            fnSmall = smallMol.getFileName()
            fnBase = os.path.splitext(os.path.split(fnSmall)[1])[0]
            fnSmallDir = self._getExtraPath(fnBase)
            fnDlg = os.path.join(fnSmallDir,fnBase+".dlg")
            if os.path.exists(fnDlg):
                args = " -d . -b -o bestDock.txt"
                self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                            pwchem_plugin.getADTPath('Utilities24/summarize_results4.py') + args,
                            cwd=fnSmallDir)
                args = " -d . -o bestCluster.txt"
                self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                            pwchem_plugin.getADTPath('Utilities24/summarize_results4.py') + args,
                            cwd=fnSmallDir)
                args = " -f %s.dlg -o %s_top.pdbqt"%(fnBase,fnBase)
                self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                            pwchem_plugin.getADTPath('Utilities24/write_lowest_energy_ligand.py') + args,
                            cwd=fnSmallDir)

                newSmallMol = SmallMolecule()
                newSmallMol.copy(smallMol)
                fh = open(os.path.join(fnSmallDir,"bestDock.txt"))
                lineNo = 0
                for line in fh.readlines():
                    if lineNo == 1:
                        tokens = line.split(',')
                        newSmallMol.dockingScoreLE = pwobj.Float(tokens[4].strip())
                        newSmallMol.ligandEfficiency = pwobj.Float(tokens[-1].strip())
                        newSmallMol.smallMoleculeFilePose = pwobj.String(os.path.join(fnSmallDir,"%s_top.pdbqt"%fnBase))
                    lineNo += 1
                fh.close()
                outputSetBest.append(newSmallMol)

                fh = open(os.path.join(fnSmallDir, "bestCluster.txt"))
                lineNo = 0
                for line in fh.readlines():
                    if lineNo >= 1:
                        newSmallMol = SmallMolecule()
                        newSmallMol.copy(smallMol)
                        newSmallMol.cleanObjId()
                        tokens = line.split(',')
                        newSmallMol.dockingScoreLE = pwobj.Float(tokens[2].strip())
                        newSmallMol.ligandEfficiency = pwobj.Float(tokens[-1].strip())
                        outputSet.append(newSmallMol)
                    lineNo += 1
                fh.close()

        self._defineOutputs(outputSmallMolecules=outputSet)
        self._defineSourceRelation(self.inputGrid, outputSet)
        self._defineSourceRelation(self.inputLibrary, outputSet)

        self._defineOutputs(outputSmallMoleculesBest=outputSetBest)
        self._defineSourceRelation(self.inputGrid, outputSetBest)
        self._defineSourceRelation(self.inputLibrary, outputSetBest)

    def _citations(self):
        return ['Morris2009']
