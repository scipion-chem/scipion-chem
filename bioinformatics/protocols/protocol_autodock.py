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

import glob
import os

from pwem.protocols import EMProtocol
from pyworkflow.protocol.params import PointerParam, IntParam, FloatParam
from bioinformatics import Plugin as bioinformatics_plugin
from pyworkflow.utils.path import makePath, createLink

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
        # self._insertFunctionStep('createOutputStep', prerequisites=dockSteps)

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

        self.runJob(bioinformatics_plugin.getMGLPath('bin/pythonsh'),
                    bioinformatics_plugin.getADTPath('Utilities24/prepare_dpf42.py')+args)

        fnSmallLocal = os.path.split(fnSmall)[1]
        createLink(fnSmall,os.path.join(fnSmallDir,fnSmallLocal))
        createLink(fnReceptor,os.path.join(fnSmallDir,"atomStruct.pdbqt"))
#        createLink(os.path.join(fnGridDir,"atomStruct.maps.fld"),os.path.join(fnSmallDir,"atomStruct.maps.fld"))
#        createLink(os.path.join(fnGridDir,"atomStruct.maps.xyz"),os.path.join(fnSmallDir,"atomStruct.maps.xyz"))
#        for fn in glob.glob(os.path.join(fnGridDir,"*map")):
#            createLink(fn, os.path.join(fnSmallDir, os.path.split(fn)[1]))

        args = " -r atomStruct.pdbqt -l %s -o library.gpf"%fnSmallLocal
        self.runJob(bioinformatics_plugin.getMGLPath('bin/pythonsh'),
                    bioinformatics_plugin.getADTPath('Utilities24/prepare_gpf4.py') + args,
                    cwd=fnSmallDir)

        args = "-p library.gpf -l library.glg"
        self.runJob(bioinformatics_plugin.getAutodockPath("autogrid4"), args, cwd=fnSmallDir)

        args = "-p %s.dpf -l %s.dlg"%(fnBase,fnBase)
        self.runJob(bioinformatics_plugin.getAutodockPath("autodock4"), args, cwd=fnSmallDir)

    def _citations(self):
        return ['Morris2009']
