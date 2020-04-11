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
import sys

from Bio import SeqIO

from pyworkflow.protocol.params import (EnumParam, PointerParam, StringParam)
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from pwem.objects.data import Sequence
from pwem.convert.atom_struct import AtomicStructHandler
from bioinformatics.objects import SetOfDatabaseID

class ProtBioinformaticsCQuark(EMProtocol):
    """Query CQuark server (https://zhanglab.ccmb.med.umich.edu/C-QUARK/) with an aminoacid sequence.
    The server returns several proposals of 3D structure. """
    _label = 'cquark'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSeq', PointerParam, pointerClass="SetOfDatabaseID, Sequence",
                       label='Sequence:', allowsNull=False,
                       help="Input FASTA sequence, the set must have a single structure")

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        outFileName = self._getExtraPath("sequence.fasta")
        if isinstance(self.inputSeq.get(),Sequence):
            fh=open(outFileName,"w")
            fh.write(">%s\n"%self.inputSeq.get().getSeqName())
            fh.write("%s\n"%self.inputSeq.get().getSequence())
            fh.close()
        elif isinstance(self.inputSeq.get(),SetOfDatabaseID):
            obj=self.inputSeq.get().getFirstItem()
            createLink(obj._uniprotFile.get(),outFileName)
            fh=open(outFileName)
        record = SeqIO.read(outFileName, "fasta")

        summary = "C-Quark does not allow programmatic access. Go to: https://zhanglab.ccmb.med.umich.edu/C-QUARK/\n"
        summary +="You will have to register and submit the sequence yourself. Use the job ID 'ScipionRun %s'\n"% \
                  os.path.split(self._getPath())[1][0:6]
        summary +="When you get an answer in your email, open the Analyze Results before 3 days and provide the results URL\n\n\n"
        summary +=str(record.seq)+"\n"
        fh=open(self._getPath("summary.txt"),'w')
        fh.write(summary)

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = []
        if isinstance(self.inputSeq.get(),SetOfDatabaseID):
            if len(self.inputSeq.get())!=1:
                errors.append("The input list can only have a single sequence")
            obj=self.inputSeq.get().getFirstItem()
            if not hasattr(obj,"_uniprotFile"):
                errors.append("The input list does not have a sequence file")
        from pyworkflow.utils.which import which
        if which("wget") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install wget, yum install wget or equivalent in your system")
        return errors

    def _summary(self):
        summary = []
        fnSummary = self._getPath("summary.txt")
        if os.path.exists(fnSummary):
            fh=open(fnSummary)
            for line in fh.readlines():
                summary.append(line.strip())
        return summary

    def _citations(self):
        return ['Zheng2019']
