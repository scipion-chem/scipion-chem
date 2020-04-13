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

from pyworkflow.protocol.params import (StringParam, PointerParam)
from pyworkflow.utils.path import createLink
from pwem.protocols import EMProtocol
from pwem.objects.data import Sequence
from bioinformatics.objects import SetOfDatabaseID

class ProtBioinformaticsRaptorX(EMProtocol):
    """This is a wrapper to http://raptorx.uchicago.edu/StructPredV2/predict/. The underlying program
       predicts the 3D structure of a protein from its aminoacid sequence"""
    _label = 'raptorX'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSeq', PointerParam, pointerClass="SetOfDatabaseID, Sequence",
                       label='Sequence:', allowsNull=False,
                       help="Input FASTA sequence, the set must have a single structure")
        form.addParam('title', StringParam,
                         label="Job title:", default="",
                         help="Put a label that helps you to identify the job")
        form.addParam('email', StringParam,
                         label="Email:", default="",
                         help="The web page will send an email to this address when the calculations are finished. "
                              "Copy the URL at the email in the Analyze Results of this protocol")

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

        title = self.title.get()
        if not title:
            title="ScipionRun%s"%os.path.split(self._getPath())[1][0:6]

        args='-F "jobname=%s" -F "email=%s" -F "seqeunces=%s" http://raptorx.uchicago.edu/StructPredV2/predict/' %\
             (title, self.email.get(), str(record.seq), self._getPath("output.html"))
        self.runJob("curl",args)


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

    def _citations(self):
        return ['Zheng2019'] # ***
