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
        form.addParam('title', StringParam,
                         label="Job title:", default="Scipion search",
                         help="Put a label that helps you to identify the job")
        form.addParam('email', StringParam,
                         label="Email:", default="",
                         help="The web page will send an email to this address when the calculations are finished. "
                              "Copy the URL at the email in the Analyze Results of this protocol")
        form.addParam('password', StringParam,
                         label="Password in C-Quark:", default="",
                         help="The server asks for a password associated to this email")

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

        args ='-d SEQUENCE="%s"'%record.seq
        args+=' -d REPLY-E-MAIL="%s"'%self.email.get()
        args+=' -d PASSWORD="%s"'%self.password.get()
        args+=' -d TARGET="%s"'%self.title.get()
        args+=" https://zhanglab.ccmb.med.umich.edu/C-QUARK/"
        args+=" > %s"%self._getExtraPath("submission_response.html")
        print(args)
        #self.runJob("curl",args)

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = []
        if self.email.get()=="":
            errors.append("The email cannot be empty")
        if self.password.get()=="":
            errors.append("Password cannot be empty")
        if isinstance(self.inputSeq.get(),SetOfDatabaseID):
            if len(self.inputSeq.get())!=1:
                errors.append("The input list can only have a single sequence")
            obj=self.inputSeq.get().getFirstItem()
            if not hasattr(obj,"_uniprotFile"):
                errors.append("The input list does not have a sequence file")
        from pyworkflow.utils.which import which
        if which("curl") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install curl, yum install curl or equivalent in your system")
        if which("wget") == "":
            errors.append("Cannot find curl in the path. Install with apt-get install wget, yum install wget or equivalent in your system")
        return errors

    def _summary(self):
        summary = []
        return summary

    def _citations(self):
        return ['Zheng2019']
