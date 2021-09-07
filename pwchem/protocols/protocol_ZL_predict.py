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

from pyworkflow.protocol.params import (EnumParam, PointerParam)
from pwem.protocols import EMProtocol
from pwchem.utils.utilsFasta import *

class ProtChemZLPredict(EMProtocol):
    """Query Zhang-Lab servers (https://zhanglab.ccmb.med.umich.edu/) with an aminoacid sequence.
    The server returns several proposals of 3D structure. """
    _label = 'ZL predict'
    _program = ""

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputSeq', PointerParam, pointerClass="SetOfDatabaseID, Sequence",
                       label='Sequence:', allowsNull=False,
                       help="Input FASTA sequence, the set must have a single structure")
        form.addParam('method', EnumParam, choices=['C-Quark [<500] aa', 'I-Tasser [10,1500] aa',
                                                    'Quark [<200] aa', 'Lomets [<1500] aa', 'CEthreader',
                                                    'Muster', 'Segmer'],
                      label='Method', default=0,
                      help='Go to https://zhanglab.ccmb.med.umich.edu for an explanation of each method')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('searchStep')

    def searchStep(self):
        seq = copyFastaSequenceAndRead(self)

        url={}
        url[0]='https://zhanglab.ccmb.med.umich.edu/C-QUARK/'
        url[1]='https://zhanglab.ccmb.med.umich.edu/I-TASSER/'
        url[2]='https://zhanglab.ccmb.med.umich.edu/QUARK/'
        url[3]='https://zhanglab.ccmb.med.umich.edu/LOMETS/'
        url[4]='https://zhanglab.ccmb.med.umich.edu/CEthreader/'
        url[5]='https://zhanglab.ccmb.med.umich.edu/MUSTER/'
        url[6]='https://zhanglab.ccmb.med.umich.edu/SEGMER/'

        summary = "Zhang Lab does not allow programmatic access. Go to: %s\n"%url[self.method.get()]
        summary +="You may have to register depending on the method and submit the sequence yourself. Use the job ID 'ScipionRun%s'\n"% \
                  os.path.split(self._getPath())[1][0:6]
        summary +="When you get an answer in your email, open the Analyze Results before 3 days and provide the results URL\n\n\n"

        ok=True
        if self.method.get()==0 and len(seq)>500:
            ok=False
        elif self.method.get()==1 and (len(seq)<10 or len(seq)>1500):
            ok=False
        elif self.method.get() == 2 and len(seq) > 200:
            ok=False
        elif self.method.get()==3 and len(seq)>1500:
            ok=False
        if ok:
            summary +=seq+"\n"
        else:
            summary +="The input sequence does not fit into the method requirements"
        fh=open(self._getPath("summary.txt"),'w')
        fh.write(summary)

    # --------------------------- UTILS functions ------------------
    def _validate(self):
        errors = checkInputHasFasta(self)
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
