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

from pyworkflow.protocol.params import PointerParam, EnumParam, StringParam, BooleanParam

from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import AtomicStructHandler
from pwchem import Plugin as pwchem_plugin
from pwem.objects.data import AtomStruct

class ProtBioinformaticsADTPrepare(EMProtocol):
    def _defineParamsBasic(self, form):
        choicesRepair = ['None', 'Bonds hydrogens', 'Bonds', 'Hydrogens']
        if self.typeRL=="target":
            choicesRepair.append('Check hydrogens')
        form.addParam('repair', EnumParam, choices=choicesRepair,
                      default=0, label='Repair action:',
                      help='Bonds hydrogens: build bonds and add hydrogens\n'
                           'Bonds: build a single bond from each atom with no bonds to its closest neighbor\n'
                           'Hydrogens: add hydrogens\n'
                           'Check hydrogens: add hydrogens only if there are none already')
        form.addParam('preserveCharges', EnumParam, choices=['Add gasteiger charges','Preserve input charges',
                                                             'Preserve charges of specific atoms'],
                      default=0, label='Charge handling')
        form.addParam('chargeAtoms', StringParam, default="", condition='preserveCharges==2',
                      label='Atoms to preserve charge', help='Separated by commas: Zn, Fe, ...')
        form.addParam('nphs', BooleanParam, default=True,
                      label='Merge charges and remove non-polar hydrogens')
        form.addParam('lps', BooleanParam, default=True,
                      label='Merge charges and remove lone pairs')
        form.addParam('waters', BooleanParam, default=True,
                      label='Remove water residues')
        if self.typeRL=="target":
            form.addParam('nonstdres', BooleanParam, default=True,
                          label='Remove chains composed entirely of non-standard residues')
            form.addParam('nonstd', BooleanParam, default=False,
                          label='Remove non-standard residues from all chains')

    # --------------------------- INSERT steps functions --------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def callPrepare(self, prog, args):
        if self.repair.get()==1:
            args+=' -A bonds_hydrogens'
        elif self.repair.get()==2:
            args += ' -A bonds'
        elif self.repair.get()==3:
            args += ' -A hydrogens'
        elif self.repair.get()==4:
            args += ' -A checkhydrogens'

        if self.preserveCharges.get()==1:
            args+=" -C"
        elif self.preserveCharges.get()==2:
            for atom in self.chargeAtoms.get().split(','):
                args+=" -p %s"%atom.strip()

        cleanup=" -U "
        first=True
        if self.nphs.get():
            cleanup+=" nphs"
            first=False
        if self.lps.get():
            if not first:
                cleanup+="_"
            cleanup+="lps"
            first=False
        if self.waters.get():
            if not first:
                cleanup+="_"
            cleanup+="waters"
            first=False
        if self.typeRL=="target":
            if self.nonstdres.get():
                if not first:
                    cleanup+="_"
                cleanup+="nonstdres"
        if cleanup!="-U ":
            args+=cleanup

        if self.typeRL=="target":
            if self.nonstd.get():
                args+=" -e"

        self.runJob(pwchem_plugin.getMGLPath('bin/pythonsh'),
                    pwchem_plugin.getADTPath('Utilities24/%s.py'%prog)+args)

    def createOutput(self):
        fnOut = self._getExtraPath('atomStruct.pdbqt')
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputStructure, target)

class ProtBioinformaticsADTPrepareReceptor(ProtBioinformaticsADTPrepare):
    """Prepare receptor using Autodocking Tools from MGL"""
    _label = 'target preparation ADT'
    _program = ""

    def _defineParams(self, form):
        self.typeRL="target"
        form.addSection(label='Input')
        form.addParam('inputStructure', PointerParam, pointerClass="AtomStruct",
                      label='Atomic Structure:', allowsNull=False,
                      help='It must be in pdb,mol2,pdbq,pdbqs,pdbqt format, you may use Schrodinger convert to change it')
        ProtBioinformaticsADTPrepare._defineParamsBasic(self, form)

    def preparationStep(self):
        if self.inputStructure.get().getFileName().endswith('.cif'):
            fnIn = self._getTmpPath("atomStructIn.pdb")
            aStruct1 = AtomicStructHandler(self.inputStructure.get().getFileName())
            aStruct1.write(fnIn)
        else:
            fnIn = self.inputStructure.get().getFileName()
        fnOut = self._getExtraPath('atomStruct.pdbqt')

        args = ' -v -r %s -o %s'%(fnIn,fnOut)
        ProtBioinformaticsADTPrepare.callPrepare(self,"prepare_receptor4",args)

    def createOutput(self):
        fnOut = self._getExtraPath('atomStruct.pdbqt')
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputStructure, target)

    def _validate(self):
        errors = []
        if not self.inputStructure.get().getFileName().endswith('.mol2') and \
           not self.inputStructure.get().getFileName().endswith('.pdb') and \
           not self.inputStructure.get().getFileName().endswith('.pdbq') and \
           not self.inputStructure.get().getFileName().endswith('.pdbqt') and \
           not self.inputStructure.get().getFileName().endswith('.pdbqs') and \
           not self.inputStructure.get().getFileName().endswith('.cif'):
            errors.append('The input structure must be either .mol2, .pdb, .pdbq, .pdbqt, .pdbqs or .cif')
            errors.append("Current name: %s"%self.inputStructure.get().getFileName())
        return errors
