# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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
import os, json

from pyworkflow.protocol import params
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol

from pwchem.utils import cleanPDB, getBaseName, removeStartWithLines
from pwchem import Plugin as pwchemPlugin

class ProtChemPrepareReceptor(EMProtocol):
    """Prepare receptor by removing HETATM atoms, waters, keeping only specific chains..."""
    _label = 'target preparation'
    _program = ""
    def defineCleanParams(self, form, w=True, h=True, hk=True, c=True):
      cGroup = form.addGroup('Clean Structure File')
      if w:
          cGroup.addParam("waters", params.BooleanParam,
                         label='Remove waters: ', default=True, important=True,
                         help='Remove all waters molecules from a pdb file')
      if h:
          cGroup.addParam("HETATM", params.BooleanParam,
                         label='Remove ligands HETATM: ', default=True, important=True,
                         help='Remove all ligands and HETATM contained in the protein')
      if hk:
          cGroup.addParam('het2keep', params.StringParam, condition='HETATM', default='',
                         label='Heteroatom residues to keep: ', expertLevel=params.LEVEL_ADVANCED,
                         help='Specify specific heteroatom residues to keep')
      if c:
          cGroup.addParam("rchains", params.BooleanParam,
                         label='Select specific chains: ', default=False, important=True,
                         help='Keep only the chains selected')

          cGroup.addParam("chain_name", params.StringParam,
                         label="Keep chains: ", important=True, condition="rchains==True",
                         help="Select the chain(s) you want to keep in the structure")
      return cGroup

    def _defineParams(self, form):
        form.addSection(label='Input')
        form.addParam('inputAtomStruct', params.PointerParam, pointerClass="AtomStruct",
                      label='Atomic Structure: ', allowsNull=False,
                      help='It must be in pdb,mol2,pdbq,pdbqs,pdbqt format, you may use Schrodinger convert to change it')
        self.defineCleanParams(form)

        pGroup = form.addGroup('PDBFixer')
        pGroup.addParam("usePDBFixer", params.BooleanParam, label='Use PDBFixer: ', default=False, important=True,
                        help='Whether to use PDBFixer to further.')
        pGroup.addParam('addAtoms', params.EnumParam, default=0, condition="usePDBFixer",
                        label="Add missing atoms: ", choices=['All', 'Heavy', 'Hydrogen', 'None'],
                        help='Use PDBFixer to add the missing atoms specified in the PDB atomic structure')
        pGroup.addParam('addRes', params.BooleanParam, default=True,
                        label="Add missing residues: ",  condition="usePDBFixer",
                        help='Use PDBFixer to add missing residues')
        pGroup.addParam('repNonStd', params.BooleanParam, default=False,
                        label="Replace non-standard residues: ", condition="usePDBFixer",
                        help='Use PDBFixer to replace nonstandard residues with standard equivalents')

        pGroup.addParam('extraClean', params.BooleanParam, default=False, condition="usePDBFixer",
                        label='Perform extra cleaning after PDBFixer: ', expertLevel=params.LEVEL_ADVANCED,
                        help='Perform extra cleaning round after PDBFixer has been used. Sometimes, PDBFixer will '
                             'point out some HETATMS that were originally labeled as ATOM.')

    def _insertAllSteps(self):
        self._insertFunctionStep('preparationStep')
        self._insertFunctionStep('createOutput')

    def preparationStep(self):
        chain_ids = None
        if self.rchains.get():
            chainJson = json.loads(self.chain_name.get())  # From wizard dictionary
            if 'chain' in chainJson:
                chain_ids = [chainJson["chain"].upper().strip()]
            elif 'model-chain' in chainJson:
                modelChains = chainJson["model-chain"].upper().strip()
                chain_ids = [x.split('-')[1] for x in modelChains.split(',')]

        het2keep = self.het2keep.get().split(', ')
        inFile = self.inputAtomStruct.get().getFileName()
        cleanFile = self.getCleanedFile()
        cleanFile = cleanPDB(inFile, cleanFile, self.waters.get(), self.HETATM.get(), chain_ids, het2keep)

        if self.usePDBFixer.get():
          addResStr = ' --add-residues' if self.addRes else ''
          repNStdStr = ' --replace-nonstandard' if self.repNonStd else ''
          addAtomsStr = self.getEnumText("addAtoms").lower()

          prepFile = self.getPreparedFile()
          args = f'{cleanFile} --add-atoms={addAtomsStr}{addResStr}{repNStdStr} --output {prepFile}'
          pwchemPlugin.runOPENBABEL(self, 'pdbfixer', args=args, cwd=self._getExtraPath())

          if self.extraClean.get():
            removeStartWithLines(prepFile, 'HETATM')

    def createOutput(self):
        fnOut = self.getPreparedFile() if self.usePDBFixer.get() else self.getCleanedFile()
        if os.path.exists(fnOut):
            target = AtomStruct(filename=fnOut)
            self._defineOutputs(outputStructure=target)
            self._defineSourceRelation(self.inputAtomStruct, target)

    def getPreparedFile(self):
        inFile = self.inputAtomStruct.get().getFileName()
        return os.path.abspath(self._getPath(f'{getBaseName(inFile)}.pdb'))

    def getCleanedFile(self):
        inFile = self.inputAtomStruct.get().getFileName()
        ext = os.path.splitext(inFile)[1]
        return os.path.abspath(self._getPath(f'{getBaseName(inFile)}{ext}'))
