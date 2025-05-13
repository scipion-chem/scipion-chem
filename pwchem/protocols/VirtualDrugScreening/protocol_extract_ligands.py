# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Aida Pinacho Pérez
# *              Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

"""
This module will extract the ligand from a complex pdb.
"""

import os, sys, json
from Bio.PDB import PDBParser, PDBIO, Select, MMCIFParser

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import *


class ResidueSelect(Select):
    def __init__(self, chain, residue):
        self.chain = chain
        self.residue = residue

    def accept_chain(self, chain):
        return chain.id == self.chain.id

    def accept_residue(self, residue):
        """ Recognition of heteroatoms - Remove water molecules """
        return residue == self.residue and isHet(residue)


class ProtExtractLigands(EMProtocol):
    _label = 'Extract Ligand from structure'

    def __init__(self, **kwargs):
        EMProtocol.__init__(self, **kwargs)

    def _cleanStructureParams(self, form):
        clean = form.addGroup("Clean atomic structure")
        clean.addParam('cleanPDB', params.BooleanParam, default=False, label='Clean PDB: ')
        clean.addParam("waters", params.BooleanParam,
                       label='Remove waters', condition='cleanPDB',
                       default=True, important=True,
                       help='Remove all waters molecules from a pdb file')

        clean.addParam("rchains", params.BooleanParam,
                       label='Remove redundant chains',
                       default=False, important=True, condition='cleanPDB',
                       help='Remove redundant chains in the proteins')

        clean.addParam("chain_name", params.StringParam,
                       label="Conserved chain",
                       important=True,
                       condition="cleanPDB and rchains==True",
                       help="Select the chain on which you want to carry out the "
                            "molecular docking. You must know the protein and structure "
                            "file that you loaded. \n\nFor example, the protein mdm2 "
                            "(4ERF) has a C1 symmetry, which indicates that its chains "
                            "are at least 95% equal, so you would write A, C or E "
                            "(names of the chains).")

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure from which you want to extract a ligand')
        form.addParam('nAtoms', params.IntParam,
                      label="Minimum number of atoms: ", default=4,
                      help='Minimum number of atoms for a ligand to be considered (used for ignoring '
                           'waters for example)')

        self._cleanStructureParams(form)

    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('createOutputStep')


    def createOutputStep(self):
        inputStructureFile = os.path.abspath(self.inputStructure.get().getFileName())
        inputBaseName, inputExt = inputStructureFile.split('.')
        if not inputExt in ['pdb', 'ent', 'cif']:
            print('Unknown format for file ', inputStructureFile)
            exit()

        if self.rchains.get():
            chain = json.loads(self.chain_name.get())  # From wizard dictionary
            chain_id = chain["chain"].upper().strip()
        else:
            chain_id = None

        cleanFile = self._getPath('{}.pdb'.format(getBaseName(inputStructureFile)))
        tmpCleanFile = self._getTmpPath('{}.pdb'.format(getBaseName(inputStructureFile)))

        # Keep only structure chains selected, with ligands for extraction
        tmpCleanFile = cleanPDB(inputStructureFile, tmpCleanFile, self.waters.get(), False, chain_id)
        ligandFiles = self.extract_ligands(tmpCleanFile)

        # Keep only structure chains selected, without ligands for reference structure
        cleanedPDB = cleanPDB(inputStructureFile, cleanFile, self.waters.get(), True, chain_id)

        outputSet = SetOfSmallMolecules().create(outputPath=self._getPath())
        outputSet.setProteinFile(cleanedPDB)
        outputSet.setDocked(True)
        for i, lFile in enumerate(ligandFiles):
            molName = getBaseName(lFile).split('-')[0]
            oMol = SmallMolecule(smallMolFilename=lFile, molName=molName, poseFile=lFile, poseId=1)
            oMol.setGridId(i+1)
            oMol.setConfId(i+1)
            outputSet.append(oMol)
        self._defineOutputs(outputSmallMolecules=outputSet)

    # --------------------------- INFO functions -----------------------------------

    def extract_ligands(self, struct_file):
        """ Extraction of the heteroatoms of .pdb files """
        struct_name = getBaseName(struct_file)
        if struct_file.endswith('.pdb') or struct_file.endswith('.ent'):
            struct = PDBParser().get_structure(struct_name, struct_file)
        elif struct_file.endswith('.cif'):
            struct = MMCIFParser().get_structure(struct_name, struct_file)
        else:
            print('Unknown format for file ', struct_file)
            exit()

        ligandFiles = []
        io = PDBIO()
        io.set_structure(struct)
        i = 1
        for model in struct:
            for chain in model:
                for residue in chain:
                    if not isHet(residue) or len(list(residue.get_atoms())) < self.nAtoms.get():
                        continue
                    print(f"saving {chain} {residue}")
                    outFile = self._getPath(f"{struct_name}_{residue.get_resname()}-{i}.pdb")
                    ligandFiles.append(outFile)
                    io.save(outFile, ResidueSelect(chain, residue))
                    i += 1

        return ligandFiles