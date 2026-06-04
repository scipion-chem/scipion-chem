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
    """
    This protocol is used to extract ligands from an atomic structure (PDB, CIF, or ENT)
    and generate a SetOfSmallMolecules compatible with downstream Scipion workflows.

    The protocol identifies heteroatoms (HETATM records) in a protein–ligand complex,
    filters them based on size criteria, and exports each ligand as an independent
    small molecule structure. Optionally, the receptor structure can be cleaned by
    removing waters and selected chains before ligand extraction.

    Core Concepts
    -------------
    Atomic Structure:
        A macromolecular structure file (PDB, CIF, or similar) containing protein
        chains, water molecules, ions, and potential bound ligands.

    Ligand Identification:
        Ligands are detected as hetero-residues (HETATM) that are not water and
        contain a minimum number of heavy atoms.

    Heavy Atoms:
        Non-hydrogen atoms used to filter out small artifacts such as waters or ions.

    Residue Selection:
        A Biopython-based selection strategy used to extract individual ligand residues
        from a full structure.

    Structure Cleaning (optional):
        Pre-processing step to remove waters and/or selected chains before ligand
        extraction.

    Workflow
    --------
    1. Input atomic structure (PDB / CIF / ENT format).
    2. Optionally clean structure:
       - Remove water molecules
       - Remove selected chains
    3. Parse structure using Biopython PDB/MMCIF parser.
    4. Identify hetero-residues (ligands).
    5. Filter ligands by minimum number of heavy atoms.
    6. Extract each ligand as an independent structure file.
    7. Save receptor (cleaned structure) separately.
    8. Build SetOfSmallMolecules object containing all extracted ligands.

    Ligand Filtering Criteria
    -------------------------
    - Must be a hetero-residue (HETATM).
    - Must not be water molecules.
    - Must contain at least `nAtoms` heavy atoms (default: 4).

    Structure Cleaning Options
    --------------------------
    Clean Structure:
        Enables preprocessing of the input structure before ligand extraction.

    Remove Waters:
        Removes crystallographic water molecules from the structure.

    Remove Chains:
        Allows selection of specific protein chains to keep in the receptor model.

    Ligand Representation
    ---------------------
    Each ligand is stored as:
    - Independent PDB file
    - SmallMolecule object with:
      - molName (derived from residue name)
      - poseFile (ligand structure file)
      - gridId and confId (incremental identifiers)

    Output
    ------
    - outputSmallMolecules:
        SetOfSmallMolecules containing all extracted ligands.

    - Protein structure (cleaned):
        Receptor structure without ligands (and optionally without waters/chains).

    Use Cases
    ---------
    - Preparation of ligand libraries from co-crystal structures
    - Extraction of bound ligands for docking validation
    - Generation of datasets for binding site analysis
    - Separation of receptor and ligand for molecular modeling workflows
    """
    _label = 'Extract Ligand from structure'

    def _cleanStructureParams(self, form):
        clean = form.addGroup("Clean atomic structure")
        clean.addParam('cleanPDB', params.BooleanParam, default=False, label='Clean structure: ')
        clean.addParam("waters", params.BooleanParam, label='Remove waters: ', condition='cleanPDB', default=True,
                       help='Remove all waters molecules from the structure')

        clean.addParam("rchains", params.BooleanParam, label='Remove chains: ', default=False, condition='cleanPDB',
                       help='Remove specific chains in the proteins')

        clean.addParam("chain_name", params.StringParam, label="Chains to keep: ", condition="cleanPDB and rchains",
                       help="Select the chain(s) you want to keep in your structure")

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputStructure', params.PointerParam,
                      label="Input structure: ", allowsNull=False,
                      important=True, pointerClass='AtomStruct',
                      help='Atom structure from which you want to extract a ligand')
        form.addParam('nAtoms', params.IntParam,
                      label="Minimum number of heavy atoms: ", default=4,
                      help='Minimum number of heavy atoms for a ligand to be considered (used for ignoring '
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

        chainId = None
        if self.cleanPDB.get():
            if self.rchains.get():
                chain = json.loads(self.chain_name.get())  # From wizard dictionary
                chainId = chain["chain"].upper().strip()

            tmpCleanFile = self._getTmpPath('{}.cif'.format(getBaseName(inputStructureFile)))
            tmpCleanFile = cleanPDB(inputStructureFile, tmpCleanFile, self.waters.get(), False, chainId)

        else:
            tmpCleanFile = inputStructureFile

        ligandFiles = self.extract_ligands(tmpCleanFile)

        cleanFile = self._getPath('{}.pdb'.format(getBaseName(inputStructureFile)))

        # Keep only structure chains selected, without ligands for reference structure
        cleanedPDB = cleanPDB(inputStructureFile, cleanFile, self.waters.get(), True, chainId)

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
        for model in struct:
            for chain in model:
                for residue in chain:
                    heavyAtoms = [atom for atom in residue if atom.element != 'H']
                    if not isHet(residue) or len(heavyAtoms) < self.nAtoms.get():
                        continue
                    print(f"saving {chain} {residue}")
                    outFile = self._getPath(f"{struct_name}_{residue.get_resname()}.pdb")
                    ligandFiles.append(outFile)
                    io.save(outFile, ResidueSelect(chain, residue))

        return ligandFiles