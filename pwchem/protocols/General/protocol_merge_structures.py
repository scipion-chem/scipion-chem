# **************************************************************************
# *
# * Authors:	Blanca Pueche (blanca.pueche@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307 USA
# *
# * All comments concerning this program package may be sent to the
# * e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

"""
"""

# Scipion em imports
from pwem.protocols import EMProtocol
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.objects import AtomStruct

from pathlib import Path

from Bio.PDB import (
    PDBParser, MMCIFParser,
    PDBIO, MMCIFIO,
    Structure, Model
)


class ProtMergeStructs(EMProtocol):
  """
    AI Generated:

    This protocol merges multiple atomic structures into a single macromolecular
    structure by combining all chains from the selected input structures into
    one model.

    The protocol preserves the coordinates and chain organization of each
    structure while automatically renaming chains when duplicate chain
    identifiers are encountered. The output is written in the same format
    (PDB or mmCIF) as the first input structure.

    Workflow
    --------
    1. Load all selected atomic structures.
    2. Parse each structure (PDB or mmCIF).
    3. Extract the first model from every structure.
    4. Copy all chains into a new merged structure.
    5. Automatically rename duplicated chain identifiers to ensure uniqueness.
    6. Save the merged structure.

    Input
    -----
    - inputStructs:
        List of atomic structures (AtomStruct objects) to merge.

        Supported formats:
        - PDB
        - mmCIF

    Chain Handling
    --------------
    Each chain from the input structures is copied into the merged model.

    If two or more structures contain chains with the same identifier, the
    protocol automatically assigns a new unique chain ID using the sequence:

        A-Z : a-z : 0-9

    This guarantees that every chain in the merged structure has a unique
    identifier.

    Output
    ------
    - outputStructure:
        Single merged AtomStruct containing all chains from the input
        structures.

        The output format matches that of the first input structure:
        - PDB if the first input is a PDB file.
        - mmCIF if the first input is a CIF/mmCIF file.

    Use Cases
    ---------
    - Building multimeric assemblies from individual chain predictions
    - Combining independently predicted protein chains
    - Preparing structures for docking or interface analysis
    - Creating complete complexes from separate structural models
    - Merging structures prior to visualization or downstream structural
      comparison

    Notes
    -----
    - Only the first model of each input structure is merged.
    - Atomic coordinates are preserved exactly as provided in the input
      structures.
    - No structural alignment or coordinate transformation is performed before
      merging. Input structures should therefore already be in the desired
      relative orientation.
  """
  _label = 'Merge structures'

  def _defineParams(self, form):
    form.addSection(label=Message.LABEL_INPUT)

    group = form.addGroup('Input')
    group.addParam('inputStructs', params.MultiPointerParam, label="Input structures: ",
                   pointerClass='AtomStruct',
                   help='Structures to merge into one single structure.')

  # --------------------------- Steps functions --------------------
  def _insertAllSteps(self):
    self._insertFunctionStep(self.createOutputStep)

  def createOutputStep(self):
      structure = Structure.Structure("merged")
      model = Model.Model(0)
      structure.add(model)

      used_chain_ids = set()
      chain_pool = iter(
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
          "abcdefghijklmnopqrstuvwxyz"
          "0123456789"
      )

      outputExt = None

      for ptr in self.inputStructs:
          atomStruct = ptr.get()
          structFile = atomStruct.getFileName()

          ext = Path(structFile).suffix.lower()

          if outputExt is None:
              outputExt = ext

          if ext == ".pdb":
              parser = PDBParser(QUIET=True)
          elif ext in [".cif", ".mmcif"]:
              parser = MMCIFParser(QUIET=True)
          else:
              raise Exception(f"Unsupported structure format: {structFile}")

          print(f"Reading {structFile}")

          s = parser.get_structure("tmp", structFile)

          try:
              modelIn = next(s.get_models())
          except StopIteration:
              raise Exception(f"No models found in structure: {structFile}")

          for chain in modelIn:
              newChain = chain.copy()

              if newChain.id in used_chain_ids:
                  while True:
                      newId = next(chain_pool)
                      if newId not in used_chain_ids:
                          newChain.id = newId
                          break

              used_chain_ids.add(newChain.id)
              model.add(newChain)

      # Write output in the same format as the first input
      if outputExt == ".pdb":
          outFile = self._getPath("merged_struct.pdb")
          io = PDBIO()
      else:
          outFile = self._getPath("merged_struct.cif")
          io = MMCIFIO()

      io.set_structure(structure)
      io.save(outFile)

      output = AtomStruct(filename=outFile)
      self._defineOutputs(outputStructure=output)


  # --------------------------- INFO functions -----------------------------------
