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
from pwchem.objects import SmallMolecule, SetOfSmallMolecules

from pathlib import Path
import os

from Bio.PDB import (
    PDBParser, MMCIFParser,
    PDBIO, MMCIFIO,
    Structure, Model
)

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import OPENBABEL_DIC


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

    group.addParam(
        'inputLigands',
        params.MultiPointerParam,
        label="Input ligands: ",
        pointerClass='SetOfSmallMolecules',
        allowsNull=True,
        help='Optional sets of ligands to add to the merged structure.'
    )

  # --------------------------- Steps functions --------------------
  def _insertAllSteps(self):
    self._insertFunctionStep(self.createOutputStep)

  def createOutputStep(self):
      structure = Structure.Structure("merged")
      model = Model.Model(0)
      structure.add(model)

      usedChainIds = set()
      chain_pool = iter(
          "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
          "abcdefghijklmnopqrstuvwxyz"
          "0123456789"
      )

      def getUniqueChainId(preferredId=None):
          """
          Return a chain ID that is not already present in the merged model.
          """
          if preferredId and preferredId not in usedChainIds:
              chainId = preferredId
          else:
              while True:
                  try:
                      chainId = next(chain_pool)
                  except StopIteration:
                      raise RuntimeError(
                          "No more unique chain IDs available."
                      )

                  if chainId not in usedChainIds:
                      break

          usedChainIds.add(chainId)
          return chainId

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
              raise Exception(
                  f"Unsupported structure format: {structFile}"
              )

          print(f"[ProtMergeStructs] Structure: {structFile}")

          s = parser.get_structure("tmp", structFile)

          try:
              modelIn = next(s.get_models())
          except StopIteration:
              raise Exception(
                  f"No models found in structure: {structFile}"
              )

          for chain in modelIn:
              newChain = chain.copy()
              newChain.id = getUniqueChainId(chain.id)
              model.add(newChain)

      if self.inputLigands:
          for ptr in self.inputLigands:

              ligandSet = ptr.get()

              if ligandSet is None:
                  continue

              for ligand in ligandSet:
                  if ligandSet.isDocked():
                      ligandFile = ligand.getPoseFile()
                  else:
                      ligandFile = ligand.getFileName()

                  if not ligandFile:
                      continue

                  ligandFile = os.path.abspath(ligandFile)

                  ext = Path(ligandFile).suffix.lower()

                  if ext in [".sdf", ".mol2"]:
                      convDir = self._getExtraPath("ligands")
                      os.makedirs(convDir, exist_ok=True)

                      ligandName = Path(ligandFile).stem

                      args = (
                          f' -i "{ligandFile}"'
                          f' -of pdb'
                          f' -o "{ligandName}"'
                          f' -od "{os.path.abspath(convDir)}"'
                      )

                      print(
                          f"[ProtMergeStructs] Converting ligand: "
                          f"{ligandFile}"
                      )

                      pwchemPlugin.runScript(
                          self,
                          'obabel_IO.py',
                          args,
                          env=OPENBABEL_DIC,
                          cwd=convDir
                      )

                      ligandFile = os.path.join(
                          convDir,
                          f"{ligandName}.pdb"
                      )
                      ligandFile = self._fixLigandPDBAtomNames(
                          ligandFile
                      )

                  ext = Path(ligandFile).suffix.lower()

                  if ext == ".pdb":
                      parser = PDBParser(QUIET=True)

                  elif ext in [".cif", ".mmcif"]:
                      parser = MMCIFParser(QUIET=True)

                  else:
                      raise ValueError(
                          f"Unsupported ligand format: {ligandFile}"
                      )

                  print(
                      f"[ProtMergeStructs] Adding ligand: "
                      f"{ligandFile}"
                  )

                  s = parser.get_structure(
                      "ligand",
                      ligandFile
                  )

                  try:
                      modelIn = next(s.get_models())

                  except StopIteration as e:
                      raise ValueError(
                          f"No models found in ligand: {ligandFile}"
                      ) from e

                  for chain in modelIn:

                      newChain = chain.copy()
                      newChain.id = getUniqueChainId()

                      for residue in newChain:
                          _, resseq, icode = residue.id

                          residue.resname = "LIG"

                          residue.id = (
                              "H_LIG",
                              resseq,
                              icode
                          )

                      model.add(newChain)

                      print(
                          f"Added ligand chain {newChain.id}"
                      )

      if outputExt == ".pdb":
          outFile = self._getPath(
              "merged_struct.pdb"
          )
          io = PDBIO()
      else:
          outFile = self._getPath(
              "merged_struct.cif"
          )
          io = MMCIFIO()

      io.set_structure(structure)
      io.save(outFile)

      print(
          f"Merged structure written to: "
          f"{outFile}"
      )

      output = AtomStruct(
          filename=outFile
      )
      self._defineOutputs(
          outputStructure=output
      )


  # --------------------------- INFO functions -----------------------------------

  def _fixLigandPDBAtomNames(self, pdbFile):
      fixedFile = self._getTmpPath(
          f"{Path(pdbFile).stem}_fixed.pdb"
      )

      atomCounters = {}

      with open(pdbFile) as f_in, open(fixedFile, "w") as f_out:
          for line in f_in:
              if line.startswith(("ATOM  ", "HETATM")):
                  element = line[76:78].strip()

                  if not element:
                      element = line[12:16].strip()[0]

                  atomCounters[element] = atomCounters.get(element, 0) + 1
                  atomName = f"{element}{atomCounters[element]}"

                  atomName = atomName[:4]

                  line = line[:12] + f"{atomName:>4}" + line[16:]

              f_out.write(line)

      return fixedFile

