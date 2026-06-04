# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

from pyworkflow.protocol import params
import pyworkflow.object as pwobj

from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules, SmallMoleculesLibrary
from pwchem.utils import getFilteredOutput


class ProtExtractInteractingMols(EMProtocol):
  """
  This protocol is used to define structural regions (ROIs) based on ligand–receptor
  contacts extracted from docked small molecules.

  The protocol identifies residues or atoms in the receptor that are within a given
  distance threshold of ligand atoms and groups them into structurally coherent regions
  of interest (ROIs). These ROIs can be clustered spatially and optionally mapped to
  surface coordinates to better represent binding pockets.

  Core Concepts
  -------------
  Small Molecules:
      Docked ligands provided as a SetOfSmallMolecules, each associated with
      a receptor structure.

  Contact Definition:
      A ligand–receptor contact is defined when any ligand atom is within a
      user-specified distance threshold (Å) from a receptor atom.

  Contact Representation:
      Contacts can be stored at either atom-level or residue-level resolution
      for both ligand and receptor.

  ROI (Region of Interest):
      A spatial cluster of receptor contact points representing a putative
      binding pocket or interaction region.

  Surface Mapping:
      Optionally maps contact coordinates to the nearest solvent-accessible
      surface points using MSMS-generated surface data.

  Workflow
  --------
  1. Input a SetOfSmallMolecules with associated receptor structure.
  2. Convert ligand files into a unified format (PDB or MAE conversion if needed).
  3. Parse receptor and ligand atomic structures.
  4. Compute pairwise distances between ligand and receptor atoms.
  5. Identify atom–atom contacts using a distance threshold.
  6. Store contacts in structured text files (per ligand).
  7. Aggregate contacts across ligands and extract unique receptor contact points.
  8. Optionally filter contacts by chain selection.
  9. Convert contact points into 3D coordinates.
  10. Optionally map coordinates to molecular surface points.
  11. Cluster coordinates into ROIs using hierarchical spatial clustering.
  12. Generate structural ROI objects and store results.

  Contact Analysis
  ----------------
  - Distance threshold (threshold):
      Maximum distance between ligand and receptor atoms to consider a contact.

  - Ligand contact level:
      Defines whether ligand contacts are stored at atom or residue level.

  - Receptor contact level:
      Defines whether receptor contacts are stored at atom or residue level.

  - Contact dictionary:
      Maps ligand atoms → list of receptor atoms in contact.

  Clustering Logic
  ----------------
  - Spatial clustering is performed using pairwise Euclidean distances.
  - Clusters are defined using a maximum intra-cluster distance (maxIntraDistance).
  - Optionally, surface mapping replaces raw coordinates with nearest surface points.
  - Clusters can be further refined depending on ligand grouping or selection mode.

  Surface Mapping
  ---------------
  When enabled:
  - Contact coordinates are projected onto solvent-accessible surface points.
  - A maximum depth threshold (maxDepth) defines allowable mapping distance.
  - MSMS is used to compute molecular surface representation.

  Output
  ------
  - outputStructROIs:
      SetOfStructROIs containing clustered structural regions.

  - StructROI objects:
      Each ROI includes:
          - Pocket coordinates
          - Contact annotations
          - Optional surface-mapped representation
          - Optional PDB/hetatm export

  Use Cases
  ---------
  - Identification of ligand binding pockets
  - Extraction of interaction hotspots from docking experiments
  - Structural comparison of binding regions across ligands
  - Preparation of regions for downstream docking or drug design workflows
  """
  _label = 'extract interacting molecules'

  def __init__(self, **kwargs):
    EMProtocol.__init__(self, **kwargs)
    self.stepsExecutionMode = params.STEPS_PARALLEL

  def _defineParams(self, form):
    form.addSection(label='Input')
    iGroup = form.addGroup('Input')
    iGroup.addParam('inputSequences', params.PointerParam, pointerClass="SetOfSequencesChem",
                    label='Input protein sequences: ',
                    help="Set of protein sequences containing interacting information (from ConPLex from example)")

    mGroup = form.addGroup('Filters')
    mGroup.addParam('chooseSeq', params.StringParam, label='Extract interacting molecules for protein: ', default='All',
                    help='Extract the interacting molecules for the selected protein sequence. '
                         'Use the wizard to get a list of the protein sequences in the input')
    mGroup.addParam('chooseMol', params.StringParam, label='Extract only from this subset of molecules: ',
                    default='All', expertLevel=params.LEVEL_ADVANCED,
                    help='Extract the only interacting molecules for the selected protein sequence in this subset. '
                         'Use the wizard to get a list of the interacting molecules in the input')
    mGroup.addParam('chooseScore', params.StringParam, label='Filter by this score: ', default='',
                    help='Choose to filter by the score provided by a specific tool.')
    mGroup.addParam('scThres', params.FloatParam, label='Score threshold: ', default=0.3,
                    help='Score threshold to filter interacting molecules below it for the selected protein sequences')

  def _insertAllSteps(self):
    self._insertFunctionStep(self.createOutputStep)

  def defineMolsOutput(self, intDic, molNames, seqName):
    outMols = SetOfSmallMolecules().create(outputPath=self._getPath())
    scoreType = self.chooseScore.get()
    score = f'score_{scoreType}'

    for mol in self.getInputMols():
      molName = mol.getMolName()
      if molName in molNames:
        if seqName:
          mol._interactScore = pwobj.Float(intDic[seqName][molName][score])
        outMols.append(mol)

    self._defineOutputs(outputSmallMolecules=outMols)

  def defineLibraryOutput(self, intMols, intDic, molNames, seqName):
    inFile, oFile = intMols.getFileName(), self._getPath('outputLibrary.smi')
    with open(inFile) as fIn:
      with open(oFile, 'w') as fO:
        for line in fIn:
          smi, smiName = line.split()[0].strip(), line.split()[1].strip()
          if smiName in molNames:
            if seqName:
              score = intDic[seqName][smiName]
            fO.write(f'{smi}\t{smiName}\t{score}\n')

    intMols.setFileName(oFile)
    intMols.calculateLength()
    self._defineOutputs(outputLibrary=intMols)

  def createOutputStep(self):
    inSeqs = self.inputSequences.get()

    filtSeqNames, filtMolNames = self.chooseSeq.get().strip().split(','), self.chooseMol.get().strip().split(',')
    filtScoreType = self.chooseScore.get()

    molNames = getFilteredOutput(inSeqs, filtSeqNames, filtMolNames, filtScoreType, self.scThres.get())[2]

    seqName = None
    if len(filtSeqNames) == 1 and filtSeqNames[0].strip() != 'All':
      seqName = filtSeqNames[0].strip()

    intMols = self.getInputMols()
    intDic = inSeqs.getInteractScoresDic()

    if isinstance(intMols, SetOfSmallMolecules):
      self.defineMolsOutput(intDic, molNames, seqName)

    elif isinstance(intMols, SmallMoleculesLibrary):
      self.defineLibraryOutput(intMols, intDic, molNames, seqName)


  ############## UTILS ########################
  def getInputMols(self):
    return self.inputSequences.get().getInteractMols()

  def getScoreOptions(self):
      return self.inputSequences.get().getScoreTypes()

