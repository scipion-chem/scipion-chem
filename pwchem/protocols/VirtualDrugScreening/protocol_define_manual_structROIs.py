# -*- coding: utf-8 -*-
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
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************


"""
This protocol is used to manually define structural regions from coordinates, residues or docked small molecules

"""
import os, json
from scipy.spatial import distance
from Bio.PDB.ResidueDepth import ResidueDepth, get_surface, min_dist, residue_depth
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.object import String
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from pwchem import Plugin
from pwchem.constants import MGL_DIC

COORDS, RESIDUES, LIGANDS = 0, 1, 2

class ProtDefineStructROIs(EMProtocol):
    """
    Defines a set of structural ROIs from a set of coordinates / residues / predocked ligands
    """
    _label = 'Define structural ROIs'
    typeChoices = ['Coordinates', 'Residues', 'Ligand']

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      allowsNull=False, label="Input AtomStruct: ",
                      help='Select the AtomStruct object where the structural ROIs will be defined')
        group.addParam('origin', params.EnumParam, default=RESIDUES,
                      label='Extract ROIs from: ', choices=self.typeChoices,
                      help='The ROIs will be defined from a set of elements of this type')

        group.addParam('inCoords', params.TextParam, default='', width=70,
                      condition='origin=={}'.format(COORDS), label='Input coordinates: ',
                      help='Input coordinates to define the structural ROI. '
                           'The coordinates will be mapped to surface points closer than maxDepth '
                           'and points closer than maxIntraDistance will be considered the same ROI.'
                           '\nCoordinates in format: (1,2,3);(4,5,6);...')

        group.addParam('chain_name', params.StringParam,
                      allowsNull=False, label='Chain of interest', condition='origin=={}'.format(RESIDUES),
                      help='Specify the chain of the residue of interest')
        group.addParam('resPosition', params.StringParam, condition='origin=={}'.format(RESIDUES),
                      allowsNull=False, label='Residues of interest',
                      help='Specify the residue to define a region of interest.\n'
                           'You can either select a single residue or a range '
                           '(it will take into account the first and last residues selected)')
        group.addParam('addResidue', params.LabelParam,
                      label='Add defined residue', condition='origin=={}'.format(RESIDUES),
                      help='Here you can define a residue which will be added to the list of residues below.')
        group.addParam('inResidues', params.TextParam, width=70, default='',
                      condition='origin=={}'.format(RESIDUES), label='Input residues: ',
                      help='Input residues to define the structural ROI. '
                           'The coordinates of the residue atoms will be mapped to surface points closer than maxDepth '
                           'and points closer than maxIntraDistance will be considered the same pocket')

        group.addParam('extLig', params.BooleanParam, default=True, condition='origin=={}'.format(LIGANDS),
                       label='Input is a external molecule? ',
                       help='Whether the ligand is docked in a external molecule or in the AtomStruct itself')

        group.addParam('inSmallMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      condition='origin=={} and extLig'.format(LIGANDS),
                      label='Input molecules: ', help='Predocked molecules which will define the protein pockets')

        group.addParam('ligName', params.StringParam,
                       condition='origin=={} and extLig'.format(LIGANDS),
                       label='Ligand name: ', help='Specific ligand of the set')

        group.addParam('molName', params.StringParam,
                       condition='origin=={} and not extLig'.format(LIGANDS),
                       label='Molecule name: ', help='Name of the HETATM molecule in the AtomStruct')

        group = form.addGroup('Distances')
        group.addParam('maxDepth', params.FloatParam, default='3.0',
                      label='Maximum atom depth (A): ',
                      help='Maximum atom distance to the surface to be considered and mapped')
        group.addParam('maxIntraDistance', params.FloatParam, default='2.0',
                      label='Maximum distance between pocket points (A): ',
                      help='Maximum distance between two pocket atoms to considered them same pocket')



    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        # Insert processing steps
        self._insertFunctionStep('getSurfaceStep')
        self._insertFunctionStep('definePocketsStep')
        self._insertFunctionStep('defineOutputStep')

    def getSurfaceStep(self):
        parser = PDBParser()
        pdbFile = self.getProteinFileName()
        if self.origin.get() == LIGANDS and not self.extLig:
            pdbFile = clean_PDB(self.getProteinFileName(), os.path.abspath(self._getExtraPath('cleanPDB.pdb')),
                                waters=True, HETATM=True)

        structure = parser.get_structure(self.getProteinName(), pdbFile)
        self.structModel = structure[0]
        self.structSurface = get_surface(self.structModel,
                                         MSMS=Plugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))

    def definePocketsStep(self):
        originCoords = self.getOriginCoords()
        surfaceCoords = self.mapSurfaceCoords(originCoords)

        self.coordsClusters = self.clusterSurfaceCoords(surfaceCoords)

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for i, clust in enumerate(self.coordsClusters):
            pocketFile = self.createPocketFile(clust, i)
            pocket = StructROI(pocketFile, self.getProteinFileName())
            if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                pocket._maeFile = String(os.path.abspath(inpStruct.getFileName()))
            pocket.calculateContacts()
            outPockets.append(pocket)

        if len(outPockets) > 0:
            outPockets.buildPDBhetatmFile()
            self._defineOutputs(outputStructROIs=outPockets)



    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------

    def getProteinFileName(self):
        inpStruct = self.inputAtomStruct.get()
        inpFile = inpStruct.getFileName()
        if inpFile.endswith('.cif'):
          inpPDBFile = self._getExtraPath(os.path.basename(inpFile).replace('.cif', '.pdb'))
          cifToPdb(inpFile, inpPDBFile)

        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpPDBFile = self._getExtraPath(os.path.basename(inpFile).replace(inpStruct.getExtension(), '.pdb'))
            inpStruct.convert2PDB(outPDB=inpPDBFile)

        else:
          inpPDBFile = inpFile

        return inpPDBFile

    def getProteinName(self):
      return os.path.splitext(os.path.basename(self.getProteinFileName()))[0]

    def getOriginCoords(self):
        oCoords = []
        if self.origin.get() == COORDS:
            coordsStr = self.inCoords.get().split(';')
            oCoords = [eval(c) for c in coordsStr]

        elif self.origin.get() == RESIDUES:
            residuesStr = self.inResidues.get().strip().split('\n')
            for rStr in residuesStr:
                resDic = json.loads(rStr)
                chainId, resIdxs = resDic['chain'], resDic['index']
                idxs = [int(resIdxs.split('-')[0]), int(resIdxs.split('-')[1])]
                for resId in range(idxs[0], idxs[1] + 1):
                    residue = self.structModel[chainId][resId]
                    atoms = residue.get_atoms()
                    for a in atoms:
                      oCoords.append(list(a.get_coord()))

        elif self.origin.get() == LIGANDS:
            if self.extLig:
                if not self.ligName.get():
                    mols = self.inSmallMols.get()
                else:
                    for mol in self.inSmallMols.get():
                        if self.ligName.get() == mol.__str__():
                            mols = [mol]
                            break

                for mol in mols:
                    curPosDic = mol.getAtomsPosDic(onlyHeavy=False)
                    oCoords += list(curPosDic.values())

            else:
                oCoords = self.getLigCoords(self.getProteinFileName(), self.molName.get())

        return oCoords

    def mapSurfaceCoords(self, oCoords):
        sCoords = []
        for coord in oCoords:
            closerSCoords = self.closerSurfaceCoords(coord)
            for cCoord in closerSCoords:
                cCoord = list(cCoord)
                if not cCoord in sCoords:
                  sCoords.append(cCoord)

        return sCoords


    def closerSurfaceCoords(self, coord):
        distances = distance.cdist([coord], self.structSurface)
        closestIndexes = distances < self.maxDepth.get()
        return list(self.structSurface[closestIndexes[0]])

    def clusterSurfaceCoords(self, surfCoords):
        clusters = []
        for coord in surfCoords:
            newClusters = []
            newClust = [coord]
            for clust in clusters:
                merge = False
                for cCoord in clust:
                    dist = calculateDistance(coord, cCoord)
                    if dist < self.maxIntraDistance.get():
                        merge = True
                        break

                if merge:
                    newClust += clust
                else:
                    newClusters.append(clust)

            newClusters.append(newClust)
            clusters = newClusters.copy()
        return clusters

    def createPocketFile(self, clust, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(clust):
                f.write(writePDBLine(['HETATM', str(j), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
        return outFile

    def getLigCoords(self, ASFile, ligName):
        if ASFile.endswith('.pdb') or ASFile.endswith('.ent'):
          pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
          parser = PDBParser().get_structure(pdb_code, ASFile)
        elif ASFile.endswith('.cif'):
          pdb_code = os.path.basename(os.path.splitext(ASFile)[0])
          parser = MMCIFParser().get_structure(pdb_code, ASFile)
        else:
          print('Unknown AtomStruct file format')
          return

        coords = []
        for model in parser:
            for chain in model:
                for residue in chain:
                    if residue.resname == ligName:
                        for atom in residue:
                            coords.append(list(atom.get_coord()))
        return coords




