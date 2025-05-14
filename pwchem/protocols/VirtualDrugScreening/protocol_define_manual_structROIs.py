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
import os, json, math, sys
from scipy.spatial import distance
from scipy.cluster.hierarchy import linkage, fcluster
from Bio.PDB.ResidueDepth import get_surface
from Bio.PDB.PDBParser import PDBParser

from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol
from pwem.convert import cifToPdb
from pwem.objects import Pointer

from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import *
from pwchem import Plugin
from pwchem.constants import MGL_DIC

COORDS, RESIDUES, LIGANDS, PPI, NRES = 0, 1, 2, 3, 4

class ProtDefineStructROIs(EMProtocol):
    """
    Defines a set of structural ROIs from a set of coordinates / residues / predocked ligands
    """
    _label = 'Define structural ROIs'
    typeChoices = ['Coordinates', 'Residues', 'Ligand', 'Protein-Protein Interface', 'Near Residues']

    # -------------------------- DEFINE param functions ----------------------
    def _defineClusterParams(self, group, condition='True'):
      group.addParam('maxIntraDistance', params.FloatParam, default='2.0', condition=condition,
                     label='Maximum distance between pocket points (A): ',
                     help='Maximum distance between two pocket atoms to considered them same pocket')

      group.addParam('surfaceCoords', params.BooleanParam, default=True,
                     label='Map coordinates to surface? ',
                     help='Whether to map the input coordinates (from the residues, coordinates, or ligand) to the '
                          'closest surface coordinates or use them directly.')
      group.addParam('maxDepth', params.FloatParam, default='3.0',
                     label='Maximum atom depth (A): ', condition='surfaceCoords',
                     help='Maximum atom distance to the surface to be considered and mapped')
      return group

    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                      allowsNull=True, label="Input AtomStruct: ",
                      help='Select the AtomStruct object where the structural ROIs will be defined.'
                           'It only can be empty if the input is a SetOfSmallMolecules docked (with AtomStruct)')

        group = form.addGroup('Origin')
        group.addParam('origin', params.EnumParam, default=RESIDUES,
                      label='Extract ROIs from: ', choices=self.typeChoices,
                      help='The ROIs will be defined from a set of elements of this type')

        line = group.addLine('Coordinates: ', condition='origin=={}'.format(COORDS),
                             help='Coordinates of the defiend ROI')
        line.addParam('coordX', params.IntParam, default=10, label='X: ')
        line.addParam('coordY', params.IntParam, default=10, label='Y: ')
        line.addParam('coordZ', params.IntParam, default=10, label='Z: ')

        group.addParam('chain_name', params.StringParam,
                       allowsNull=False, label='Chain of interest:', condition='origin in [{}, {}]'.
                       format(RESIDUES, PPI),
                       help='Specify the chain of the residue of interest')
        group.addParam('resPosition', params.StringParam, condition='origin=={}'.format(RESIDUES),
                       allowsNull=False, label='Residues of interest',
                       help='Specify the residue to define a region of interest.\n'
                            'You can either select a single residue or a range '
                            '(it will take into account the first and last residues selected)')

        group.addParam('extLig', params.BooleanParam, default=True, condition='origin=={}'.format(LIGANDS),
                       label='Input is a external molecule? ',
                       help='Whether the ligand is docked in a external molecule or in the AtomStruct itself')

        group.addParam('inSmallMols', params.PointerParam, pointerClass='SetOfSmallMolecules',
                      condition='origin=={} and extLig'.format(LIGANDS), allowsNull=True,
                      label='Input molecules: ', help='Predocked molecules which will define the protein pockets')

        group.addParam('ligName', params.StringParam,
                       condition='origin=={} and extLig'.format(LIGANDS),
                       label='Ligand name: ', help='Specific ligand of the set')

        group.addParam('molName', params.StringParam,
                       condition='origin=={} and not extLig'.format(LIGANDS),
                       label='Molecule name: ', help='Name of the HETATM molecule in the AtomStruct')
        group.addParam('remMol', params.BooleanParam, default=False, label='Remove molecule from structure: ',
                       condition='origin=={} and not extLig'.format(LIGANDS),
                        help='Whether to remove the molecule from the output structure in which the ROIs are defined')

        group.addParam('chain_name2', params.StringParam,
                       allowsNull=False, label='Chain of interest 2: ', condition='origin=={}'.format(PPI),
                       help='Select the second protein chain to determine its interface with the first chain')
        group.addParam('ppiDistance', params.FloatParam, default='3.0',
                       label='Interface distance (A): ', condition='origin=={}'.format(PPI),
                       help='Maximum distance between two chain atoms to considered them part of the interface')

        group.addParam('resNRes', params.StringParam, default='', label='Residue pattern: ',
                       condition='origin=={}'.format(NRES),
                       help='Define a ROI by specifying the number of residues of a specific type that must be near '
                            'each other. Use the 3 letters code, comma separated. '
                            '\ne.g: 3 cysteines: "CYS, CYS, CYS"')
        group.addParam('resDistance', params.FloatParam, default='5.0', condition='origin=={}'.format(NRES),
                       label='Maximum distance between the residues (A): ',
                       help='Maximum distance between the center of mass of two residues to consider them part '
                            'of the ROI.')
        group.addParam('linkNRes', params.EnumParam, default=0, label='Clustering linkage type for near residues: ',
                       choices=['Single', 'Complete', 'Average', 'Centroid', 'Median', 'Ward'],
                       condition='origin=={}'.format(NRES),
                       help='Define the type of linkage oof the clustering that will define the Near Residues ROI.\n'
                            'Single linkage will group residues where each of them is closer than the threshold to at '
                            'least one other residue.\nComplete linkage will require all residue pairs to be closer '
                            'than the threshold, which will lead to more "globular" regions')

        group.addParam('addROI', params.LabelParam, label='Add defined ROI: ',
                       help='Use this wizard to store the defined ROI in the list')


        group.addParam('inROIs', params.TextParam, width=70, default='',
                       label='List of input ROIs: ',
                       help='Inputs to define the structural ROI.\n'
                            'The coordinates of the interface atoms will be mapped to surface points closer than '
                            'maxDepth and points closer than maxIntraDistance will be considered the same pocket')

        group = form.addGroup('Pocket definition')
        group = self._defineClusterParams(group)

        form.addSection(label='Input Pointers')
        form.addParam('inputPointerLabels', params.LabelParam, important=True,
                      label='Records of inputs. Do not modify manually',
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')
        form.addParam('inputPointers', params.MultiPointerParam, pointerClass="SetOfSmallMolecules",
                      label='Input Pointers: ', allowsNull=True,
                      help='This is a list of the input pointer to keep track of the inputs received.\n'
                           'It is automatically updated with the first section wizards.\n'
                           'Manual modification (adding inputs from the lens) will have no actual impact on the '
                           'protocol performance')

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
            pdbFile = cleanPDB(self.getProteinFileName(), os.path.abspath(self._getExtraPath('cleanPDB.pdb')),
                                waters=True, hetatm=True)

        structure = parser.get_structure(self.getProteinName(), pdbFile)
        self.structModel = structure[0]

        self.structSurface = get_surface(self.structModel,
                                         MSMS=Plugin.getProgramHome(MGL_DIC, 'MGLToolsPckgs/binaries/msms'))

    def definePocketsStep(self):
        self.coordsClusters = []
        for roiStr in self.inROIs.get().split('\n'):
            if roiStr.strip():
                pocketCoords = self.getOriginCoords(roiStr)
                if self.surfaceCoords:
                    pocketCoords = self.mapSurfaceCoords(pocketCoords)

                self.coordsClusters += clusterSurfaceCoords(pocketCoords, self.maxIntraDistance.get())

    def defineOutputStep(self):
        inpStruct = self.inputAtomStruct.get()
        pdbFile = self.getProteinFileName()
        newPDBFile = self._getPath(os.path.basename(pdbFile))
        newPDBFile = self.checkLigands(pdbFile, newPDBFile)

        outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for i, clust in enumerate(self.coordsClusters):
            pocketFile = createPocketFile(clust, i, outDir=self._getExtraPath())
            pocket = StructROI(pocketFile, newPDBFile)
            if str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
                pocket._maeFile = String(os.path.relpath(inpStruct.getFileName()))
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
        if not self.getProteinFileName(inProtocol=False):
            errors.append('You must specify an input structure')

        for roiStr in self.inROIs.get().split('\n'):
          if roiStr.strip():
            jSonStr = ':'.join(roiStr.split(':')[1:]).strip()
            jDic = json.loads(jSonStr)
            roiKey = roiStr.split()[1].strip()
            if roiKey == 'Near_Residues:':
              residueTypes = jDic['residues'].split(',')
              residueTypes = [res.strip().upper() for res in residueTypes]
              for res in residueTypes:
                  if not res in RESIDUES3TO1:
                      errors.append('Cannot recognize residue {}'.format(res))

        return errors

    # --------------------------- UTILS functions -----------------------------------
    def checkLigands(self, pdbFile, newPDBFile):
      ligToRem = []
      for roiStr in self.inROIs.get().split('\n'):
        if roiStr.strip():
          jSonStr = ':'.join(roiStr.split(':')[1:]).strip()
          jDic = json.loads(jSonStr)
          roiKey = roiStr.split()[1].strip()
          if roiKey == 'Ligand:' and jDic['remove'] == 'True':
            ligToRem.append(jDic['molName'])

      if ligToRem:
        newPDBFile = cleanPDB(pdbFile, newPDBFile, het2rem=ligToRem)
      else:
        shutil.copy(pdbFile, newPDBFile)
      return newPDBFile

    def getProteinFileName(self, inProtocol=True):
        inpStruct = self.inputAtomStruct.get()
        if inpStruct:
            inpFile = inpStruct.getFileName()
        else:
            inpFile = self.inSmallMols.get().getProteinFile()

        if not inpFile:
            return None

        if inpFile.endswith('.cif'):
          if inProtocol:
              inpPDBFile = self._getExtraPath(os.path.basename(inpFile).replace('.cif', '.pdb'))
          else:
              inpPDBFile = self.getProject().getTmpPath(os.path.basename(inpFile).replace('.cif', '.pdb'))
          cifToPdb(inpFile, inpPDBFile)

        elif str(type(inpStruct).__name__) == 'SchrodingerAtomStruct':
            inpPDBFile = self.getProject().getTmpPath(os.path.basename(inpFile).
                                                      replace(inpStruct.getExtension(), '.pdb'))
            inpStruct.convert2PDB(outPDB=inpPDBFile)

        else:
          inpPDBFile = inpFile

        return inpPDBFile

    def getProteinName(self):
      return os.path.splitext(os.path.basename(self.getProteinFileName()))[0]

    def getPrevPointersIds(self, prevPointers):
      ids = []
      for p in prevPointers:
        ids.append(p.get().getObjId())
      return ids

    def getNewPointers(self):
      prevPointers = self.inputPointers
      prevIds = self.getPrevPointersIds(prevPointers)
      newObj = self.inSmallMols.get()
      newId = newObj.getObjId()
      if newId not in prevIds:
        newIndex = len(prevPointers)
        prevPointers.append(Pointer(newObj))
      else:
        newIndex = prevIds.index(newId)
      return newIndex, prevPointers

    def getDefinedROILine(self):
      if self.origin.get() == COORDS:
        roiDef = f'Coordinate: {{"X": {self.coordX.get()}, "Y": {self.coordY.get()}, "Z": {self.coordZ.get()}}}\n'

      elif self.origin.get() == RESIDUES:
        chainDic, resDic = json.loads(self.chain_name.get()), json.loads(self.resPosition.get())
        roiDef = f'Residues: {{"model": {chainDic["model"]}, "chain": "{chainDic["chain"]}", ' \
                 f'"index": "{resDic["index"]}", "residues": "{resDic["residues"]}"}}\n'

      elif self.origin.get() == LIGANDS:
        if not self.extLig.get():
          roiDef = f'Ligand: {{"molName": "{self.molName.get()}", "remove": "{self.remMol.get()}"}}\n'
        else:
          newIndex, _ = self.getNewPointers()
          roiDef = f'Ext-Ligand: {{"pointerIdx": "{newIndex}", "ligName": "{self.ligName.get()}"}}\n'

      elif self.origin.get() == PPI:
        chain1Dic, chain2Dic = json.loads(self.chain_name.get()), json.loads(self.chain_name2.get())
        iDist = self.ppiDistance.get()

        c1, c2 = '{}-{}'.format(chain1Dic['model'], chain1Dic['chain']), \
                 '{}-{}'.format(chain2Dic['model'], chain2Dic['chain'])

        roiDef = f'PPI: {{"chain1": "{c1}", "chain2": "{c2}", "interDist": "{iDist}"}}\n'

      elif self.origin.get() == NRES:
        resTypes, resDist = self.resNRes.get(), self.resDistance.get()
        resLink = self.getEnumText("linkNRes")
        roiDef = f'Near_Residues: {{"residues": "{resTypes}", "distance": "{resDist}", "linkage": "{resLink}"}}\n'
      return roiDef

    def getOriginCoords(self, roiStr):
        oCoords = []
        jSonStr = ':'.join(roiStr.split(':')[1:]).strip()
        jDic = json.loads(jSonStr)
        roiKey = roiStr.split()[1].strip()
        if roiKey == 'Coordinate:':
            oCoords = [[jDic[c] for c in jDic]]

        elif roiKey == 'Residues:':
            chainId, resIdxs = jDic['chain'], jDic['index']
            idxs = [int(resIdxs.split('-')[0]), int(resIdxs.split('-')[1])]
            for resId in range(idxs[0], idxs[1] + 1):
                residue = self.structModel[chainId][resId]
                atoms = residue.get_atoms()
                for a in atoms:
                  oCoords.append(list(a.get_coord()))

        elif roiKey == 'Ligand:':
            oCoords = getLigCoords(self.getProteinFileName(), jDic['molName'])

        elif roiKey == 'Ext-Ligand:':
            mol = None
            molSet = self.inputPointers[int(jDic['pointerIdx'])].get()
            for mol in molSet:
              if jDic['ligName'] == mol.__str__():
                break
            if mol:
              curPosDic = mol.getAtomsPosDic(onlyHeavy=False)
              oCoords += list(curPosDic.values())

        elif roiKey == 'PPI:':
            chain1Id, chain2Id = jDic['chain1'].split('-')[1], jDic['chain2'].split('-')[1]
            ppiDist = float(jDic['interDist'])
            chain1, chain2 = self.structModel[chain1Id], self.structModel[chain2Id]

            print('Checking interface between chains "{}" and "{}"'.format(chain1Id, chain2Id))
            sys.stdout.flush()
            for atom1 in chain1.get_atoms():
                for atom2 in chain2.get_atoms():
                    coord1, coord2 = list(atom1.get_coord()), list(atom2.get_coord())
                    dist = math.dist(coord1, coord2)
                    if dist <= ppiDist:
                        oCoords.append(coord1)
                        oCoords.append(coord2)
                        break

        elif roiKey == 'Near_Residues:':
            residueTypes, resDist, resLink = jDic['residues'].split(','), jDic['distance'], jDic['linkage']
            residueTypes = [res.strip().upper() for res in residueTypes]
            cMassResidues = {}
            for res in self.structModel.get_residues():
                if res.get_resname() in residueTypes:
                    cMassResidues[res] = res.center_of_mass()

            if len(cMassResidues) > 0:
                linked = linkage(list(cMassResidues.values()), resLink.lower())
                clusterIdxs = fcluster(linked, resDist, 'distance')

                clusters = {}
                for i, res in enumerate(cMassResidues):
                    resIdx = clusterIdxs[i]
                    if resIdx in clusters:
                        clusters[resIdx] += [res]
                    else:
                        clusters[resIdx] = [res]

                for _, clust in clusters.items():
                    resNameList = [res.get_resname() for res in clust]
                    if self.isSublist(residueTypes, resNameList):
                        oCoords += list(a.get_coord() for res in clust for a in res.get_atoms())
                        self.addToPatternFile(clust)
                        for res in clust:
                            print(res, cMassResidues[res])

        return oCoords

    def getPatternFile(self):
        return os.path.abspath(self._getExtraPath('patternFile.txt'))
    def addToPatternFile(self, resList):
        pFile, openMode = self.getPatternFile(), 'a'
        if not os.path.exists(pFile):
            openMode = 'w'

        with open(pFile, openMode) as f:
            resIds = ['{}_{}'.format(res.get_resname(), res.get_id()[1]) for res in resList]
            f.write('{}\n'.format(','.join(resIds)))

    def mapSurfaceCoords(self, oCoords):
        sCoords = []
        for coord in oCoords:
            closerSCoords = self.closerSurfaceCoords(coord)
            for cCoord in closerSCoords:
                cCoord = list(cCoord)
                if not cCoord in sCoords:
                  sCoords.append(cCoord)

        return sCoords

    def isSublist(self, subList, lst):
      isSub = True
      countSub = {a: subList.count(a) for a in subList}
      for item, count in countSub.items():
          if not item in lst or lst.count(item) < countSub[item]:
              isSub = False
              break
      return isSub

    def closerSurfaceCoords(self, coord):
        distances = distance.cdist([coord], self.structSurface)
        closestIndexes = distances < self.maxDepth.get()
        return list(self.structSurface[closestIndexes[0]])

    def createPocketFile(self, clust, i):
        outFile = self._getExtraPath('pocketFile_{}.pdb'.format(i))
        with open(outFile, 'w') as f:
            for j, coord in enumerate(clust):
                f.write(writePDBLine(['HETATM', str(j), 'APOL', 'STP', 'C', '1', *coord, 1.0, 0.0, '', 'Ve']))
        return outFile




