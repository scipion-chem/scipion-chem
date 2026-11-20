# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:   Blanca Pueche (blanca.pueche@cnb.csic.es)
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

import os
import json
import csv
import numpy as np

from Bio.PDB import PDBParser, MMCIFParser

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SmallMolecule, SetOfSmallMolecules
from pwchem.utils import pdbqt2other
from pwem.convert import cifToPdb
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC, OPENBABEL_DIC


class ProtContactLigands(EMProtocol):
    """

    """
    _label = 'Filter contact ligands'

    def _defineParams(self, form):
        form.addSection(label='Input')

        form.addParam('inputLigandSet', params.PointerParam,
                      label="Input docked ligands: ", allowsNull=False,
                      important=True, pointerClass='SetOfSmallMolecules',
                      help='Set of docked ligands to filter.')

        form.addParam('recChain', params.StringParam,
                      label="Receptor chain: ", default='',
                      help='Chain of the receptor that we want to evaluate.')
        form.addParam('recResidue', params.StringParam,
                      label="Receptor residue:", default='',
                      help="Select the receptor residue to evaluate." )
        form.addParam('recAtom', params.StringParam,
                      label="Receptor atom:", default='',
                      help="Select the receptor atom to evaluate. Enter the atom name displayed in PyMOL, such as OG, CA, or OD1, not the sequential label number.")
        form.addParam('ligAtom', params.IntParam,
                      label="Ligand atom: ", default=0,
                      help='Atom of the ligand that we want to evaluate.')

        form.addParam('cutoff', params.FloatParam,
                      label="Distance cutoff (Å): ", default=4,
                      help='Distance cutoff for evaluation.')

    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('filterLigandsStep')
        self._insertFunctionStep('createOutputStep')

    def filterLigandsStep(self):
        """Calculate atom-to-atom distances and retain passing ligands."""

        inputSet = self.inputLigandSet.get()
        receptorFile = inputSet.getProteinFile()

        if not receptorFile or not os.path.isfile(receptorFile):
            raise ValueError(
                'The input ligand set has no accessible receptor file.'
            )

        cutoff = float(self.cutoff.get())
        ligandAtomIndex = int(self.ligAtom.get())

        if ligandAtomIndex < 1:
            raise ValueError('Ligand atom index must be >= 1.')

        if cutoff <= 0:
            raise ValueError('Distance cutoff must be positive.')

        referenceFile = self._getReferenceLigandFile(inputSet)

        if not os.path.isfile(referenceFile):
            raise ValueError(
                'Reference ligand file does not exist: {}'.format(
                    referenceFile
                )
            )

        self.info(
            'Reference ligand for atom mapping: {}'.format(referenceFile)
        )

        receptorFile = self.convertFormat(receptorFile)
        receptor = self._loadStructure(receptorFile)
        targetResidue = self._getReceptorResidue(receptor)

        receptorAtomCoordinates = (
            self._getReceptorAtomCoordinates(targetResidue)
        )

        self._passingMolecules = []
        reportFile = self._getPath('contact_distances.csv')

        with open(reportFile, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)

            writer.writerow([
                'Molecule',
                'Pose file',
                'Reference ligand file',
                'Receptor chain',
                'Receptor residue',
                'Receptor atom',
                'Reference heavy-atom index',
                'Mapped pose atom index',
                'Mapped ligand element',
                'Ligand atom coordinates (A)',
                'Receptor atom coordinates (A)',
                'Distance (A)',
                'Cutoff (A)',
                'Retained',
                'Error'
            ])

            for molecule in inputSet:
                molName = molecule.getMolName()
                poseFile = ''
                ligandCoordinates = None
                ligandElement = ''
                poseAtomIndex = ''
                distance = None
                errorMessage = ''
                keep = False

                try:
                    poseFile = self._getPoseFile(molecule)

                    (
                        ligandCoordinates,
                        ligandElement,
                        poseAtomIndex
                    ) = self._getMappedLigandAtomCoordinates(
                        referenceFile,
                        poseFile
                    )

                    distance = self._calculateDistance(
                        ligandCoordinates,
                        receptorAtomCoordinates
                    )

                    keep = distance <= cutoff

                    if keep:
                        self._passingMolecules.append(molecule.clone())

                    print(
                        '{} | reference heavy atom {} | '
                        'pose atom {} ({}) at {} | '
                        'receptor atom {} at {} | '
                        'distance={:.3f} A | retained={}'.format(
                            os.path.basename(poseFile),
                            ligandAtomIndex,
                            poseAtomIndex,
                            ligandElement,
                            np.round(ligandCoordinates, 3).tolist(),
                            self.recAtom.get(),
                            np.round(
                                receptorAtomCoordinates, 3
                            ).tolist(),
                            distance,
                            keep
                        )
                    )

                except Exception as exc:
                    errorMessage = str(exc)
                    keep = False

                    self.warning(
                        'Could not evaluate {}: {}'.format(
                            molName, errorMessage
                        )
                    )

                writer.writerow([
                    molName,
                    poseFile,
                    referenceFile,
                    self.recChain.get(),
                    self.recResidue.get(),
                    self.recAtom.get(),
                    ligandAtomIndex,
                    poseAtomIndex,
                    ligandElement,
                    (
                        np.round(ligandCoordinates, 3).tolist()
                        if ligandCoordinates is not None else ''
                    ),
                    np.round(
                        receptorAtomCoordinates, 3
                    ).tolist(),
                    (
                        '{:.3f}'.format(distance)
                        if distance is not None else ''
                    ),
                    cutoff,
                    'Yes' if keep else 'No',
                    errorMessage
                ])

        self.info(
            'Contact filtering completed. Results saved to {}.'.format(
                reportFile
            )
        )

    def createOutputStep(self):
        """Create the output set containing the retained molecules."""

        inputSet = self.inputLigandSet.get()

        outputSet = SetOfSmallMolecules().create(
            outputPath=self._getPath()
        )

        outputSet.copyInfo(inputSet)
        outputSet.setProteinFile(inputSet.getProteinFile())
        outputSet.setDocked(True)

        for molecule in self._passingMolecules:
            outputSet.append(molecule)

        self._defineOutputs(outputSmallMolecules=outputSet)

        self._defineSourceRelation(self.inputLigandSet, outputSet)

        self.info(
            'Contact filtering completed: {} molecules retained out of {}.'
            .format(
                len(self._passingMolecules),
                inputSet.getSize()
            )
        )

        # --------------------------- INFO functions -----------------------

    def _summary(self):
        return []

    def _methods(self):
        return []

    def _validate(self):
        validations = []

        molSet = self.inputLigandSet.get()

        if not molSet.isDocked():
            validations.append(
                '{} is not docked yet.'.format(molSet)
            )

        if self.cutoff.get() <= 0:
            validations.append('Cutoff must be positive.')

        if self.ligAtom.get() < 1:
            validations.append(
                'Ligand atom index must be greater than or equal to 1.'
            )

        if not self.recAtom.get().strip():
            validations.append(
                'Please specify a receptor atom name.'
            )

        return validations

    def _warnings(self):
        return []

        # --------------------------- Utility functions --------------------

    @staticmethod
    def _loadStructure(filename):
        """Read a PDB or mmCIF structure."""

        ext = os.path.splitext(filename)[1].lower()

        if ext in ('.pdb', '.ent'):
            parser = PDBParser(QUIET=True)

        elif ext in ('.cif', '.mmcif'):
            parser = MMCIFParser(QUIET=True)

        else:
            raise ValueError(
                'Unsupported structure format: {}'.format(filename)
            )

        return parser.get_structure(
            os.path.basename(filename),
            filename
        )

    def _getReceptorResidue(self, structure):
        """Find the selected receptor residue using the wizard outputs."""

        try:
            chainData = json.loads(self.recChain.get())
            residueData = json.loads(self.recResidue.get())

            chainId = str(chainData['chain'])
            residueIndex = str(
                residueData['index']
            ).split('-')[0]

            residueName = residueData.get('residue', '')

        except (TypeError, ValueError, KeyError) as exc:
            raise ValueError(
                'Could not interpret the selected receptor chain or '
                'residue. Check the selections in the wizard. {}'.format(
                    exc
                )
            )

        matches = []

        for model in structure:
            if chainId not in model:
                continue

            chain = model[chainId]

            for residue in chain:
                heteroFlag, resseq, icode = residue.id

                if (
                        str(resseq) == residueIndex
                        and heteroFlag == ' '
                        and (
                        not residueName
                        or residue.resname == residueName
                )
                ):
                    matches.append(residue)

        if not matches:
            raise ValueError(
                'Residue {} in chain {} was not found in the receptor.'
                .format(residueIndex, chainId)
            )

        return matches[0]

    def _getReceptorAtomCoordinates(self, residue):
        """Return coordinates of the selected heavy atom in the residue."""

        atomName = self.recAtom.get().strip()

        if not atomName:
            raise ValueError(
                'Please specify a receptor atom name.'
            )

        matches = [
            atom for atom in residue.get_atoms()
            if atom.get_name().strip() == atomName
        ]

        if len(matches) != 1:
            availableAtoms = [
                atom.get_name().strip()
                for atom in residue.get_atoms()
            ]

            raise ValueError(
                'Expected exactly one atom named "{}" in residue {}, '
                'but found {}. Available atoms: {}'.format(
                    atomName,
                    residue.resname,
                    len(matches),
                    ', '.join(availableAtoms)
                )
            )

        atom = matches[0]
        element = (atom.element or '').strip().upper()

        if element in ('H', 'D'):
            raise ValueError(
                'Selected receptor atom "{}" is hydrogen or deuterium.'
                .format(atomName)
            )

        return np.asarray(atom.coord, dtype=float)

    @staticmethod
    def _getPoseFile(molecule):
        """Obtain the pose file associated with a SmallMolecule."""

        poseFile = molecule.getPoseFile()

        if poseFile and os.path.isfile(poseFile):
            return poseFile

        # Fall back to the molecule's main file.
        molFile = molecule.getFileName()

        if molFile and os.path.isfile(molFile):
            return molFile

        raise ValueError(
            'Could not locate a valid pose file for molecule {}.'
            .format(molecule.getMolName())
        )

    def _getLigandAtomCoordinates(self, poseFile):
        """
        Obtain the coordinates of the selected ligand heavy atom.

        The ligand atom is identified by its 1-based position among
        heavy atoms, excluding hydrogen and deuterium.

        Supported formats:
            PDB, ENT, CIF, mmCIF, PDBQT, SDF and MOL (V2000).

        IMPORTANT:
            The pose file must contain the ligand alone, or at least
            list the ligand atoms separately in a way consistent with
            the atom order used by the PyMOL wizard.
        """

        ext = os.path.splitext(poseFile)[1].lower()
        atomIndex = int(self.ligAtom.get())

        if atomIndex < 1:
            raise ValueError(
                'Ligand atom index must be >= 1.'
            )

        # PDB and mmCIF: use the heavy-atom order in the structure file.
        if ext in ('.pdb', '.ent', '.cif', '.mmcif'):
            structure = self._loadStructure(poseFile)

            heavyAtoms = [
                atom for atom in structure.get_atoms()
                if (atom.element or '').strip().upper() not in ('H', 'D')
            ]

            if atomIndex > len(heavyAtoms):
                raise ValueError(
                    'Heavy-atom index {} exceeds the number of heavy '
                    'atoms in {} ({} heavy atoms found).'.format(
                        atomIndex,
                        poseFile,
                        len(heavyAtoms)
                    )
                )

            atom = heavyAtoms[atomIndex - 1]

            return np.asarray(atom.coord, dtype=float)

        # PDBQT: parse atom coordinates and atom types directly.
        elif ext == '.pdbqt':
            heavyAtomCoordinates = []

            with open(poseFile, 'r') as handle:
                for line in handle:
                    if not line.startswith(('ATOM  ', 'HETATM')):
                        continue

                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except (ValueError, IndexError):
                        raise ValueError(
                            'Could not parse coordinates in {}: {}'.format(
                                poseFile, line.rstrip()
                            )
                        )

                    fields = line.split()

                    if not fields:
                        continue

                    # The final PDBQT field is the AutoDock atom type.
                    atom_type = fields[-1].strip().upper()

                    # Examples of hydrogen types include H and HD.
                    if atom_type in ('H', 'HD', 'D'):
                        continue

                    heavyAtomCoordinates.append([x, y, z])

            if atomIndex > len(heavyAtomCoordinates):
                raise ValueError(
                    'Heavy-atom index {} exceeds the number of heavy '
                    'atoms in {} ({} heavy atoms found).'.format(
                        atomIndex,
                        poseFile,
                        len(heavyAtomCoordinates)
                    )
                )

            return np.asarray(
                heavyAtomCoordinates[atomIndex - 1],
                dtype=float
            )

        # SDF/MOL: parse the V2000 atom block.
        elif ext in ('.sdf', '.mol'):
            with open(poseFile, 'r') as handle:
                lines = handle.readlines()

            if len(lines) < 4:
                raise ValueError(
                    'Invalid or incomplete MOL/SDF file: {}'.format(
                        poseFile
                    )
                )

            counts_line = lines[3]

            if 'V3000' in counts_line:
                raise ValueError(
                    'V3000 MOL/SDF files are not supported by this '
                    'parser: {}. Convert the file to V2000 first.'
                    .format(poseFile)
                )

            try:
                n_atoms = int(counts_line[0:3])
            except ValueError:
                raise ValueError(
                    'Could not read the V2000 atom count from {}.'
                    .format(poseFile)
                )

            atom_lines = lines[4:4 + n_atoms]

            if len(atom_lines) != n_atoms:
                raise ValueError(
                    'Atom block is incomplete in {}.'.format(poseFile)
                )

            heavyAtomCoordinates = []

            for i, line in enumerate(atom_lines, start=1):
                try:
                    x = float(line[0:10])
                    y = float(line[10:20])
                    z = float(line[20:30])
                    element = line[31:34].strip().upper()

                except (ValueError, IndexError):
                    raise ValueError(
                        'Could not parse an atom record in {}: {}'.format(
                            poseFile, line.rstrip()
                        )
                    )

                # Print the original SDF atom index, including hydrogens.
                print(
                    'SDF INDEX={:<3} ELEMENT={:<2} COORDS={}'.format(
                        i,
                        element,
                        np.round([x, y, z], 3).tolist()
                    )
                )

                if element in ('H', 'D'):
                    continue

                heavyAtomCoordinates.append([x, y, z])

            if atomIndex > len(heavyAtomCoordinates):
                raise ValueError(
                    'Heavy-atom index {} exceeds the number of heavy '
                    'atoms in {} ({} heavy atoms found).'.format(
                        atomIndex,
                        poseFile,
                        len(heavyAtomCoordinates)
                    )
                )

            return np.asarray(
                heavyAtomCoordinates[atomIndex - 1],
                dtype=float
            )

        else:
            raise ValueError(
                'Unsupported ligand pose format: {}'.format(poseFile)
            )

    @staticmethod
    def _calculateDistance(ligandCoordinates, receptorCoordinates):
        """Calculate the Euclidean distance between two atoms in Angstroms."""

        differences = receptorCoordinates - ligandCoordinates

        return float(np.linalg.norm(differences))

    def convertFormat(self, filename):
        """Convert PDBQT or CIF receptor files to PDB when necessary."""

        ext = os.path.splitext(filename)[1].lower()

        if ext == '.pdbqt':
            outputFile = os.path.abspath(
                self._getExtraPath(
                    os.path.splitext(
                        os.path.basename(filename)
                    )[0] + '.pdb'
                )
            )

            pdbqt2other(self, filename, outputFile)

            return outputFile

        elif ext in ('.cif', '.mmcif'):
            outputFile = os.path.abspath(
                self._getExtraPath(
                    os.path.splitext(
                        os.path.basename(filename)
                    )[0] + '.pdb'
                )
            )

            cifToPdb(filename, outputFile)

            return outputFile

        return filename

    def _getMappedLigandAtomCoordinates(self, referenceFile, poseFile):
        referenceFile = self._getRDKitCompatibleLigandFile(referenceFile)
        poseFile = self._getRDKitCompatibleLigandFile(poseFile)

        outputFile = self._getTmpPath('mapped_ligand_atom.json')

        args = (
            '--reference "{}" '
            '--pose "{}" '
            '--atom-index {} '
            '--output "{}"'
        ).format(
            os.path.abspath(referenceFile),
            os.path.abspath(poseFile),
            int(self.ligAtom.get()),
            os.path.abspath(outputFile)
        )

        pwchemPlugin.runScript(
            self,
            'rdkit_map_ligand_atom.py',
            args,
            env=RDKIT_DIC,
            popen=True,
            cwd=os.path.dirname(os.path.abspath(poseFile))
        )

        # Wait for the external script to finish before reading its output.
        import time

        timeout = 60
        elapsed = 0

        while not os.path.isfile(outputFile) and elapsed < timeout:
            time.sleep(0.2)
            elapsed += 0.2

        if not os.path.isfile(outputFile):
            raise RuntimeError(
                'RDKit atom-mapping script did not produce its output file: '
                '{}'.format(outputFile)
            )

        with open(outputFile, 'r') as handle:
            result = json.load(handle)

        return (
            np.asarray(result['coordinates'], dtype=float),
            result['element'],
            result['pose_atom_index']
        )

    @staticmethod
    def _getReferenceLigandFile(inputSet):
        """Get the ligand file used as reference for atom numbering."""

        for molecule in inputSet:
            poseFile = molecule.getPoseFile()

            if poseFile and os.path.isfile(poseFile):
                return poseFile

            molFile = molecule.getFileName()

            if molFile and os.path.isfile(molFile):
                return molFile

        raise ValueError('Could not find a reference ligand file.')

    def _getRDKitCompatibleLigandFile(self, ligandFile):
        """Return an SDF file for CIF inputs, otherwise the original file."""

        ext = os.path.splitext(ligandFile)[1].lower()

        if ext not in ('.cif', '.mmcif'):
            return ligandFile

        outputFile = os.path.abspath(
            self._getTmpPath(
                os.path.splitext(os.path.basename(ligandFile))[0] + '.sdf'
            )
        )

        if os.path.isfile(outputFile) and os.path.getsize(outputFile) > 0:
            return outputFile

        args = '-i "{}" -of sdf -o "{}" -od "{}" --overWrite'.format(
            os.path.abspath(ligandFile),
            os.path.splitext(os.path.basename(outputFile))[0],
            os.path.dirname(outputFile)
        )

        pwchemPlugin.runScript(
            self,
            'obabel_IO.py',
            args,
            env=OPENBABEL_DIC,
            cwd=os.path.dirname(os.path.abspath(ligandFile))
        )

        if not os.path.isfile(outputFile) or os.path.getsize(outputFile) == 0:
            raise RuntimeError(
                'Could not convert CIF ligand to SDF: {}'.format(ligandFile)
            )

        return outputFile