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

from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import pdbqt2other, getBaseName, parseAtomStruct, \
    calculateDistance, convertToSdf, getChainIds
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import RDKIT_DIC


class ProtContactLigands(EMProtocol):
    """
            Filters docked small molecules based on the distance between a selected
            ligand atom and a selected receptor atom.

            The protocol evaluates a specific atom-to-atom contact for each docked
            ligand in an input SetOfSmallMolecules. The ligand atom is identified
            by its 1-based position among the heavy atoms of a reference ligand.
            RDKit is used to map this atom onto the corresponding atom in each
            docked pose, allowing the protocol to retrieve its coordinates even
            when the atom ordering differs between the reference molecule and the
            docked structures.

            Molecules are retained when the distance between the mapped ligand atom
            and the selected receptor atom is equal to or below the specified
            cutoff. The evaluation results are also recorded in a CSV report.

            Workflow
            --------
            1. Receive a SetOfSmallMolecules containing docked ligand poses and
               their associated receptor structure.
            2. Retrieve the receptor structure and validate the filtering
               parameters.
            3. Select the first ligand pose as the reference for ligand atom
               mapping.
            4. Load the receptor structure, converting it from PDBQT to PDB when
               necessary.
            5. Identify the selected receptor chain and residue.
            6. Retrieve the coordinates of the selected receptor atom.
            7. Iterate over the input docked ligands.
            8. Convert the reference ligand and each docked pose to SDF format
               when necessary.
            9. Use RDKit to map the selected reference ligand heavy atom onto the
               corresponding atom in each docked pose.
            10. Retrieve the mapped ligand atom coordinates and calculate the
                Euclidean distance to the selected receptor atom.
            11. Retain ligands whose atom-to-atom distance is equal to or below
                the specified cutoff.
            12. Write the evaluation results for all ligands to a CSV report.
            13. Create an output SetOfSmallMolecules containing the retained
                docked ligands.

            Input
            -----
            - inputLigandSet:
                SetOfSmallMolecules containing the docked ligand poses to
                evaluate.

                The input set must be marked as docked and must contain an
                accessible receptor structure associated with the docking results.

            Parameters
            ----------
            - Receptor chain:
                Chain identifier of the receptor containing the residue to
                evaluate.

                The chain is selected using the chain-selection wizard.

            - Receptor residue:
                Residue in the selected receptor chain containing the atom to
                evaluate.

                The residue is selected using the residue-selection wizard. Its
                residue index and, when available, residue name are used to
                identify it in the receptor structure.

            - Receptor atom:
                Name of the receptor atom whose coordinates will be used for
                contact evaluation.

                The atom name must correspond to the name displayed in PyMOL,
                such as ``OG``, ``CA`` or ``OD1``. It must be an actual atom name,
                not a sequential atom label.

            - Ligand atom:
                1-based index of the selected ligand heavy atom in the reference
                ligand.

                Hydrogen and deuterium atoms are excluded from the indexing.
                The selected atom is mapped onto the corresponding atom in each
                docked pose using RDKit.

                The default value is 0, but a valid index starting from 1 must be
                provided.

            - Distance cutoff:
                Maximum allowed distance, in Å, between the selected ligand atom
                and the selected receptor atom.

                A ligand is retained when its calculated distance is equal to or
                below this cutoff.

                The default value is 4 Å.

            Atom Mapping and Distance Evaluation
            ------------------------------------
            The protocol uses the first molecule in the input set as the reference
            for identifying the selected ligand atom.

            Before mapping, the reference ligand and each docked pose are converted
            to SDF format using the available molecular conversion utilities.

            The RDKit atom-mapping script identifies the atom corresponding to the
            selected reference heavy-atom index in each docked pose and returns:

            - The mapped atom coordinates.
            - The mapped atom element.
            - The mapped atom index in the pose.

            The protocol then calculates the Euclidean distance between the mapped
            ligand atom and the selected receptor atom.

            The filtering criterion is:

            - Retained: distance <= cutoff.
            - Discarded: distance > cutoff.

            Molecules that cannot be evaluated, for example because atom mapping
            fails or the required pose file is unavailable, are discarded and
            reported with the corresponding error in the CSV file.

            Output
            ------
            - outputSmallMolecules:
                SetOfSmallMolecules containing the docked ligands that satisfy the
                specified atom-to-atom distance criterion.

                The output set preserves the input set information and associated
                receptor structure. It is marked as docked, and each retained
                molecule is copied from the input set.

            - contact_distances.csv:
                CSV report containing the evaluation results for every input
                molecule, including:

                - Molecule name
                - Docked pose file
                - Reference ligand file
                - Receptor chain
                - Receptor residue
                - Receptor atom
                - Reference ligand heavy-atom index
                - Mapped atom index in the docked pose
                - Mapped ligand atom element
                - Ligand atom coordinates
                - Receptor atom coordinates
                - Calculated distance
                - Distance cutoff
                - Retention status
                - Error message, when applicable

            Summary
            -------
            The protocol filters docked ligands according to a selected
            ligand–receptor atom contact.

            The output contains only the molecules whose mapped ligand atom is
            within the specified distance cutoff of the selected receptor atom.
            Detailed evaluation results are available in
            ``contact_distances.csv``.

            Use Cases
            ---------
            - Filtering docked ligands based on a specific protein–ligand contact
            - Identifying poses that preserve interactions with selected
              active-site residues
            - Prioritizing docking poses according to a known interaction
              constraint
            - Selecting ligands that position a specific functional group near a
              target receptor atom
            - Applying atom-level geometric constraints to docking results
            - Preparing filtered ligand collections for downstream analysis

            Notes
            -----
            The input must be a docked SetOfSmallMolecules with an accessible
            associated receptor structure.

            The first molecule in the input set is used as the reference ligand
            for atom mapping. The selected ligand atom index therefore refers to
            this reference molecule and must identify the same chemical atom in
            the other docked poses.

            Atom mapping is performed using RDKit and does not rely on the
            sequential atom ordering being identical between the reference ligand
            and the docked poses.

            The ligand atom index is 1-based and excludes hydrogen and deuterium
            atoms. The receptor atom is identified by its atom name within the
            selected residue.

            The selected receptor atom must be a unique heavy atom in the selected
            residue. Hydrogen and deuterium atoms cannot be used for the distance
            evaluation.

            The protocol evaluates one ligand atom against one receptor atom.
            It does not calculate the minimum distance between the complete
            ligand and receptor, nor does it evaluate all possible interactions
            between the two molecules.

            The distance cutoff must be positive. A ligand is retained when its
            calculated distance is exactly equal to the cutoff as well as when it
            is shorter.

            Individual ligand evaluation errors do not stop the entire filtering
            procedure. Molecules that cannot be evaluated are discarded, and
            their errors are recorded in the CSV report.
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
                      help="Select the receptor residue to evaluate.")
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
        self._insertFunctionStep(self.filterLigandsStep)
        self._insertFunctionStep(self.createOutputStep)

    def filterLigandsStep(self):
        """Calculate atom-to-atom distances and retain passing ligands."""

        inputSet, referenceFile, cutoff, ligandAtomIndex = (
            self._prepareFilteringInputs()
        )

        receptorCoordinates = self._prepareReceptorCoordinates(inputSet)

        self._passingMolecules = []
        reportFile = self._getPath('contact_distances.csv')

        with open(reportFile, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            self._writeReportHeader(writer)

            for molecule in inputSet:
                result = self._evaluateMolecule(
                    molecule,
                    referenceFile,
                    receptorCoordinates,
                    cutoff,
                    ligandAtomIndex
                )
                self._writeMoleculeResult(
                    writer,
                    molecule,
                    result,
                    referenceFile,
                    receptorCoordinates,
                    cutoff,
                    ligandAtomIndex
                )

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
    def _prepareFilteringInputs(self):
        """Validate the input set and filtering parameters."""

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

        referenceFile = self._getPoseFile(next(iter(inputSet)))

        self.info(
            'Reference ligand for atom mapping: {}'.format(referenceFile)
        )

        return inputSet, referenceFile, cutoff, ligandAtomIndex

    def _prepareReceptorCoordinates(self, inputSet):
        """Load the receptor and retrieve the selected atom coordinates."""

        receptorFile = self.convertFormat(inputSet.getProteinFile())
        receptor = parseAtomStruct(receptorFile)
        targetResidue = self._getReceptorResidue(receptor)

        return self._getReceptorAtomCoordinates(targetResidue)

    def _evaluateMolecule(
            self,
            molecule,
            referenceFile,
            receptorCoordinates,
            cutoff,
            ligandAtomIndex
    ):
        """Map the reference ligand atom and evaluate its distance."""

        result = {
            'poseFile': '',
            'ligandCoordinates': None,
            'ligandElement': '',
            'poseAtomIndex': '',
            'distance': None,
            'error': '',
            'keep': False
        }

        try:
            poseFile = self._getPoseFile(molecule)
            result['poseFile'] = poseFile

            (
                result['ligandCoordinates'],
                result['ligandElement'],
                result['poseAtomIndex']
            ) = self._getMappedLigandAtomCoordinates(
                referenceFile,
                poseFile
            )

            result['distance'] = calculateDistance(
                result['ligandCoordinates'],
                receptorCoordinates
            )

            result['keep'] = result['distance'] <= cutoff

            if result['keep']:
                self._passingMolecules.append(molecule.clone())

            self._logMoleculeEvaluation(
                poseFile,
                result,
                receptorCoordinates,
                ligandAtomIndex
            )

        except Exception as exc:
            result['error'] = str(exc)
            result['keep'] = False

            self.warning(
                'Could not evaluate {}: {}'.format(
                    molecule.getMolName(), result['error']
                )
            )

        return result

    def _logMoleculeEvaluation(
            self,
            poseFile,
            result,
            receptorCoordinates,
            ligandAtomIndex
    ):
        """Print the distance evaluation for a molecule."""

        print(
            '{} | reference heavy atom {} | '
            'pose atom {} ({}) at {} | '
            'receptor atom {} at {} | '
            'distance={:.3f} A | retained={}'.format(
                os.path.basename(poseFile),
                ligandAtomIndex,
                result['poseAtomIndex'],
                result['ligandElement'],
                np.round(result['ligandCoordinates'], 3).tolist(),
                self.recAtom.get(),
                np.round(receptorCoordinates, 3).tolist(),
                result['distance'],
                result['keep']
            )
        )

    def _writeReportHeader(self, writer):
        """Write the contact-distance report header."""

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

    def _writeMoleculeResult(
            self,
            writer,
            molecule,
            result,
            referenceFile,
            receptorCoordinates,
            cutoff,
            ligandAtomIndex
    ):
        """Write one molecule's evaluation to the CSV report."""

        ligandCoordinates = result['ligandCoordinates']
        distance = result['distance']

        writer.writerow([
            molecule.getMolName(),
            result['poseFile'],
            referenceFile,
            self.recChain.get(),
            self.recResidue.get(),
            self.recAtom.get(),
            ligandAtomIndex,
            result['poseAtomIndex'],
            result['ligandElement'],
            (
                np.round(ligandCoordinates, 3).tolist()
                if ligandCoordinates is not None else ''
            ),
            np.round(receptorCoordinates, 3).tolist(),
            '{:.3f}'.format(distance) if distance is not None else '',
            cutoff,
            'Yes' if result['keep'] else 'No',
            result['error']
        ])

    def _getReceptorResidue(self, structure):
        """Find the selected receptor residue using the wizard outputs."""

        try:
            chainId = getChainIds(self.recChain.get())[0]
            residueData = json.loads(self.recResidue.get())

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
                heteroFlag, resseq, _ = residue.id

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

    def convertFormat(self, filename):
        """Convert PDBQT receptor files to PDB; other formats are read directly."""
        if filename.lower().endswith('.pdbqt'):
            outFile = os.path.abspath(self._getExtraPath(getBaseName(filename) + '.pdb'))
            return pdbqt2other(self, filename, outFile)
        return filename

    def _getMappedLigandAtomCoordinates(self, referenceFile, poseFile):
        referenceFile = convertToSdf(self, referenceFile)
        poseFile = convertToSdf(self, poseFile)

        # The script writes 'mapped_ligand_atom.json' into its cwd.
        outputDirectory = os.path.abspath(self._getTmpPath())
        os.makedirs(outputDirectory, exist_ok=True)
        outputFile = os.path.join(outputDirectory, 'mapped_ligand_atom.json')

        # Remove any result from a previous run so it is never reused.
        if os.path.isfile(outputFile):
            os.remove(outputFile)

        args = (
            '--reference "{}" '
            '--pose "{}" '
            '--atom-index {}'
        ).format(
            os.path.abspath(referenceFile),
            os.path.abspath(poseFile),
            int(self.ligAtom.get())
        )

        pwchemPlugin.runScript(
            self,
            'rdkit_map_ligand_atom.py',
            args,
            env=RDKIT_DIC,
            popen=True,
            cwd=outputDirectory
        )

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