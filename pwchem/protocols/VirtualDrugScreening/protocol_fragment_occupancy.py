# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Joaqiun Algorta (joaquin.algorta@cnb.csic.es)
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
Score docked poses by how much of a set of reference fragments' occupancy they reproduce.
"""
import csv, json

from pwem.convert import cifToPdb
from pyworkflow.object import Float, String
from pyworkflow.protocol import params
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSmallMolecules
from pwchem.utils import os
from pwchem import Plugin, RDKIT_DIC

AGG = ['union', 'mean', 'max', 'count']
MIN_CA_PAIRS = 30


class ProtChemFragmentOccupancy(EMProtocol):
    """
    AI Generated:

This protocol scores docked ligand poses by how much of the occupancy of a set
of reference molecules (typically crystallographic fragments) they reproduce,
and optionally how many of the reference interactions they recover.

The scores are written back as new attributes of the input molecules, so they
can be used as additional voters in the Rank docking score or consensus
protocols.

Overview
--------
The protocol takes a SetOfSmallMolecules with the docked poses and a second
SetOfSmallMolecules with the reference fragments. Each reference defines a
region of the pocket (a hotspot) that good binders are expected to fill.

The fragments usually come from different crystal structures, each in its own
coordinate frame, so they are first carried into the docking frame.

Occupancy Score
---------------
For a pose P and a reference fragment F:

    coverage(P, F) = fraction of the F heavy atoms with a P heavy atom
                     within the cutoff distance

This is an asymmetric containment measure normalised by the fragment: a larger
molecule that engulfs a small fragment scores ~1.

Per pose, the coverages are aggregated over the references:

- union: fraction of the pooled reference-atom cloud covered.
  Rewards spanning the consensus hotspot, not just filling one fragment.
- mean: mean of the per-fragment coverages.
- max: best single-fragment coverage.
- count: number of fragments covered above a threshold.

Interaction Score (optional)
----------------------------
Interaction fingerprints are computed for every reference and every pose
against the docking receptor, using only elements (no hydrogens needed):

- polar: N/O - N/O contacts within 3.5 A (H-bond / salt bridge proxy)
- hydrophobic: C/S/halogen - C contacts within 4.0 A

The score is the fraction of the reference interactions (per receptor residue
and type) the pose reproduces, aggregated with the same mode as the occupancy.

Frame Alignment
---------------
Unless the references are already in the docking frame:

- Each reference protein is superposed onto every chain of the docking
  receptor by Calpha Kabsch, matching residues by chain and residue number
- The same rigid transform is applied to the fragment coordinates
- For multimeric receptors, each fragment is placed in every equivalent
  pocket and the best one is kept
- Per-fragment RMSDs are reported in extra/alignment_report.csv (sub-A
  values are expected for the same protein)

File Handling
-------------
- CIF structures are converted to PDB
- Maestro files (.mae/.maegz, Glide output) are converted to PDB with RDKit
- PDB, PDBQT, SDF, MOL and MOL2 are read directly

Workflow
--------
1. Collect the pose, reference, reference protein and receptor files
2. Convert CIF and Maestro files to PDB
3. Align each reference protein onto the receptor chains and place the fragments
4. Compute per-fragment coverages (and interaction recoveries) for every pose
5. Aggregate them per pose
6. Assign the scores to the output molecules

Output
------
- outputSmallMolecules:
    A copy of the input set enriched with:

    - _fragOcc:
        Aggregated occupancy score. Higher values indicate poses that
        reproduce more of the reference occupancy.
    - _fragOccInt (if the interaction score is enabled):
        Aggregated interaction recovery. Higher is better.
    - _fragOcc_<refName>, _fragOccInt_<refName>, _fragOccBest (if the
      per-reference columns are stored):
        Per-reference coverages / recoveries and the best covered reference.

Key Features
------------
- Fragment-based rescoring independent of docking scoring functions
- Automatic superposition of references coming from different structures
- Support for multimeric receptors with equivalent pockets
- Optional hydrogen-free interaction recovery score
- Per-reference columns to see which fragment each pose reproduces

Use Cases
---------
- Rescoring docking poses with fragment-screening (e.g. XChem) data
- Prioritizing hits that grow from or link known fragments
- Adding a structure-based voter to consensus ranking

Notes
-----
- A strict 1.0 A cutoff was validated as the most discriminant; cutoffs
  of 2 A or more saturate and lose discrimination
- Reference proteins must share residue numbering with the docking receptor
- Large natural ligands (e.g. cofactors) can be kept as columns but excluded
  from the aggregate
"""
    _label = 'fragment occupancy scoring'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label='Input')
        iGroup = form.addGroup('Input')
        iGroup.addParam('inputSmallMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Input docked molecules: ', help='Set of docked poses to be scored.')
        iGroup.addParam('referenceMolecules', params.PointerParam, pointerClass='SetOfSmallMolecules',
                        label='Reference molecules (fragments): ',
                        help='Fragments whose occupancy defines the hotspot to reproduce.')
        iGroup.addParam('excludeReferences', params.StringParam, default='', expertLevel=params.LEVEL_ADVANCED,
                        label='Exclude fragments from aggregate: ',
                        help='Fragments (comma-separated reference molecule names) that are scored and kept as per-reference '
                             'columns, but left out of the aggregated score. Useful for large natural ligands or '
                             'cofactors (e.g. NAD, XMP) that would dominate the hotspot. Use the wizard to select '
                             'one or several of the reference molecules.')
        iGroup.addParam('includeHydrogens', params.BooleanParam, default=False, expertLevel=params.LEVEL_ADVANCED,
                        label='Include hydrogens: ',
                        help='Use the hydrogen atoms (of poses and references) for the occupancy. By default only '
                             'heavy atoms are used, since hydrogens depend on the protonation and preparation of '
                             'each structure. The interaction score always uses heavy atoms.')

        group = form.addGroup('Occupancy metric')
        group.addParam('aggregation', params.EnumParam, choices=AGG, default=0,
                       label='Aggregation across references: ',
                       help='union: fraction of the pooled reference-atom cloud covered (default, most '
                            'robust). mean/max: over per-fragment coverages. '
                            'count: number of fragments covered >= threshold.')
        group.addParam('cutoff', params.FloatParam, default=1.0, label='Heavy-atom cutoff (A): ',
                       help='A reference atom is "covered" if a pose heavy atom is within this distance. '
                            '1.0 A validated as best; >= 2 A saturates.')
        group.addParam('coverageThreshold', params.FloatParam, default=0.5, condition='aggregation==3',
                       label='Per-fragment coverage threshold: ', help='Only for the "count" aggregation.')
        group.addParam('perReferenceColumns', params.BooleanParam, default=True,
                       label='Store per-reference coverage: ',
                       help='Also store each reference coverage as _fragOcc_<refName> and the name of '
                            'the best covered aggregated reference as _fragOccBest.')

        group.addParam('interactionScore', params.BooleanParam, default=False,
                       label='Also score interaction recovery: ',
                       help='Also compute the fraction of the reference interactions with the docking receptor '
                            'that each pose reproduces, stored as _fragOccInt (and _fragOccInt_<refName> '
                            'with the per-reference columns). Interactions are element-based, so no hydrogens '
                            'are needed: polar (N/O - N/O <= 3.5 A) and hydrophobic (C/S/halogen - C <= 4.0 A) '
                            'contacts per receptor residue. Uses the same aggregation as the occupancy.')

        group = form.addGroup('Frame alignment')
        group.addParam('assumeSameFrame', params.BooleanParam, default=False,
                       label='References already in the docking frame? ',
                       help='If Yes, use reference coordinates directly. If No, each reference protein is '
                            'Calpha-superposed onto the docking receptor and the transform applied to '
                            'the fragment.')
        group.addParam('placeInAllChains', params.BooleanParam, default=True,
                       condition='assumeSameFrame==False', expertLevel=params.LEVEL_ADVANCED,
                       label='Place fragments in all receptor chains: ',
                       help='For multimeric receptors; take the max coverage over equivalent pockets.')

        form.addParallelSection(threads=4, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.scoreStep)
        self._insertFunctionStep(self.createOutputStep)

    def scoreStep(self):
        """Write the manifest the geometry script consumes and run it in the RDKit env."""
        inSet, refSet = self.inputSmallMolecules.get(), self.referenceMolecules.get()
        poses = [{'objId': mol.getObjId(), 'molName': mol.getMolName(),
                  'poseFile': self.getPdbFile(mol.getPoseFile() or mol.getFileName())} for mol in inSet]
        refs = [{'molName': mol.getMolName(), 'molFile': self.getPdbFile(mol.getPoseFile() or mol.getFileName()),
                 'proteinFile': self.getPdbFile(mol.getProteinFile() or refSet.getProteinFile())} for mol in refSet]

        manifest = {
            'poses': poses, 'refs': refs, 'receptor': self.getPdbFile(inSet.getProteinFile()),
            'cutoff': self.cutoff.get(), 'aggregation': AGG[self.aggregation.get()],
            'coverageThreshold': self.coverageThreshold.get(),
            'assumeSameFrame': self.assumeSameFrame.get(), 'placeInAllChains': self.placeInAllChains.get(),
            'includeHydrogens': self.includeHydrogens.get(), 'exclude': self.getExcluded(),
            'minCaPairs': MIN_CA_PAIRS,
            'interactions': self.interactionScore.get(),
            'outCsv': os.path.abspath(self.getScoresFile()),
            'alignReport': os.path.abspath(self.getAlignFile()),
        }
        if os.path.exists(self.getMaeDir()):
            # Maestro files (Glide output), staged by getPdbFile, to PDB with RDKit in a single call
            args = f' --multiFiles -iD "{self.getMaeDir()}" --pattern "*" -of pdb --keepHs --outputDir "{self.getMaeDir()}"'
            Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=self.getMaeDir())

        manifestFile = os.path.abspath(self._getExtraPath('manifest.json'))
        with open(manifestFile, 'w') as f:
            json.dump(manifest, f)

        Plugin.runScript(self, 'fragment_occupancy_score.py', manifestFile, env=RDKIT_DIC,
                         cwd=self._getExtraPath())

    def createOutputStep(self):
        scores = {}
        with open(self.getScoresFile()) as f:
            for row in csv.DictReader(f):
                scores[int(row.pop('objId'))] = row

        attr = '_fragOcc'
        inSet = self.inputSmallMolecules.get()
        outSet = SetOfSmallMolecules.createCopy(inSet, self._getPath(), copyInfo=True)
        for mol in inSet:
            newMol = mol.clone()
            row = scores.get(mol.getObjId(), {})
            setattr(newMol, attr, Float(row.get('fragOcc') or None))
            if self.interactionScore.get():
                setattr(newMol, attr + 'Int', Float(row.get('fragInt') or None))
            if self.perReferenceColumns.get():
                setattr(newMol, attr + 'Best', String(row.get('bestFrag') or None))
                for col, val in row.items():
                    if col.startswith('cov_'):
                        setattr(newMol, attr + '_' + col[4:], Float(val or None))
                    elif col.startswith('int_'):
                        setattr(newMol, attr + 'Int_' + col[4:], Float(val or None))
            outSet.append(newMol)

        self._defineOutputs(outputSmallMolecules=outSet)
        self._defineSourceRelation(self.inputSmallMolecules, outSet)
        self._defineSourceRelation(self.referenceMolecules, outSet)

    # --------------------------- INFO functions -----------------------------------
    def _validate(self):
        errors = []
        if not self.inputSmallMolecules.get().getProteinFile() and \
                (not self.assumeSameFrame.get() or self.interactionScore.get()):
            errors.append('The input molecules have no receptor (needed for the alignment and interactions).')
        if not self.assumeSameFrame.get():
            refSet = self.referenceMolecules.get()
            noProt = [mol.getMolName() for mol in refSet if not (mol.getProteinFile() or refSet.getProteinFile())]
            if noProt:
                errors.append('These reference molecules have no protein file to align: {}. '
                              'Set "References already in the docking frame?" if they need no '
                              'alignment.'.format(', '.join(noProt)))
        return errors

    def _summary(self):
        summary = ['Fragment occupancy scoring ({} aggregation, cutoff {:.1f} A) -> attribute _fragOcc.'.
                   format(AGG[self.aggregation.get()], self.cutoff.get())]
        if os.path.exists(self.getAlignFile()):
            with open(self.getAlignFile()) as f:
                rows = list(csv.DictReader(f))
            bad = [r['fragment'] for r in rows if int(r['nCalpha']) < MIN_CA_PAIRS]
            summary.append('Alignment RMSD range: {:.2f}-{:.2f} A (see extra/alignment_report.csv).'.
                           format(*[f(float(r['rmsd']) for r in rows) for f in (min, max)]) if rows else
                           'No reference could be aligned to the receptor.')
            if bad:
                summary.append('WARNING: fewer than {} matched Calpha for: {}'.
                               format(MIN_CA_PAIRS, ', '.join(sorted(set(bad)))))
        return summary

    def _methods(self):
        return ['Each docked pose was scored by the fraction of reference-fragment heavy-atom occupancy '
                'reproduced within {:.1f} A ({} aggregation), after Calpha superposition of each reference '
                'structure onto the docking receptor.'.format(self.cutoff.get(), AGG[self.aggregation.get()])]

    # --------------------------- UTILS functions -----------------------------------
    def getScoresFile(self):
        return self._getExtraPath('occupancy_scores.csv')

    def getAlignFile(self):
        return self._getExtraPath('alignment_report.csv')

    def getExcluded(self):
        return [x.strip() for x in self.excludeReferences.get().split(',') if x.strip()]

    def getMaeDir(self):
        return os.path.abspath(self._getTmpPath('maestro'))

    def getPdbFile(self, inFile):
        """The geometry script cannot read CIF or Maestro (proteins or molecules). CIF is converted to PDB here,
        Maestro is staged in getMaeDir and converted with RDKit by scoreStep. Returns the PDB to use."""
        if not inFile:
            return None
        inFile = os.path.abspath(inFile)
        base, ext = os.path.splitext(inFile)
        # Named after the full path: different inputs may share a basename
        uniqueName = base.strip('/').replace('/', '_')
        if ext.lower() == '.cif':
            pdbFile = os.path.abspath(self._getTmpPath(uniqueName + '.pdb'))
            if not os.path.exists(pdbFile):
                cifToPdb(inFile, pdbFile)
            return pdbFile
        if ext.lower() in ('.mae', '.maegz'):
            os.makedirs(self.getMaeDir(), exist_ok=True)
            maeFile = os.path.join(self.getMaeDir(), uniqueName + ext)
            if not os.path.exists(maeFile):
                os.symlink(inFile, maeFile)
            return os.path.join(self.getMaeDir(), uniqueName + '.pdb')
        return inFile
