# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Judith Maestro (scipionchem@cnb.csic.es)
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
This protocol gets a residue/atom contact map computed by COCADA (https://github.com/LBS-UFMG/COCaDA),
either by importing a CSV already computed (e.g. downloaded from the COCADA web server,
https://bioinfo.dcc.ufmg.br/cocada-web/public/) or by running COCADA itself on an input AtomStruct,
and defines one structural ROI per interaction type found in the file (hydrogen bonds, hydrophobic,
salt bridges...), plus one extra ROI grouping the whole interface, so the different bond types can
be inspected / visualized independently. COCADA computes intra- and inter-molecular contacts on any
kind of biomolecular structure, not just two-chain protein complexes: single-chain proteins,
multi-chain complexes, protein-peptide, protein-DNA and protein-RNA structures are all supported
(see the example entries in the COCADA documentation).
"""
import os, csv, shutil, glob, json
from collections import defaultdict

import numpy as np

from pyworkflow.object import String
from pyworkflow.protocol import params
from pyworkflow.protocol.constants import LEVEL_ADVANCED
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem import Plugin
from pwchem.objects import SetOfStructROIs, StructROI
from pwchem.utils import parseAtomStruct, cifFromASFile, createPocketFile, natural_sort, \
    createColorVectors, RESIDUES1TO3

ALL_LABEL = 'All'
RUN_MODE_IMPORT, RUN_MODE_COMPUTE = 0, 1
COND_IMPORT = 'cocadaMode==%d' % RUN_MODE_IMPORT
COND_COMPUTE = 'cocadaMode==%d' % RUN_MODE_COMPUTE

# COCADA interaction type codes (https://bioinfo.dcc.ufmg.br/cocada-web/public/documentation/),
# only used for the ROI comment. Fall back to the raw code if one is not present here.
COCADA_TYPE_NAMES = {
    'HB': 'Hydrogen Bond', 'HY': 'Hydrophobic', 'SB': 'Salt Bridge',
    'AT': 'Attractive', 'RE': 'Repulsive', 'AS': 'Aromatic Stacking', 'DS': 'Disulfide Bond',
    'uAT': 'Uncertain Attractive', 'uRE': 'Uncertain Repulsive', 'uSB': 'Uncertain Salt Bridge',
}

# Side-chain atoms defining the aromatic ring(s) COCADA collapses into its "RNG" pseudo-atom.
# Trp's indole is approximated by its 6-membered ring; COCADA's own convention is not documented.
RING_ATOMS = {
    'PHE': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'TYR': ['CG', 'CD1', 'CD2', 'CE1', 'CE2', 'CZ'],
    'HIS': ['CG', 'ND1', 'CD2', 'CE1', 'NE2'],
    'TRP': ['CD2', 'CE2', 'CE3', 'CZ2', 'CZ3', 'CH2'],
}


class ProtCocadaInteractions(EMProtocol):
    """
    Gets a COCADA residue-residue contacts CSV (columns: Chain1,Res1,ResName1,Atom1,Chain2,
    Res2,ResName2,Atom2,Distance,Type) -- either imported or freshly computed by running COCADA
    on the input structure -- and defines a SetOfStructROIs with one ROI per interaction type
    plus one ROI for the whole interface. The contacts can come from any structure supported by
    COCADA (single-chain, multi-chain complex, protein-peptide, protein-DNA, protein-RNA...), not
    only a two-chain protein-protein complex.
    """
    _label = 'COCADA interactions'

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('inputAtomStruct', params.PointerParam, pointerClass='AtomStruct',
                       label='Input AtomStruct: ', allowsNull=False,
                       help='Structure the COCADA contacts are/were calculated on (single chain, '
                            'multi-chain complex, protein-peptide, protein-DNA, protein-RNA...). '
                            'If importing a CSV, its chain and residue numbering must match this '
                            'structure.')
        group.addParam('cocadaMode', params.EnumParam, choices=['Import existing CSV', 'Compute from structure (run COCADA)'],
                       default=RUN_MODE_COMPUTE, display=params.EnumParam.DISPLAY_HLIST,
                       label='Contacts source: ',
                       help='Either import a COCADA contacts CSV already computed elsewhere (e.g. '
                            'downloaded from the COCADA web server, '
                            'https://bioinfo.dcc.ufmg.br/cocada-web/public/), or run COCADA '
                            '(https://github.com/LBS-UFMG/COCaDA) on the input AtomStruct to '
                            'generate it.')
        group.addParam('cocadaCsv', params.PathParam, label='COCADA contacts file: ', allowsNull=False,
                       condition=COND_IMPORT,
                       help='CSV file downloaded from the COCADA web server '
                            '(https://bioinfo.dcc.ufmg.br/cocada-web/public/) with the columns: '
                            'Chain1,Res1,ResName1,Atom1,Chain2,Res2,ResName2,Atom2,Distance,Type')

        group = form.addGroup('COCADA parameters', condition=COND_COMPUTE)
        group.addParam('phValue', params.FloatParam, default=7.4,
                       label='pH value: ',
                       help='pH used to define electrostatic (attractive/repulsive/salt bridge) '
                            'contacts (0.1-14). Default is 7.4 (COCADA\'s own default). Maps to '
                            "COCADA's -ph flag.")
        group.addParam('interchainOnly', params.BooleanParam, default=False,
                       label='Interchain contacts only? ',
                       help="Only calculate contacts between different chains. Maps to COCADA's "
                            "-inter flag.")
        group.addParam('chainsParam', params.StringParam, allowsNull=True,
                       label='Chains: ',
                       help='Restrict the analysis to these chains, e.g. "A,B". Leave empty to '
                            "analyze all chains. Use the wizard or type the chain letters "
                            "directly. Maps to COCADA's -c flag.")
        group.addParam('regionParam', params.StringParam, allowsNull=True, expertLevel=LEVEL_ADVANCED,
                       label='Residue region: ',
                       help='Restrict the analysis to a region of residues, either a range '
                            '("10-19") or a list ("10,32,65"). Leave empty to analyze all '
                            "residues. Maps to COCADA's -r flag.")

        distGroup = form.addGroup('Custom contact distances', condition=COND_COMPUTE,
                                  expertLevel=LEVEL_ADVANCED,
                                  help='Distance ranges (in Angstroms) defining each contact type, '
                                       "prefilled with COCADA's own defaults. Change them if needed. "
                                       "Maps to COCADA's -d flag.")
        line = distGroup.addLine('Salt bridge: ', help='Equally charged atoms AND hydrogen bonding')
        line.addParam('saltBridgeMin', params.FloatParam, label='Min', default=0)
        line.addParam('saltBridgeMax', params.FloatParam, label='Max', default=3.9)

        line = distGroup.addLine('Hydrophobic: ', help='Hydrophobic atom pair')
        line.addParam('hydrophobicMin', params.FloatParam, label='Min', default=2.0)
        line.addParam('hydrophobicMax', params.FloatParam, label='Max', default=4.5)

        line = distGroup.addLine('Hydrogen bond: ', help='Acceptor and Donor atom pair')
        line.addParam('hydrogenBondMin', params.FloatParam, label='Min', default=0)
        line.addParam('hydrogenBondMax', params.FloatParam, label='Max', default=3.9)

        line = distGroup.addLine('Repulsive: ', help='Equally charged atoms')
        line.addParam('repulsiveMin', params.FloatParam, label='Min', default=2.0)
        line.addParam('repulsiveMax', params.FloatParam, label='Max', default=6.0)

        line = distGroup.addLine('Attractive: ', help='Differently charged atoms')
        line.addParam('attractiveMin', params.FloatParam, label='Min', default=3.9)
        line.addParam('attractiveMax', params.FloatParam, label='Max', default=6.0)

        line = distGroup.addLine('Disulfide bond: ', help='Cys:SG atom pair')
        line.addParam('disulfideBondMin', params.FloatParam, label='Min', default=0)
        line.addParam('disulfideBondMax', params.FloatParam, label='Max', default=2.8)

        line = distGroup.addLine('Aromatic stacking: ',
                                 help='Centroids of two aromatic rings in parallel or '
                                      'perpendicular orientation')
        line.addParam('aromaticMin', params.FloatParam, label='Min', default=2.0)
        line.addParam('aromaticMax', params.FloatParam, label='Max', default=5.0)

        form.addParallelSection(threads=1, mpi=1)

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep(self.convertInputStep)
        if self.cocadaMode.get() == RUN_MODE_COMPUTE:
            self._insertFunctionStep(self.runCocadaStep)
        self._insertFunctionStep(self.defineOutputStep)

    def convertInputStep(self):
        inFile = os.path.abspath(self.inputAtomStruct.get().getFileName())
        cifFromASFile(inFile, self.getCifFile())

    def runCocadaStep(self):
        outDir = self.getCocadaOutDir()
        os.makedirs(outDir, exist_ok=True)

        args = '-f {} -o {}'.format(self.getCifFile(), outDir)
        if self.numberOfThreads.get() > 1:
            args += ' -m {}'.format(self.numberOfThreads.get())
        if self.interchainOnly.get():
            args += ' -inter'
        if self.chainsParam.get():
            args += ' -c {}'.format(self.getChainsArg())
        if self.regionParam.get():
            args += ' -r {}'.format(self.regionParam.get().strip())
        if self.phValue.get() is not None:
            args += ' -ph {}'.format(self.phValue.get())
        args += ' -d {}'.format(self.getDistancesArg())

        Plugin.runCocada(self, args)

        if not self.getComputedCsvPath():
            raise RuntimeError('COCADA did not produce a contacts CSV in {}'.format(outDir))

    def defineOutputStep(self):
        self.structModel = parseAtomStruct(self.getCifFile())[0]
        rowsByType = self.parseCocadaCsv()

        groups = list(rowsByType.items())
        if len(groups) > 1:
            allRows = [row for rows in rowsByType.values() for row in rows]
            groups.append((ALL_LABEL, allRows))

        outPockets = SetOfStructROIs(filename=self._getPath('StructROIs.sqlite'))
        for i, (bondType, rows) in enumerate(groups):
            pocket = self.buildROIFromRows(rows, i + 1, bondType)
            if pocket is not None:
                outPockets.append(pocket)

        if len(outPockets) > 0:
            outPockets.buildPDBhetatmFile()

            csvPath = self.getCocadaCsvPath()
            rawCsv = self._getExtraPath(os.path.basename(csvPath))
            shutil.copyfile(csvPath, rawCsv)

            summaryCsv = self.buildResiduePairsCsv(rowsByType)
            outPockets.setInteractingResiduesFile(summaryCsv)

            outPockets._cocadaContactsFile = String(os.path.relpath(rawCsv))

            bondsPml = self.buildBondsPml(rowsByType)
            outPockets._cocadaBondsPml = String(os.path.relpath(bondsPml))

            self._defineOutputs(outputStructROIs=outPockets)
            self._defineSourceRelation(self.inputAtomStruct, outPockets)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        if self.getOutputsSize() > 0:
            summary.append('Got {} structural ROIs from the COCADA contacts, one per '
                           'interaction type found ({}).'.
                           format(self.outputStructROIs.getSize(), self.getCocadaCsvPath()))
        return summary

    def _validate(self):
        errors = []
        if self.cocadaMode.get() == RUN_MODE_IMPORT:
            if not os.path.exists(self.cocadaCsv.get()):
                errors.append('The COCADA contacts file does not exist: {}'.format(self.cocadaCsv.get()))
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def getCifFile(self):
        '''Default path of the cif file converted from the input AtomStruct in convertInputStep.
        Used instead of storing the path in self between steps, since plain attributes set in one
        step are not persisted if the protocol is stopped and resumed later.'''
        return os.path.abspath(self._getExtraPath('inputStructure.cif'))

    def getCocadaOutDir(self):
        '''Default directory where COCADA writes its output when computing the contacts'''
        return os.path.abspath(self._getExtraPath('cocada_output'))

    def getComputedCsvPath(self):
        '''Find the contacts CSV produced by runCocadaStep in its default output directory
        (glob because COCADA names the file after the input, not with a fixed name). Recomputed
        on demand instead of cached in self, for the same reason as getCifFile.'''
        matches = glob.glob(os.path.join(self.getCocadaOutDir(), '*_contacts.csv'))
        return matches[0] if matches else None

    def getChainsArg(self):
        '''Parse the chainsParam value into a comma-separated list of chain ids for COCADA's -c
        flag. Accepts either a plain "A,B" string typed by the user or the JSON dictionary
        produced by the chain-selection wizard (SelectMultiChainWizard), which encodes chains as
        {"chain": "A", ...} (single selection) or {"model-chain": "0-A, 0-B"} (multiple)'''
        value = self.chainsParam.get().strip()
        try:
            chainJson = json.loads(value)
        except ValueError:
            return value

        if 'chain' in chainJson:
            return chainJson['chain'].upper().strip()
        elif 'model-chain' in chainJson:
            chainIds = [mc.strip().split('-')[-1] for mc in chainJson['model-chain'].upper().split(',')]
            return ','.join(chainIds)
        return value

    def getDistancesArg(self):
        '''Build the 14 comma-separated distance values (min,max per contact type, in the order
        COCADA's own -d flag expects) from the individual min/max form params'''
        order = ['saltBridge', 'hydrophobic', 'hydrogenBond', 'repulsive', 'attractive',
                'disulfideBond', 'aromatic']
        values = []
        for name in order:
            values.append(str(getattr(self, name + 'Min').get()))
            values.append(str(getattr(self, name + 'Max').get()))
        return ','.join(values)

    def getCocadaCsvPath(self):
        '''Return the COCADA contacts csv path: the one imported by the user, or the one just
        computed by runCocadaStep, depending on the selected cocadaMode'''
        if self.cocadaMode.get() == RUN_MODE_IMPORT:
            return self.cocadaCsv.get()
        return self.getComputedCsvPath()

    def parseCocadaCsv(self):
        '''Group the rows of the COCADA csv by interaction Type'''
        rowsByType = defaultdict(list)
        with open(self.getCocadaCsvPath()) as f:
            reader = csv.DictReader(f)
            for row in reader:
                rowsByType[row['Type'].strip()].append(row)
        return rowsByType

    def buildROIFromRows(self, rows, idx, bondType):
        '''Build a StructROI from the atoms involved in a list of COCADA contact rows'''
        coords, atomIds, resIds = [], set(), set()
        for row in rows:
            for chainKey, resKey, atomKey in (('Chain1', 'Res1', 'Atom1'), ('Chain2', 'Res2', 'Atom2')):
                chainId, resNum = row[chainKey].strip(), row[resKey].strip()
                coord, serial = self.getAtomCoord(chainId, resNum, row[atomKey].strip())
                if coord is not None:
                    coords.append(coord)
                    resIds.add('{}_{}'.format(chainId, resNum))
                    if serial is not None:
                        atomIds.add(serial)

        if not coords:
            self.warning('No atoms of the input structure could be matched for interaction type {}, '
                         'skipping ROI'.format(bondType))
            return None

        pocketFile = self._getExtraPath('pocketFile_COCADA_{}.cif'.format(bondType))
        createPocketFile(coords, idx, pocketFile)

        # Build with proteinFile=None so the constructor does not auto-recalculate contacts by
        # distance; the exact COCADA contacts are set explicitly right below instead.
        pocket = StructROI(pocketFile, proteinFile=None, pClass='COCADA')
        pocket.setProteinFile(self.getCifFile())
        pocket.setContactAtoms('-'.join(natural_sort(list(atomIds))))
        pocket.setContactResidues('-'.join(natural_sort(list(resIds))))
        pocket.setVolume(pocket.getPocketVolume())

        typeName = COCADA_TYPE_NAMES.get(bondType, bondType)
        pocket.setObjLabel('COCADA_{} ({} contacts)'.format(bondType, len(rows)))
        pocket.setObjComment('{}: {} contact(s) imported from COCADA'.format(typeName, len(rows)))
        return pocket

    def getAtomCoord(self, chainId, resNum, atomName):
        '''Return (coord, serialNumber) for a COCADA atom reference. COCADA reports aromatic
        (AT/AS) contacts using a "RNG" pseudo-atom (the ring centroid) instead of a real atom,
        in which case serialNumber is None since it is not an actual structure atom.'''
        try:
            residue = self.structModel[chainId][int(resNum)]
        except KeyError:
            self.warning('Could not find residue {}:{} in the input structure'.format(chainId, resNum))
            return None, None

        if atomName == 'RNG':
            ringAtomNames = RING_ATOMS.get(residue.get_resname())
            if ringAtomNames is None:
                self.warning('Unknown aromatic ring atoms for residue {}:{} ({}), skipping RNG contact'.
                             format(chainId, resNum, residue.get_resname()))
                return None, None
            ringCoords = [residue[aName].get_coord() for aName in ringAtomNames if aName in residue]
            if not ringCoords:
                return None, None
            return list(np.mean(ringCoords, axis=0)), None

        try:
            atom = residue[atomName]
            return list(atom.get_coord()), str(atom.serial_number)
        except KeyError:
            self.warning('Could not find atom {}:{}.{} in the input structure'.
                         format(chainId, resNum, atomName))
            return None, None

    def buildResiduePairsCsv(self, rowsByType):
        '''Aggregate the per-atom COCADA rows into one row per interacting residue pair (minimum
        distance among its atom pairs), in the format expected by ViewerGeneralStructROIs'
        "Residue interaction view" (Chain1,Residue1,Chain2,Residue2,Min distance), the same one
        the PPI branch of "Define structural ROIs" produces.'''
        pairDist, pairTypes = {}, defaultdict(set)
        for bondType, rows in rowsByType.items():
            for row in rows:
                chain1, chain2 = row['Chain1'].strip(), row['Chain2'].strip()
                resRep1 = self.getResidueRepr(row['ResName1'], row['Res1'])
                resRep2 = self.getResidueRepr(row['ResName2'], row['Res2'])
                key = (chain1, resRep1, chain2, resRep2)
                dist = float(row['Distance'])
                pairDist[key] = min(dist, pairDist.get(key, dist))
                pairTypes[key].add(bondType)

        outCsv = self._getExtraPath('interacting_residues.csv')
        with open(outCsv, 'w', newline='') as f:
            # The viewer (ViewerGeneralStructROIs.load_interaction_data) parses this csv with
            # np.genfromtxt, which naively splits on the delimiter and does not honor csv
            # quoting: the "Types" field must not itself contain a comma, hence the '+' join.
            writer = csv.writer(f)
            writer.writerow(['Chain1', 'Residue1', 'Chain2', 'Residue2', 'Min distance', 'Types'])
            for key, dist in sorted(pairDist.items(), key=lambda kv: kv[1]):
                chain1, resRep1, chain2, resRep2 = key
                writer.writerow([chain1, resRep1, chain2, resRep2, '{:.2f}'.format(dist),
                                 '+'.join(sorted(pairTypes[key]))])
        return outCsv

    def getResidueRepr(self, resNameCsv, resNum):
        '''Build a "RESNAME:RESNUM" residue representation, as used by the PPI branch of
        "Define structural ROIs", from the (1-letter) residue code reported by COCADA'''
        resName3 = RESIDUES1TO3.get(resNameCsv.strip().upper(), resNameCsv.strip())
        return '{}:{}'.format(resName3, resNum.strip())

    def buildBondsPml(self, rowsByType):
        '''Build a PyMol script drawing one dashed distance line per COCADA contact, grouped and
        colored by interaction type, so each type can be shown/hidden independently in PyMol'''
        colors = dict(zip(rowsByType.keys(), createColorVectors(len(rowsByType))))
        pseudoatoms, pseudoCache = [], {}
        colorDefs, distances, groups = [], [], []

        for bondType, rows in rowsByType.items():
            colorName = 'cocada_col_{}'.format(bondType)
            colorDefs.append('set_color {} = {}'.format(colorName, colors[bondType]))

            distNames = []
            for i, row in enumerate(rows, start=1):
                sel1 = self.getPymolSelector(row['Chain1'].strip(), row['Res1'].strip(),
                                             row['Atom1'].strip(), pseudoatoms, pseudoCache)
                sel2 = self.getPymolSelector(row['Chain2'].strip(), row['Res2'].strip(),
                                             row['Atom2'].strip(), pseudoatoms, pseudoCache)
                if sel1 is None or sel2 is None:
                    continue

                distName = 'd_{}_{}'.format(bondType, i)
                distances.append('distance {}, {}, {}'.format(distName, sel1, sel2))
                distances.append('set dash_color, {}, {}'.format(colorName, distName))
                distNames.append(distName)

            if distNames:
                groups.append('group COCADA_{}, {}'.format(bondType, ' '.join(distNames)))

        sticksSel = self.getInterfaceStickSelection(rowsByType)

        lines = [
            'load {}'.format(self.getCifFile()),
            'hide everything',
            'show cartoon',
            'color grey80',
        ]
        if sticksSel:
            lines.append('show sticks, {}'.format(sticksSel))
        lines += pseudoatoms + colorDefs + distances + ['hide labels'] + groups + ['zoom visible']

        pmlFile = self._getExtraPath('cocada_bonds.pml')
        with open(pmlFile, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        return pmlFile

    def getPymolSelector(self, chainId, resNum, atomName, pseudoatomLines, pseudoCache):
        '''Return a PyMol selection string identifying a COCADA atom reference. For the "RNG"
        aromatic-ring pseudo-atom, create (once per residue) a PyMol pseudoatom at the same
        centroid coordinate used to build the StructROIs, and return its object name instead.'''
        if atomName != 'RNG':
            return 'chain {} and resi {} and name {}'.format(chainId, resNum, atomName)

        cacheKey = (chainId, resNum)
        if cacheKey not in pseudoCache:
            coord, _ = self.getAtomCoord(chainId, resNum, atomName)
            if coord is None:
                pseudoCache[cacheKey] = None
            else:
                pseudoName = 'ring_{}_{}'.format(chainId, resNum)
                pseudoatomLines.append('pseudoatom {}, pos={}'.format(pseudoName, coord))
                pseudoCache[cacheKey] = pseudoName
        return pseudoCache[cacheKey]

    def getInterfaceStickSelection(self, rowsByType):
        '''Return a PyMol selection covering all the residues involved in any COCADA contact,
        grouped by chain, so they can be shown as sticks alongside the cartoon representation'''
        resNumsByChain = defaultdict(set)
        for rows in rowsByType.values():
            for row in rows:
                resNumsByChain[row['Chain1'].strip()].add(row['Res1'].strip())
                resNumsByChain[row['Chain2'].strip()].add(row['Res2'].strip())

        chainSels = []
        for chain, resNums in resNumsByChain.items():
            resNums = natural_sort(list(resNums))
            chainSels.append('(chain {} and resi {})'.format(chain, '+'.join(resNums)))
        return ' or '.join(chainSels)
