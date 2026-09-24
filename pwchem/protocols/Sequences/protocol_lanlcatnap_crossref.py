# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Enzo Sierra (enzogael57@gmail.com)
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
This protocol is used to cross-reference peptide candidates against known
HIV broadly-neutralizing-antibody (bnAb) linear epitopes, using the local
LANL Immunology DB (+ optional CATNAP neutralization data).

Folded into pwchem core (previously the standalone scipion-chem-lanlcatnap
plugin): wraps no external tool/binary, only two local reference databases
the user downloads manually once (LANL/CATNAP's own terms of use do not
clearly permit redistribution, and hiv.lanl.gov has no stable download URL
/ explicitly discourages automated traffic -- verified live, so this stays
the one exception to auto-installing everything else). A lightweight
protocol with no dedicated conda env belongs in pwchem's own Sequences
protocols, not in a standalone plugin.
"""

import csv
import os
import re
from typing import List, Optional

import pandas as pd
from pwchem.objects import SetOfSequenceROIs
from pwem.protocols import EMProtocol
from pyworkflow.object import Boolean, Integer, String
from pyworkflow.protocol import params
from pyworkflow.utils import Message

from pwchem import Plugin as pwchemPlugin
from pwchem.constants import CATNAP_ABS_PATH, LANL_AB_ALL_PATH

_OUTPUT_COLUMNS = [
    'sequence', 'antibody_name', 'epitope_sequence', 'match_length', 'epitope_name',
    'hxb2_location', 'neutralizing', 'antibody_type', 'subtype', 'binding_region',
    'catnap_mean_ic50', 'catnap_n_viruses',
]

_AA_ONLY = re.compile(r'[A-Zx]+')

_REQUIRED_LANL_COLUMNS = [
    'Antibody name (alias)', 'Epitope', 'Epitope name', 'HXB2 protein location',
    'Neutralizing', 'Antibody type', 'Binding region', 'Subtype',
]

_REQUIRED_CATNAP_COLUMNS = ['Name', 'Mean panel IC50', '# of viruses tested']

DEFAULT_MIN_OVERLAP = 6


class LANLCATNAPParseError(Exception):
    """The LANL/CATNAP reference file does not match the expected format."""


def loadBnabEpitopes(lanlAbAllPath: str) -> pd.DataFrame:
    """Parse 'ab_all.csv', keeping only records with a linear epitope (single-chain AA sequence)."""
    with open(lanlAbAllPath, encoding='utf-8', errors='replace', newline='') as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = list(reader)

    idx = {name: i for i, name in enumerate(header)}
    missing = [c for c in _REQUIRED_LANL_COLUMNS if c not in idx]
    if missing:
        raise LANLCATNAPParseError(f"'{lanlAbAllPath}' does not have the expected columns: missing {missing}.")

    records = []
    for row in rows:
        epitope = row[idx['Epitope']].strip()
        if not epitope or not _AA_ONLY.fullmatch(epitope):
            continue  # conformational epitope, composite notation ('A + B'), or empty -- out of scope
        records.append({
            'antibody_name': row[idx['Antibody name (alias)']].strip(),
            'epitope_sequence': epitope.upper(),
            'epitope_name': row[idx['Epitope name']].strip(),
            'hxb2_location': row[idx['HXB2 protein location']].strip(),
            'neutralizing': row[idx['Neutralizing']].strip(),
            'antibody_type': row[idx['Antibody type']].strip(),
            'binding_region': row[idx['Binding region']].strip(),
            'subtype': row[idx['Subtype']].strip(),
        })
    return pd.DataFrame.from_records(records)


def loadCatnapPotency(catnapAbsPath: str) -> pd.DataFrame:
    """Parse 'abs_YYYY-MM-DD.txt' (CATNAP) to append neutralization potency/breadth per antibody."""
    raw = pd.read_csv(catnapAbsPath, sep='\t', dtype=str)
    missing = [c for c in _REQUIRED_CATNAP_COLUMNS if c not in raw.columns]
    if missing:
        raise LANLCATNAPParseError(f"'{catnapAbsPath}' does not have the expected columns: missing {missing}.")

    result = raw[_REQUIRED_CATNAP_COLUMNS].copy()
    result.columns = ['antibody_name_norm', 'catnap_mean_ic50', 'catnap_n_viruses']
    result['antibody_name_norm'] = result['antibody_name_norm'].str.strip().str.upper()
    result['catnap_mean_ic50'] = pd.to_numeric(result['catnap_mean_ic50'], errors='coerce')
    result['catnap_n_viruses'] = pd.to_numeric(result['catnap_n_viruses'], errors='coerce')
    return result.drop_duplicates(subset='antibody_name_norm', keep='first')


def longestCommonSubstringLen(a: str, b: str) -> int:
    """Length of the longest common substring between 'a' and 'b' (DP O(len(a)*len(b)))."""
    if not a or not b:
        return 0
    prev = [0] * (len(b) + 1)
    best = 0
    for i in range(1, len(a) + 1):
        curr = [0] * (len(b) + 1)
        for j in range(1, len(b) + 1):
            if a[i - 1] == b[j - 1]:
                curr[j] = prev[j - 1] + 1
                best = max(best, curr[j])
        prev = curr
    return best


def buildBnabRecord(seq: str, ref, matchLen: int, potencyDf: Optional[pd.DataFrame]) -> dict:
    """Builds one _OUTPUT_COLUMNS row for a (candidate, reference epitope) match, enriched
    with CATNAP potency data when a matching antibody entry exists in potencyDf."""
    record = {
        'sequence': seq,
        'antibody_name': ref.antibody_name,
        'epitope_sequence': ref.epitope_sequence,
        'match_length': matchLen,
        'epitope_name': ref.epitope_name,
        'hxb2_location': ref.hxb2_location,
        'neutralizing': ref.neutralizing,
        'antibody_type': ref.antibody_type,
        'subtype': ref.subtype,
        'binding_region': ref.binding_region,
        'catnap_mean_ic50': pd.NA,
        'catnap_n_viruses': pd.NA,
    }
    if potencyDf is not None:
        hit = potencyDf[potencyDf['antibody_name_norm'] == ref.antibody_name.strip().upper()]
        if not hit.empty:
            record['catnap_mean_ic50'] = hit.iloc[0]['catnap_mean_ic50']
            record['catnap_n_viruses'] = hit.iloc[0]['catnap_n_viruses']
    return record


def queryBnabCrossref(
    sequences: List[str], lanlAbAllPath: str, catnapAbsPath: Optional[str] = None,
    minOverlap: int = DEFAULT_MIN_OVERLAP,
) -> pd.DataFrame:
    """Cross-reference 'sequences' against known bnAb linear epitopes (LANL Immunology DB).

    Does NOT do structural alignment or HXB2 coordinate mapping: compares
    candidate linear sequences directly against LANL's reported epitopes by
    longest-common-substring overlap.

    Args:
        sequences: Candidate peptides/sequences to evaluate.
        lanlAbAllPath: Path to LANL's 'ab_all.csv'.
        catnapAbsPath: Path to CATNAP's 'abs_YYYY-MM-DD.txt' (optional).
        minOverlap: Minimum substring overlap length to report a match.
            For reference epitopes SHORTER than this, the full epitope
            match is required instead (never a looser threshold than the
            epitope itself).

    Returns:
        DataFrame with one row per (candidate, reference epitope) pair
        that overlaps enough, columns _OUTPUT_COLUMNS. Empty if no match
        or 'sequences' is empty.
    """
    if not sequences:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    bnabDf = loadBnabEpitopes(lanlAbAllPath)
    if bnabDf.empty:
        return pd.DataFrame(columns=_OUTPUT_COLUMNS)

    potencyDf = loadCatnapPotency(catnapAbsPath) if catnapAbsPath else None

    rows = []
    for seq in sequences:
        seqUpper = seq.upper()
        for ref in bnabDf.itertuples(index=False):
            requiredOverlap = min(minOverlap, len(ref.epitope_sequence))
            matchLen = longestCommonSubstringLen(seqUpper, ref.epitope_sequence)
            if matchLen < requiredOverlap:
                continue
            rows.append(buildBnabRecord(seq, ref, matchLen, potencyDf))

    return pd.DataFrame(rows, columns=_OUTPUT_COLUMNS) if rows else pd.DataFrame(columns=_OUTPUT_COLUMNS)


class ProtLANLCATNAPCrossref(EMProtocol):
    """
    AI Generated:

    Cross-references every input ROI's peptide against known HIV
    broadly-neutralizing-antibody (bnAb) linear epitopes (LANL Immunology
    DB, optionally enriched with CATNAP neutralization potency, Yoon et
    al. 2015), and annotates (does NOT filter) each ROI with a summary of
    its best match. Purely informative for any pipeline: only relevant
    when the input is HIV Env, but harmless (zero matches) otherwise.

    Output
    ------
    outputROIs: the same SetOfSequenceROIs as the input, annotated with
    '_bnabMatchCount' (int, total matching reference epitopes),
    '_bnabNeutralizingMatch' (bool, True if >=1 match is a confirmed
    neutralizing antibody), '_bnabBestAntibody'/'_bnabBestMatchLength'
    (the single longest match, name and length). The full detail (every
    match, not just the best one) is persisted to
    'extra/bnab_crossref.csv'.
    """

    _label = 'lanl-catnap bnab crossref'

    def _defineParams(self, form):
        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputROIs', params.PointerParam, pointerClass='SetOfSequenceROIs',
                       label='Sequence ROIs: ',
                       help='Peptide candidates to cross-reference against known bnAb epitopes.')
        form.addParam('minOverlap', params.IntParam, default=DEFAULT_MIN_OVERLAP,
                       label='Min. substring overlap (aa): ',
                       help='Minimum substring overlap to report a match. Reference epitopes '
                            'shorter than this still require their own full length as the match.')

    def _insertAllSteps(self):
        self._insertFunctionStep(self.crossrefStep)
        self._insertFunctionStep(self.createOutputStep)

    # ---------------------------------- Steps -----------------------------------

    def _getCrossrefPath(self):
        return self._getExtraPath('bnab_crossref.csv')

    def _getRois(self):
        # Iterating a Scipion SetOfXXX reuses the same Python object per row
        # (the underlying sqlite cursor): each item must be cloned when
        # materialized into a list, or all N references end up pointing to
        # the cursor's last state.
        return [roi.clone() for roi in self.inputROIs.get()]

    def crossrefStep(self):
        rois = self._getRois()
        sequences = [roi.getROISequence() for roi in rois]
        if not sequences:
            return

        catnapPath = pwchemPlugin.getVar(CATNAP_ABS_PATH) or None
        crossrefDf = queryBnabCrossref(
            sequences, lanlAbAllPath=pwchemPlugin.getVar(LANL_AB_ALL_PATH),
            catnapAbsPath=catnapPath, minOverlap=self.minOverlap.get(),
        )
        crossrefDf.to_csv(self._getCrossrefPath(), index=False)

    def createOutputStep(self):
        rois = self._getRois()
        crossrefDf = pd.read_csv(self._getCrossrefPath()) if os.path.isfile(self._getCrossrefPath()) else pd.DataFrame()

        outROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))
        for roi in rois:
            matches = crossrefDf[crossrefDf['sequence'] == roi.getROISequence()] if not crossrefDf.empty else crossrefDf

            roi._bnabMatchCount = Integer(len(matches))
            roi._bnabNeutralizingMatch = Boolean(bool((matches['neutralizing'] == 'yes').any()) if len(matches) else False)
            if len(matches):
                best = matches.sort_values('match_length', ascending=False).iloc[0]
                roi._bnabBestAntibody = String(best['antibody_name'])
                roi._bnabBestMatchLength = Integer(int(best['match_length']))
            else:
                roi._bnabBestAntibody = String('')
                roi._bnabBestMatchLength = Integer(0)
            outROIs.append(roi)

        if len(outROIs) > 0:
            self._defineOutputs(outputROIs=outROIs)
            self._defineSourceRelation(self.inputROIs, outROIs)

    # ---------------------------------- Validation -------------------------------

    def _validate(self):
        errors = []
        lanlPath = pwchemPlugin.getVar(LANL_AB_ALL_PATH)
        if not lanlPath or not os.path.isfile(lanlPath):
            errors.append(
                f"LANL_AB_ALL_PATH is not set or does not exist: '{lanlPath}'. Download 'ab_all.csv' from "
                "https://www.hiv.lanl.gov/content/immunology/ and set LANL_AB_ALL_PATH in scipion.conf."
            )
        catnapPath = pwchemPlugin.getVar(CATNAP_ABS_PATH)
        if catnapPath and not os.path.isfile(catnapPath):
            errors.append(f"CATNAP_ABS_PATH is set but does not exist: '{catnapPath}'.")
        return errors

    def _summary(self):
        summary = []
        if self.isFinished():
            outROIs = getattr(self, 'outputROIs', None)
            if outROIs is not None:
                nMatch = sum(1 for roi in outROIs if roi._bnabMatchCount.get() > 0)
                nNeutralizing = sum(1 for roi in outROIs if roi._bnabNeutralizingMatch.get())
                summary.append(f'{nMatch}/{len(outROIs)} candidate(s) match >= 1 known bnAb epitope '
                               f'({nNeutralizing} with a confirmed neutralizing antibody).')
        return summary
