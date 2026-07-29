# **************************************************************************
# *
# * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# *          Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
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
# * e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

# General imports
import os

# Scipion em imports
from pyworkflow.protocol import params

# Plugin imports
from pwchem import Plugin
from pwchem.constants import RDKIT_DIC
from pwchem.protocols.VirtualDrugScreening.protocol_ligand_shape_distance import ProtocolShapeDistances

scriptName = 'fingerprint_similarity.py'

def parseResults(outputFile):
    scores = {}
    with open(outputFile) as f:
        heads = f.readline().strip().split('\t')
        for h in heads[1:]:
            scores[h] = []

        for line in f:
            sline = line.split("\t")
            for i, h in enumerate(scores):
                simil = float(sline[i+1])
                scores[h].append(1-simil)
    return scores

class ProtocolFingerprintDistance(ProtocolShapeDistances):
    """
    AI Generated:

    This protocol computes the 2D fingerprint-based distance between a reference molecule
    and a set of small molecules.

    It is used to quantify structural similarity in terms of molecular fingerprints,
    providing a simple numerical measure of how chemically similar compounds are.

    Core Concept
    ------------
    Molecular fingerprints:
        Binary or hashed representations of molecular structure capturing the presence
        of chemical substructures (e.g. rings, functional groups).

    Similarity metrics:
        Fingerprints are compared using similarity coefficients (Tanimoto or Dice),
        which produce values between 0 (completely different) and 1 (identical).

    Distance transformation:
        Similarity is converted into distance as:
            distance = 1 - similarity

    Fingerprint Types
    ------------------
    Morgan:
        Circular fingerprints capturing local atomic environments.

    MACCS:
        Key-based structural fingerprints based on predefined substructures.

    Both:
        Runs both fingerprint types for comparison.

    Similarity Coefficients
    -----------------------
    Tanimoto:
        Measures overlap between fingerprint sets; widely used in chemoinformatics.

    Dice:
        Similar to Tanimoto but gives slightly different weighting to shared features.

    Workflow
    --------
    1. Input a set of small molecules.
    2. Select a reference molecule.
    3. Compute molecular fingerprints for all compounds.
    4. Compare each molecule to the reference using the selected coefficient.
    5. Convert similarity values into distances.
    6. Store results as tab-separated output files.

    Output
    ------
    - TSV file containing fingerprint distances for each molecule.
    - Summary file aggregating all computed distances.

    Use Cases
    ---------
    - Ligand-based virtual screening
    - Chemical library diversity analysis
    - Lead compound similarity assessment
    - Clustering of small molecule datasets
    """
    _label = 'Fingerprint distance'
    _fpChoices = ['Morgan', 'MACCS', 'Both']
    _coefChoices = ['Tanimoto', 'Dice', 'Both']
    stepsExecutionMode = params.STEPS_PARALLEL

    ##### -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label='Input')
        form = self.addInputParams(form)

        group = form.addGroup('Reference molecule')
        group.addParam('inputRefSmallMolecules', params.PointerParam,
                       pointerClass='SetOfSmallMolecules', allowsNull=False,
                       label="Input Reference Small Molecules: ",
                       help='Select the molecules where the reference molecule is stored')

        group.addParam('inputReferenceMolecule', params.StringParam,
                       label="Reference molecule: ",
                       help='Model molecule to compare others')

        group = form.addGroup('Descriptors')
        group.addParam('fpChoice', params.EnumParam, default=0, label='Fingerprint type: ',
                       choices=self._fpChoices,
                       help="Chosen fingerprint type to perform the filtering")

        group.addParam('coefChoice', params.EnumParam, default=0, label='Similarity coefficient: ',
                        choices=self._coefChoices,
                        help="Chosen fingerprint type to perform the filtering")

        form.addParallelSection(threads=4)

    # --------------------------- STEPS functions ------------------------------

    def mainStep(self, it):
        paramsPath = self.writeParamsFile(it)
        Plugin.runScript(self, scriptName, paramsPath, env=RDKIT_DIC, cwd=self._getPath())

    # --------------- INFO functions -------------------------
    def _citations(self):
        return []

    def _methods(self):
        return []

    def _summary(self):
        """ Summarize what the protocol has done"""
        summary = []
        if os.path.exists(self._getPath("all_distances.csv")):
            summary.append("The distance results for each metric and molecule can be found at {}".
                           format(self._getPath("all_distances.csv")))
        else:
            summary.append("The protocol has not finished.")
        return summary

    # --------------------------- UTILS functions -----------------------------------
    def writeParamsFile(self, it):
        paramsPath = os.path.abspath(self._getExtraPath(f'inputParams_{it}.txt'))
        molFiles = self.getInputMolFiles(it)

        with open(paramsPath, 'w') as f:
            f.write(f'outputPath: results_{it}.tsv\n')
            f.write('referenceFile: {}\n'.format(self.getRefFile()))
            f.write('ligandFiles: {}\n'.format(' '.join(molFiles)))

            f.write('fingerprint: {}\n'.format(self.getEnumText('fpChoice')))
            f.write('coefficient: {}\n'.format(self.getEnumText('coefChoice')))

        return paramsPath

    def getParseFunc(self):
        parseFunc, args = parseResults, []
        return parseFunc, args

    def getCoefficient(self):
        return self.getEnumText('coefChoice')

