# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
This protocol is used to select the regions of the protein that are extracellular, transmembrane and intracellular

"""
import json
import requests
from fontTools.varLib.plot import stops
from pyworkflow.protocol import params
from pyworkflow.utils import Message
from pwem.protocols import EMProtocol

from pwchem.objects import SetOfSequenceROIs, Sequence,SequenceROI


class ProtGetSequenceLocation(EMProtocol):
    """
    AI Generated:

This protocol is used to retrieve and annotate functional regions of protein sequences
based on UniProt topology information.

It automatically identifies extracellular, transmembrane, and intracellular regions
by querying the UniProtKB REST API and mapping annotated features back onto the
input protein sequence as Sequence ROIs (Regions Of Interest).

The resulting regions can be used to study protein topology, membrane organization,
and functional compartmentalization of proteins.

Workflow
--------
1. Input handling:
   - A protein sequence is provided as a Scipion Sequence object.
   - The user specifies whether the sequence originates from UniProtKB.
   - If not, a manual UniProt accession ID is required.

2. Data retrieval:
   - The protocol queries the UniProt REST API:
     https://rest.uniprot.org/uniprotkb/
   - It retrieves JSON-formatted feature annotations including:
     - Transmembrane regions
     - Topological domains (extracellular/intracellular regions)

3. Feature extraction:
   - The JSON response is parsed to extract annotated sequence regions.
   - Each feature provides start/end positions and a biological label.
   - Duplicate or repeated annotations are handled by indexing.

4. Region mapping:
   - Each annotated region is mapped back onto the original sequence.
   - Sequence slices are extracted according to UniProt coordinates.
   - Region type (e.g., transmembrane, extracellular) is preserved.

5. ROI creation:
   - Each region is converted into a SequenceROI object.
   - ROIs include:
     - Start and end positions
     - Region type annotation
     - Derived subsequence
     - Updated description with topology information

6. Output generation:
   - All ROIs are stored in a SetOfSequenceROIs database.
   - The output preserves linkage to the original input sequence.
   - Each ROI represents a biologically meaningful structural domain.

Key Concepts
------------
UniProtKB:
    A comprehensive protein sequence and functional annotation database.

Topology Features:
    Annotated regions describing membrane orientation:
    - Extracellular regions (outside the cell)
    - Transmembrane helices (within lipid bilayer)
    - Intracellular regions (cytoplasmic side)

Sequence ROI:
    A subregion of a protein sequence defined by start and end positions.

REST API:
    Web service used to retrieve protein annotations in JSON format.

Output
------
outputROIs:
    A SetOfSequenceROIs containing all detected topological regions.

Each ROI corresponds to a biologically annotated segment of the protein sequence.

Use Cases
---------
- Protein topology analysis
- Membrane protein characterization
- Functional domain annotation
- Guiding mutagenesis or epitope design
- Structural and functional protein studies
"""
    _label = 'Cellular location of sequence regions'
    BASE_URL = "https://rest.uniprot.org/uniprotkb/"

    # -------------------------- DEFINE param functions ----------------------
    def _defineParams(self, form):
        """ """
        form.addSection(label=Message.LABEL_INPUT)
        group = form.addGroup('Input')
        group.addParam('uniprotSeq', params.BooleanParam, deafult=True,
                       label='Sequence was downloaded from UniProtKB: ',
                       help='Choose if the sequence was downloaded from UniProtKB.')
        group.addParam('inputSequence', params.PointerParam, pointerClass='Sequence',
                       allowsNull=True, label="Input sequence: ",
                       help='Select the sequence to find the cellular locations. In case the sequence WAS NOT downloaded from UniProt, enter manual UniProt ID.')
        group.addParam('uniprotID', params.StringParam, condition= 'not uniprotSeq',
                       allowsNull=True, label="Input UniProt ID: ",
                       help='Select the UniProt ID for the sequence of interest.')

    # --------------------------- STEPS functions ------------------------------
    def _insertAllSteps(self):
        self._insertFunctionStep('defineOutputStep')


    def defineOutputStep(self):
        # og sequence data
        inputSeq = self.inputSequence.get()
        seqStr = inputSeq.getSequence()
        if self.uniprotSeq.get():
            seqId = inputSeq.getId()
        else:
            seqId = self.uniprotID.get()
        seqDesc = inputSeq.getDescription() or ""
        alphabet = inputSeq.getAlphabet()
        isAmino = inputSeq.getIsAminoacids()


        data = self.downloadUniprotJson(seqId)

        # obtain each regions aa
        ranges = self.obtainRanges(data)
        if not ranges:
            self.error("No topology regions found in UniProt record!")
            return

        seqROIs = SetOfSequenceROIs(filename=self._getPath('sequenceROIs.sqlite'))

        for desc, info in ranges.items():
            start = info["start"]
            end = info["end"]
            regionType = info["type"]
            regionSeqStr = seqStr[start-1:end]
            regionSeq = Sequence(name=f"{desc}_{regionType}",
                                 sequence=regionSeqStr,
                                 alphabet=alphabet,
                                 isAminoacids=isAmino,
                                 id=f"{seqId}_{desc}",
                                 description=f"{seqDesc} | {desc} ({regionType}: {start}-{end})"
                                 )
            seqROI = SequenceROI(sequence=inputSeq, seqROI=regionSeq, type=regionType)
            seqROI.setROIIdx(start)
            seqROI.setROIIdx2(end)
            seqROIs.append(seqROI)

        self._defineOutputs(outputROIs=seqROIs)
        self._defineSourceRelation(self.inputSequence, seqROIs)

    # --------------------------- INFO functions -----------------------------------
    def _summary(self):
        summary = []
        return summary

    def _methods(self):
        methods = []
        return methods

    def _warnings(self):
        warns = []
        return warns

    def _validate(self):
        errors = []
        return errors

    # --------------------------- UTILS functions -----------------------------------
    def downloadUniprotJson(self, seqId):
        self.BASE_URL = f"{self.BASE_URL}{seqId}.json?fields=ft_topo_dom%2Cft_transmem"

        results = self._getExtraPath(f"{seqId}.json")

        response = requests.get(self.BASE_URL)
        if response.status_code == 200:
            data = response.json()
            if results:
                with open(results, 'w') as f:
                    json.dump(data, f, indent=2)
                print(f"Saved JSON to {results}")
                return data
        else:
            print(f"Failed to download JSON: {response.status_code}")
            return None

    def obtainRanges(self, data):
        features = data.get("features", [])
        regionRanges = {}
        descCount = {}

        for i, feature in enumerate(features, start=1):
            start = feature["location"]["start"]["value"]
            end = feature["location"]["end"]["value"]
            description = feature.get("description", f"Region_{i}")
            ftype = feature.get("type", "Unknown")

            # handle duplicate descriptions
            count = descCount.get(description, 0) + 1
            descCount[description] = count

            key = description if count == 1 else f"{description}_{count}"

            regionRanges[key] = {
                "start": start,
                "end": end,
                "type": ftype
            }

        return regionRanges