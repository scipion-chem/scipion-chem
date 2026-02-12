# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          Irene Sánchez Martín (100495638@alumnos.uc3m.es)
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
# **************************************************************************

import numpy as np
import os
import json
from rdkit import Chem
from rdkit.Chem import Descriptors

from pyworkflow.protocol import EMProtocol
from pyworkflow.object import Float
from pwchem.objects import SetOfSmallMolecules
from pwchem.constants import DESCRIPTOR_CATEGORIES


class ProtocolLigandCharacterization(EMProtocol):
    """
    Computes RDKit descriptors for a SetOfSmallMolecules and adds
    them as new properties accessible and visible in Scipion GUI.
    """

    _label = "Ligand characterization"

    # ------------------------------
    # STEP 1 — Compute descriptors
    # ------------------------------
    def calculateDescriptorsStep(self):
        inputSet = self.inputSmallMolecules.get()

        header = []
        property_dict = {}

        # Use all RDKit descriptor names in DESCRIPTOR_CATEGORIES
        descriptor_names = []
        for cat, props in DESCRIPTOR_CATEGORIES.items():
            descriptor_names.extend(props)

        header = descriptor_names

        for mol in inputSet:
            molName = mol.molName.get()
            sdfPath = mol.smallMoleculeFile.get()

            rdkit_mol = Chem.MolFromMolFile(sdfPath, sanitize=True)
            if rdkit_mol is None:
                continue

            row = []
            for d in descriptor_names:
                try:
                    fn = getattr(Descriptors, d)
                    value = float(fn(rdkit_mol))
                except Exception:
                    value = None
                row.append(value)

            property_dict[molName] = row

        # Save JSON
        output_json = self._getExtraPath("output_results.json")
        with open(output_json, "w") as f:
            json.dump({"header": header, "property_dict": property_dict},
                      f, indent=2)

        self.info("Descriptor JSON written.")


    # ----------------------------------------------------------
    # STEP 2 — Create output Set with descriptor properties
    # ----------------------------------------------------------
    def createOutputStep(self):
        inp = self.inputSmallMolecules.get()

        # Copy input set
        newMols = SetOfSmallMolecules.createCopy(inp, self._getPath(),
                                                 copyInfo=True)

        # Load results
        output_json = self._getExtraPath("output_results.json")
        if not os.path.exists(output_json):
            self.warning("Descriptor JSON missing.")
            self._defineOutputs(outputSmallMolecules=newMols)
            return

        with open(output_json, "r") as f:
            data = json.load(f)

        header = data.get("header", [])
        property_dict = data.get("property_dict", {})

        if not header or not property_dict:
            self.warning("Descriptor JSON empty.")
            self._defineOutputs(outputSmallMolecules=newMols)
            return

        # ---------------------------
        # Register columns in Scipion
        # ---------------------------
        for prop_name in header:
            category = next(
                (cat for cat, props in DESCRIPTOR_CATEGORIES.items()
                 if prop_name in props),
                "uncategorized"
            )
            attr_name = f"Property_{category}_{prop_name}"

            if not newMols.hasProperty(attr_name):
                newMols.addProperty(attr_name, Float, default=None)

        matched = 0
        unmatched = 0

        # ---------------------------
        # Assign descriptor values
        # ---------------------------
        for mol in newMols:
            molNameObj = getattr(mol, "molName", None)
            if molNameObj is None:
                unmatched += 1
                continue

            molName = str(molNameObj.get())
            if molName not in property_dict:
                unmatched += 1
                continue

            row = property_dict[molName]

            for i, prop_name in enumerate(header):
                category = next(
                    (cat for cat, props in DESCRIPTOR_CATEGORIES.items()
                     if prop_name in props),
                    "uncategorized"
                )
                attr_name = f"Property_{category}_{prop_name}"

                try:
                    value = float(row[i]) if row[i] is not None else None
                except:
                    value = None

                setattr(mol, attr_name, Float(value))

            matched += 1

        newMols.updateMolClass()
        self.info(f"Descriptor mapping complete: {matched} matched, {unmatched} unmatched.")

        self._defineOutputs(outputSmallMolecules=newMols)


    # ---------------------------
    # Define protocol steps
    # ---------------------------
    def _defineSteps(self):
        self._insertFunctionStep("calculateDescriptorsStep")
        self._insertFunctionStep("createOutputStep")


    # ---------------------------
    # Define protocol inputs
    # ---------------------------
    def _defineParams(self, form):
        form.addSection("Input")
        form.addParam("inputSmallMolecules", form.params.PointerParam,
                      label="Input molecules",
                      pointerClass="SetOfSmallMolecules",
                      help="Set of molecules to characterize.")


    # ---------------------------
    # Summary
    # ---------------------------
    def _summary(self):
        return ["Computes RDKit descriptors and adds them to the output set."]

