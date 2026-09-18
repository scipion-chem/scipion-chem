# **************************************************************************
# *
# * Authors:    Laura Pérez Liens (laura.perez@cnb.csic.es)
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

from pyworkflow.object import String
from pyworkflow.gui import ListTreeProviderString, dialog

from pwem.wizards import VariableWizard

from pwchem.protocols.Sequences.protocol_reference_genomes import ProtReferenceGenomes
from pwchem.protocols.Sequences.protocol_vcf import ProtVCF
from pwchem.protocols.Sequences.protocol_rnaseq_align import ProtRNASeqAlignment
from pwchem.utils.sequence_utils import getCommonSpecies


class SelectCommonSpeciesWizard(VariableWizard):
    """Wizard to select one or more common species."""

    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        commonSpecies = getCommonSpecies()

        finalList = [
            String(species)
            for species in commonSpecies
        ]

        provider = ListTreeProviderString(finalList)

        dlg = dialog.ListDialog(
            form.root,
            "Common species",
            provider,
            "Select one or more species",
            selectmode="extended"
        )

        if dlg.resultYes():
            selected = [
                obj.get()
                for obj in dlg.values
            ]

            form.setVar(
                'commonSpecies',
                ';'.join(selected)
            )


SelectCommonSpeciesWizard().addTarget(
    protocol=ProtReferenceGenomes,
    targets=['commonSpecies'],
    inputs=[],
    outputs=['commonSpecies']
)
SelectCommonSpeciesWizard().addTarget(
    protocol=ProtVCF,
    targets=['commonSpecies'],
    inputs=[],
    outputs=['commonSpecies']
)

class SelectGenomeFromSetWizard(VariableWizard):
    """Wizard to select one genome from a SetOfGenomes."""

    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParams = self.getInputOutput(form)

        genomeSet = getattr(protocol, inputParams[0]).get()

        if genomeSet is None:
            dialog.showError(
                "Input genomes",
                "Please select a SetOfGenomes first.",
                form.root
            )
            return

        if genomeSet.getSize() == 0:
            dialog.showError(
                "Input genomes",
                "The selected SetOfGenomes is empty.",
                form.root
            )
            return

        labels = []

        # Do not use list(genomeSet):
        # Scipion may reuse the same object while iterating over an EMSet.
        for index, genome in enumerate(genomeSet):
            species = (
                genome.getScientificName()
                if genome.hasScientificName()
                else "Unknown species"
            )

            assembly = (
                genome.getAssembly()
                if genome.hasAssembly()
                else "Unknown assembly"
            )

            release = (
                genome.getRelease()
                if genome.hasRelease()
                else None
            )

            label = "{} - {} ({})".format(
                index,
                species,
                assembly
            )

            if release:
                label += " - Ensembl {}".format(release)

            labels.append(String(label))

        provider = ListTreeProviderString(labels)

        dlg = dialog.ListDialog(
            form.root,
            "Available genomes",
            provider,
            "Select one reference genome"
        )

        if not dlg.resultYes() or not dlg.values:
            return

        selectedLabel = dlg.values[0].get()

        # The displayed label starts with:
        # 0 - Homo sapiens ...
        # 1 - Mus musculus ...
        selectedIndex = int(
            selectedLabel.split(" - ", 1)[0]
        )

        form.setVar(
            outputParams[0],
            selectedIndex
        )


SelectGenomeFromSetWizard().addTarget(
    protocol=ProtRNASeqAlignment,
    targets=['genomeIndex'],
    inputs=['inputGenomes'],
    outputs=['genomeIndex']
)