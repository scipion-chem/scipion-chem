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
import json, os
import tkinter as tk
from tkinter import filedialog

from pyworkflow import Config
from pyworkflow.object import String
from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.gui.gui import configureWeigths, getDefaultFont
from pyworkflow.utils import Icon
from pyworkflow.gui.tree import TreeProvider
from pwem.wizards import VariableWizard

from pwchem.protocols.Sequences.protocol_download_genomes import ProtDownloadGenomes
from pwchem.protocols.Sequences.protocol_download_vcf import ProtDownloadVCF
from pwchem.protocols.Sequences.protocol_rnaseq_align import ProtRNASeqAlignment
from pwchem.protocols.Sequences.protocol_import_genomes import ProtImportGenomes
from pwchem.protocols.Sequences.protocol_gatk import ProtGATK
from pwchem.protocols.Sequences.protocol_import_vcf import ProtImportVCF
from pwchem.utils.utilsRNA import getCommonSpecies

ALL_FILES = ('All files', '*')

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
    protocol=ProtDownloadGenomes,
    targets=['commonSpecies'],
    inputs=[],
    outputs=['commonSpecies']
)
SelectCommonSpeciesWizard().addTarget(
    protocol=ProtDownloadVCF,
    targets=['commonSpecies'],
    inputs=[],
    outputs=['commonSpecies']
)

class SelectItemFromSetWizard(VariableWizard):
    """
    Generic wizard to select one item from a Scipion Set.

    The selected item is stored as its zero-based index in the
    corresponding protocol parameter.
    """

    _targets, _inputs, _outputs = [], {}, {}

    def show(self, form, *params):
        protocol = form.protocol
        inputParams, outputParams = self.getInputOutput(form)

        inputSet = getattr(protocol, inputParams[0]).get()

        if inputSet is None:
            dialog.showError(
                'Input set',
                'Please select an input set first.',
                form.root
            )
            return

        if inputSet.getSize() == 0:
            dialog.showError(
                'Input set',
                'The selected set is empty.',
                form.root
            )
            return

        labels = []

        for index, item in enumerate(inputSet):
            label = self._getItemLabel(index, item)
            labels.append(String(label))

        provider = ListTreeProviderString(labels)

        dlg = dialog.ListDialog(
            form.root,
            'Available items',
            provider,
            'Select one item'
        )

        if not dlg.resultYes() or not dlg.values:
            return

        selectedLabel = dlg.values[0].get()

        selectedIndex = int(
            selectedLabel.split(' - ', 1)[0]
        )

        form.setVar(
            outputParams[0],
            selectedIndex
        )

    @staticmethod
    def _getItemLabel(index, item):
        """
        Build a descriptive label using the metadata available
        in the object.

        This works for Genome, VCFFile and other compatible
        Scipion objects.
        """
        parts = [str(index)]

        # Genome / VCF metadata
        for getterName in (
                'getScientificName',
                'getAssembly',
                'getRelease',
                'getSource'
        ):
            if hasattr(item, getterName):
                value = getattr(item, getterName)()

                if value:
                    parts.append(str(value))

        # Generic fallback
        if len(parts) == 1 and hasattr(item, 'getFileName'):
            value = item.getFileName()

            if value:
                parts.append(str(value))

        return ' - '.join(parts)


# -------------------------------------------------------------------------
# RNA-seq alignment
# SetOfGenomes -> selected Genome
# -------------------------------------------------------------------------

SelectItemFromSetWizard().addTarget(
    protocol=ProtRNASeqAlignment,
    targets=['genomeIndex'],
    inputs=['inputGenomes'],
    outputs=['genomeIndex']
)


# -------------------------------------------------------------------------
# GATK
# SetOfVCFFiles -> selected VCFFile
# -------------------------------------------------------------------------

SelectItemFromSetWizard().addTarget(
    protocol=ProtGATK,
    targets=['vcfIndex'],
    inputs=['knownSites'],
    outputs=['vcfIndex']
)

class ImportEntry:
    """Temporary entry used by the sequence-file import wizard."""

    def __init__(self, **kwargs):
        self.data = kwargs

    def get(self, key, default=''):
        return self.data.get(key, default)

    def set(self, key, value):
        self.data[key] = value

    def toDict(self):
        return dict(self.data)

    @classmethod
    def fromDict(cls, data):
        return cls(**data)


class ImportTreeProvider(TreeProvider):
    """Provide entries to the sequence-file import dialog."""

    def __init__(self, entries, fields):
        super().__init__()
        self.entries = entries
        self.fields = fields

    def getColumns(self):
        return [
            (
                field['label'],
                field.get('width', 150)
            )
            for field in self.fields
        ]

    def getObjects(self):
        return self.entries

    def getObjectInfo(self, entry):
        values = [
            entry.get(field['name'])
            for field in self.fields
        ]

        return {
            'key': str(id(entry)),
            'text': values[0] if values else '',
            'values': tuple(values[1:])
        }


class ImportFilesDialog(dialog.ToolbarListDialog):
    """Dialog to add, edit and remove import entries."""

    def __init__(
            self,
            parent,
            entries,
            fields,
            title,
            description,
            itemName,
            **kwargs
    ):
        self.entries = entries
        self.fields = fields
        self.itemName = itemName

        toolbarButtons = [
            dialog.ToolbarButton(
                'Add {}'.format(itemName),
                self._addItem,
                Icon.ACTION_NEW
            ),
            dialog.ToolbarButton(
                'Edit',
                self._editItem,
                Icon.ACTION_EDIT
            ),
            dialog.ToolbarButton(
                'Delete',
                self._deleteItem,
                Icon.ACTION_DELETE
            )
        ]

        super().__init__(
            parent,
            title,
            ImportTreeProvider(
                entries,
                fields
            ),
            description,
            toolbarButtons,
            allowsEmptySelection=True,
            itemDoubleClick=self._editItem,
            buttons=[
                ('Save', dialog.RESULT_YES),
                ('Cancel', dialog.RESULT_CANCEL)
            ],
            **kwargs
        )

    def _addItem(self, e=None):
        entry = ImportEntry()

        dlg = EditImportEntryDialog(
            self,
            'Add {}'.format(self.itemName),
            entry,
            self.fields
        )

        if dlg.resultYes():
            self.entries.append(entry)
            self.tree.update()

    def _editItem(self, e=None):
        selection = self.tree.getSelectedObjects()

        if selection:
            entry = selection[0]

            dlg = EditImportEntryDialog(
                self,
                'Edit {}'.format(self.itemName),
                entry,
                self.fields
            )

            if dlg.resultYes():
                self.tree.update()

    def _deleteItem(self, e=None):
        selection = self.tree.getSelectedObjects()

        if selection:
            for entry in selection:
                self.entries.remove(entry)

            self.tree.update()


class EditImportEntryDialog(dialog.Dialog):
    """Dialog to configure one import entry."""

    def __init__(
            self,
            parent,
            title,
            entry,
            fields,
            **kwargs
    ):
        self.entry = entry
        self.fields = fields
        self.variables = {}

        super().__init__(
            parent,
            title,
            **kwargs
        )

    def body(self, bodyFrame):
        bodyFrame.config(
            bg=Config.SCIPION_BG_COLOR
        )

        configureWeigths(
            bodyFrame,
            1,
            1
        )

        for row, field in enumerate(self.fields):
            value = self.entry.get(
                field['name']
            )

            if field.get('type') == 'file':
                variable = self._addFileEntry(
                    bodyFrame,
                    row,
                    field['label'],
                    value,
                    field.get('fileTypes')
                )
            else:
                variable = self._addEntry(
                    bodyFrame,
                    row,
                    field['label'],
                    value
                )

            self.variables[
                field['name']
            ] = variable

    @staticmethod
    def _addEntry(
            parent,
            row,
            label,
            value
    ):
        tk.Label(
            parent,
            text=label,
            bg=Config.SCIPION_BG_COLOR
        ).grid(
            row=row,
            column=0,
            sticky='w',
            padx=(15, 10),
            pady=8
        )

        var = tk.StringVar(
            value=value
        )

        tk.Entry(
            parent,
            width=60,
            font=getDefaultFont(),
            textvariable=var
        ).grid(
            row=row,
            column=1,
            columnspan=2,
            sticky='we',
            padx=(5, 15),
            pady=8
        )

        return var

    @staticmethod
    def _addFileEntry(
            parent,
            row,
            label,
            value,
            fileTypes=None
    ):
        tk.Label(
            parent,
            text=label,
            bg=Config.SCIPION_BG_COLOR
        ).grid(
            row=row,
            column=0,
            sticky='w',
            padx=(15, 10),
            pady=8
        )

        var = tk.StringVar(
            value=value
        )

        tk.Entry(
            parent,
            width=50,
            font=getDefaultFont(),
            textvariable=var
        ).grid(
            row=row,
            column=1,
            sticky='we',
            padx=(5, 5),
            pady=8
        )

        def browse():
            currentPath = var.get().strip()

            if (
                currentPath
                and os.path.isfile(currentPath)
            ):
                initialDir = os.path.dirname(
                    currentPath
                )
            else:
                initialDir = os.path.expanduser(
                    '~'
                )

            selectedPath = filedialog.askopenfilename(
                parent=parent,
                title='Select {}'.format(label),
                initialdir=initialDir,
                filetypes=fileTypes or [
                    ALL_FILES
                ]
            )

            if selectedPath:
                var.set(
                    selectedPath
                )

        tk.Button(
            parent,
            text='Browse',
            command=browse
        ).grid(
            row=row,
            column=2,
            sticky='w',
            padx=(5, 15),
            pady=8
        )

        return var

    def apply(self):
        for field in self.fields:
            name = field['name']

            self.entry.set(
                name,
                self.variables[name].get().strip()
            )


class ImportSequenceFilesWizard(VariableWizard):
    """Wizard to configure sequence-related files to import."""

    _targets, _inputs, _outputs = [], {}, {}

    GENOME_FIELDS = [
        {
            'name': 'scientificName',
            'label': 'Scientific name',
            'width': 180
        },
        {
            'name': 'assembly',
            'label': 'Assembly',
            'width': 100
        },
        {
            'name': 'release',
            'label': 'Release',
            'width': 80
        },
        {
            'name': 'fastaFile',
            'label': 'FASTA',
            'width': 200,
            'type': 'file',
            'fileTypes': [
                (
                    'FASTA files',
                    '*.fa *.fasta *.fna'
                ),
                ALL_FILES

            ]
        },
        {
            'name': 'gtfFile',
            'label': 'GTF',
            'width': 200,
            'type': 'file',
            'fileTypes': [
                (
                    'GTF files',
                    '*.gtf'
                ),
                ALL_FILES
            ]
        }
    ]

    VCF_FIELDS = [
        {
            'name': 'scientificName',
            'label': 'Scientific name',
            'width': 180
        },
        {
            'name': 'assembly',
            'label': 'Assembly',
            'width': 100
        },
        {
            'name': 'vcfFile',
            'label': 'VCF',
            'width': 200,
            'type': 'file',
            'fileTypes': [
                (
                    'VCF files',
                    '*.vcf *.vcf.gz'
                ),
                ALL_FILES
            ]
        },
        {
            'name': 'indexFile',
            'label': 'Index',
            'width': 200,
            'type': 'file',
            'fileTypes': [
                (
                    'VCF index files',
                    '*.tbi *.csi *.idx'
                ),
                ALL_FILES
            ]
        }
    ]

    def show(self, form, *params):
        protocol = form.protocol

        if isinstance(
            protocol,
            ProtImportGenomes
        ):
            parameterName = 'genomesData'
            fields = self.GENOME_FIELDS
            title = 'Import genomes'
            description = (
                'Add one or more reference genomes.'
            )
            itemName = 'genome'

        elif isinstance(
            protocol,
            ProtImportVCF
        ):
            parameterName = 'vcfsData'
            fields = self.VCF_FIELDS
            title = 'Import VCFs'
            description = (
                'Add one or more known-variant VCF files.'
            )
            itemName = 'VCF'

        else:
            return

        rawData = getattr(
            protocol,
            parameterName
        ).get()

        try:
            data = (
                json.loads(rawData)
                if rawData
                else []
            )
        except (TypeError, ValueError):
            data = []

        entries = [
            ImportEntry.fromDict(item)
            for item in data
        ]

        dlg = ImportFilesDialog(
            form.root,
            entries,
            fields,
            title,
            description,
            itemName
        )

        if dlg.resultYes():
            value = json.dumps([
                entry.toDict()
                for entry in entries
            ])

            form.setVar(
                parameterName,
                value
            )

ImportSequenceFilesWizard().addTarget(
    protocol=ProtImportGenomes,
    targets=['genomesData'],
    inputs=[],
    outputs=['genomesData']
)

ImportSequenceFilesWizard().addTarget(
    protocol=ProtImportVCF,
    targets=['vcfsData'],
    inputs=[],
    outputs=['vcfsData']
)