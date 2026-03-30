# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:     Your Name (your.email@cnb.csic.es)
# *
# * Unidad de Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 3 of the License, or
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

import os
from pyworkflow.viewer import TextView
from pyworkflow.gui.text import TextFileViewer

from pyworkflow.gui.dialog import showError
from pwem.wizards import VariableWizard
from pwchem.protocols.MolecularDynamics.protocol_trajectory_clustering import ProtocolTrajectoryClustering
from pwchem.viewers.viewers_MD import MDSystemPViewer

class OpenSystemFileWizard(VariableWizard):
    # Standard Scipion wizard storage
    _targets = []
    _inputs = {}
    _outputs = {}

    def show(self, form, *params):
        inputParam, outputParam = self.getInputOutput(form)
        system = inputParam[0]
        protocol = form.protocol

        systemFile = getattr(protocol, inputParam[0]).get().getSystemFile()

        if systemFile and os.path.exists(systemFile):
            view = TextView([systemFile], title=f"Viewing: {os.path.basename(systemFile)}", tkParent=form)
            view.show()
        else:
            showError("File Error", f"File path invalid: {systemFile}", form.window)

OpenSystemFileWizard().addTarget(
    protocol=ProtocolTrajectoryClustering,
    targets=['firstResAlign'],
    inputs=['inputMDSystem'],
    outputs=['']
)

OpenSystemFileWizard().addTarget(
    protocol=ProtocolTrajectoryClustering,
    targets=['firstResRmsd'],
    inputs=['inputMDSystem'],
    outputs=['']
)


class OpenViewerSystemFileWizard(VariableWizard):
    _targets = []
    _inputs  = {}
    _outputs = {}

    def show(self, form, *params):
        # form.protocol is the MDSystemPViewer instance when targeting a viewer
        viewer = form.protocol

        try:
            mdSystem = viewer.getMDSystem()
            systemFile = mdSystem.getSystemFile()
        except Exception as e:
            showError("Wizard Error", f"Could not retrieve system file:\n{e}", form.window)
            return
        view = TextView([systemFile], title=f"Viewing: {os.path.basename(systemFile)}", tkParent=form)
        view.show()


OpenViewerSystemFileWizard().addTarget(
    protocol=MDSystemPViewer,
    targets=['atom1'],
    inputs=['atom1'],
    outputs=['atom1'],
)