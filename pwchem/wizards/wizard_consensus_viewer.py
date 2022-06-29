# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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

"""

# Imports
from pwem.wizards.wizard import EmWizard
from pwem.wizards.wizards_3d.mask_structure_wizard import MaskStructureWizard
from ..viewers import ViewerConsensusStructROIs


class SetOutputClass(EmWizard):
    _targets = [(ViewerConsensusStructROIs, ['setClass'])]

    def show(self, form):
        viewer = form.protocol

        choices = form.widgetDict['outputSet'].param.choices
        outSetLabel = choices[form.widgetDict['outputSet'].get()]
        outSet = getattr(viewer.protocol, outSetLabel)
        outSetClass = outSet.getPocketsClass()

        form.setVar('setClass', outSetClass)



