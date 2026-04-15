# **************************************************************************
# *
# * Authors: Blanca Pueche (blanca.pueche@cnb.csic.es)
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
import json

from pyworkflow.viewer import ProtocolViewer, DESKTOP_TKINTER
from pyworkflow.protocol import params

from pwem.objects import SetOfAtomStructs
from pwem.viewers import EmPlotter
import pyworkflow.viewer as pwviewer

from pwchem.viewers import BaseInteractionViewer
from pwem.protocols import ProtSubSet


class OmniBindPredictionViewer(BaseInteractionViewer):
    _label = 'Interactions viewer'
    _targets = [SetOfAtomStructs]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def _getData(self):
        structSet = self.getStructSet()

        with open(structSet._interactScoresFile.get(), 'r') as f:
            return json.load(f)

    def _getEntityNames(self, data):
        protNames = sorted(data.keys())
        molNames = sorted(next(iter(data.values())).keys())
        scoreTypes = sorted(next(iter(next(iter(data.values())).values())).keys())
        return protNames, molNames, scoreTypes

    def _getLabels(self):
        return "Protein", "Ligand", "Affinity"

    def getStructSet(self):
        if hasattr(self.protocol, 'outputAtomStructs'):
            return self.protocol.outputAtomStructs
        return self.protocol