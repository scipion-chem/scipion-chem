# **************************************************************************
# *
# * Authors:  Carlos Oscar Sorzano (coss@cnb.csic.es)
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

from pathlib import Path
import webbrowser

from pwchem.protocols import ProtChemDali
from pyworkflow.viewer import DESKTOP_TKINTER, ProtocolViewer
from pyworkflow.protocol.params import StringParam

DALISERVER="ekhidna2.biocenter.helsinki.fi"

class ProtChemDaliViewer(ProtocolViewer):
    """ Visualize the output of protocol Dali """
    _label = 'viewer dali'
    _targets = [ProtChemDali]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        ProtocolViewer.__init__(self, **args)

    def _defineParams(self, form):
        form.addSection(label='Visualization')
        # Select the level to show
        form.addParam('url', StringParam,
                      label="URL of Dali results")

    def _getVisualizeDict(self):
        return {'url': self._viewResults}

    def getResultsDir(self):
        fnDir = self.protocol._getExtraPath(DALISERVER)
        fnBaseDir=None
        for fn in Path(fnDir).rglob('*.html'):
            if fn.name=="index.html":
                fnBaseDir=str(fn.parent)
                break
        return fnBaseDir

    def _viewResults(self, e=None):
        import os
        views = []
        fnBaseDir = self.getResultsDir()
        if not fnBaseDir:
            if os.path.exists(self.protocol._getExtraPath(DALISERVER)):
                os.system("rm -rf %s"%self.protocol._getExtraPath(DALISERVER))
            url=self.url.get()
            if not url.endswith("index.html"):
                url+="/index.html"
            os.system('cd %s; wget -r --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 5 %s'%(self.protocol._getExtraPath(),url))
            fnBaseDir = self.getResultsDir()
        if fnBaseDir:
            webbrowser.open_new_tab(os.path.join(fnBaseDir,"index.html"))

            if not hasattr(self.protocol,"outputIds"):
                for fn in Path(fnBaseDir).rglob('*.txt'):
                    self.protocol.constructOutput(str(fn), self.protocol)
        return views

