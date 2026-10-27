# **************************************************************************
# *
# * Authors:     Laura Pérez Liens (laura.perez@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
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
from pathlib import Path
import webbrowser

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER

from pwchem.objects import FastqFile


class FastqHtmlViewer(Viewer):
    """Open FASTQ HTML reports stored in a FastqFile object.

    All reports are checked before the browser is launched. Duplicate paths
    open once, and local file URIs preserve spaces, accents and URL characters.
    """

    _label = 'FASTQ HTML reports viewer'
    _targets = [FastqFile]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **kwargs):
        htmlFiles = obj.getHtmlFiles()

        if not htmlFiles:
            raise RuntimeError(
                'No HTML reports found for this FASTQ object.'
            )

        paths = [Path(path).absolute() for path in dict.fromkeys(htmlFiles)]
        for path in paths:
            if not path.is_file() or path.stat().st_size == 0:
                raise RuntimeError('HTML report is missing or empty: {}'.format(path))
            if path.suffix.lower() != '.html':
                raise RuntimeError('Report is not an HTML file: {}'.format(path))
        for path in paths:
            webbrowser.open(path.as_uri())
        return []
