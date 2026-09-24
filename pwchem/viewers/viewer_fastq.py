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
import os
import webbrowser

from pyworkflow.viewer import Viewer, DESKTOP_TKINTER

from pwchem.objects import FastqFile


class FastqHtmlViewer(Viewer):
    """Viewer to open FASTQ-related HTML reports stored in a FastqFile object."""

    _label = 'FASTQ HTML reports viewer'
    _targets = [FastqFile]
    _environments = [DESKTOP_TKINTER]

    def _visualize(self, obj, **kwargs):
        htmlFiles = obj.getHtmlFiles()

        if not htmlFiles:
            raise RuntimeError(
                'No HTML reports found for this FASTQ object.'
            )

        for filePath in htmlFiles:
            if not os.path.exists(filePath):
                raise RuntimeError(
                    f'HTML report does not exist: {filePath}'
                )

            if not filePath.endswith('.html'):
                raise RuntimeError(
                    f'Report is not an HTML file: {filePath}'
                )

            webbrowser.open(
                'file://' + os.path.abspath(filePath)
            )

        return []

