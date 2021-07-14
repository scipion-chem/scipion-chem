# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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

import pyworkflow.viewer as pwviewer
import pwem.viewers.views as views
import pwem.viewers.showj as showj
import pwchem.objects

class SetOfDatabaseIDView(views.ObjectView):
    """ Customized ObjectView for SetOfDatabaseID. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.MODE: 'metadata',
                             showj.RENDER: '_PDBLigandImage'}
        defaultViewParams.update(viewParams)
        views.ObjectView.__init__(self, project, inputid, path, other,
                                  defaultViewParams, **kwargs)

class BioinformaticsDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        pwchem.objects.SetOfDatabaseID,
        pwchem.objects.ProteinSequenceFile,
        pwchem.objects.NucleotideSequenceFile,
        pwchem.objects.SetOfSmallMolecules,
        pwchem.objects.SetOfBindingSites
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)

        # For now handle both types of SetOfTiltSeries together
        if issubclass(cls, pwchem.objects.SetOfDatabaseID):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfSmallMolecules):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfBindingSites):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.ProteinSequenceFile):
            views.append(self.textView([obj.getFileName()]))
        elif issubclass(cls, pwchem.objects.NucleotideSequenceFile):
            views.append(self.textView([obj.getFileName()]))

        return views
