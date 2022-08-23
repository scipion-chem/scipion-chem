# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
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

import os, pickle

from pyworkflow.protocol.params import LabelParam
import pyworkflow.viewer as pwviewer

from pwchem.viewers import PyMolViewer, PyMolView
from pwchem.protocols import ProtocolPharmacophoreFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_pharmacophore_generation import *
from pwchem.constants import PML_PHARM, SPHERE, LOAD_LIGAND, SPHERE_LIST, DISABLE_LIGAND
from pwchem.utils import getBaseFileName

colors = [(0, 0.9, 0), (0.9, 0, 0), (1, 0.9, 0), (0.5, 0, 1)]
colors_advanced = [(0, 1, 1), (1, 0.5, 0), (0, 0, 0.9), (0, 0.5, 1)]

feature_colors = {
    feat.lower(): color for feat, color in zip(FEATURE_LABELS_SIMPLE, colors)
}

feature_colors.update({
  feat.lower(): color for feat, color in zip(FEATURE_LABELS_ADVANCED , colors_advanced)
})


class PharmacophoreViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer pharmacophore'
  _targets = [ProtocolPharmacophoreFiltering]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    #Doking section
    form.addSection(label='Pharmacophore view')
    form.addParam('displayPharmacophore', LabelParam, label='Display docking on set result: ',
                  help='Display all ligands docked in this set')

  def _getVisualizeDict(self):
    return {
      #Docking
      'displayPharmacophore': self._viewPharm,
    }


################# DOCKING VIEWS ###################
  def _viewPharm(self, paramName=None):
      idxPath = os.path.abspath(self.protocol._getPath('cluster_index.pkl'))
      with open(idxPath, "rb") as clIdx:
        cluster_indices_sel = pickle.load(clIdx)

      cenPath = os.path.abspath(self.protocol._getPath('cluster_centers.pkl'))
      with open(cenPath, "rb") as clCenters:
        cluster_centers_sel = pickle.load(clCenters)

      # Load clusters
      cluster_indices_sel = {k.lower(): v for k, v in cluster_indices_sel.items()}
      cluster_centers_sel = {k.lower(): v for k, v in cluster_centers_sel.items()}

      pmlFile = self.buildPharmacophoreFile(cluster_indices_sel, cluster_centers_sel)
      return [PyMolView(pmlFile, cwd=self.protocol._getPath())]



#####################################################

  def buildPharmacophoreFile(self, idxDic, cenDic):
      spheresStr = ''
      for feature_type in idxDic.keys():
          featSpheres = ''
          centers = cenDic[feature_type]
          if centers:
              for i, loc in enumerate(centers):
                  sphere_radius = 1
                  if feature_type in feature_colors:
                      feature_color = feature_colors[feature_type]
                      featSpheres += SPHERE.format(*feature_color, *loc, sphere_radius)

              spheresStr += SPHERE_LIST.format(feature_type, featSpheres, feature_type, feature_type)

      protFile = os.path.abspath(self.protocol.inputLigands.get().getProteinFile())
      ligStr = self.getLigandsStr(self.getLigandsFiles(), disabled=True)

      pharStr = PML_PHARM.format(protFile, ligStr, spheresStr)

      pmlFile = os.path.abspath(self.protocol._getExtraPath('pharmacophore.pml'))
      with open(pmlFile, 'w') as f:
          f.write(pharStr)
      return pmlFile

  def getLigandsFiles(self):
      paramsFile = os.path.abspath(self.protocol._getExtraPath('inputParams.txt'))
      with open(paramsFile) as f:
          for line in f:
              if line.startswith('ligandFiles'):
                  files = line.split()[1:]
      return files

  def getLigandsStr(self, ligandFiles, disabled=False):
      ligStr = ''
      for file in ligandFiles:
          ligBase = getBaseFileName(file)
          ligStr += LOAD_LIGAND.format(file, ligBase)
          if disabled:
            ligStr += DISABLE_LIGAND.format(ligBase)
      return ligStr


