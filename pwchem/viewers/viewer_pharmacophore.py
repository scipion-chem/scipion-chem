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

import os

from pyworkflow.protocol.params import *
import pyworkflow.viewer as pwviewer

from pwchem.viewers import PyMolView
from pwchem.objects import PharmacophoreChem
from pwchem.protocols import ProtocolPharmacophoreGeneration, ProtocolPharmacophoreFiltering
from pwchem.protocols.VirtualDrugScreening.protocol_pharmacophore_generation import *
from pwchem.constants import *
from pwchem.utils import getBaseName

colors = [(0, 0.9, 0), (0.9, 0, 0), (1, 0.9, 0), (0.5, 0, 1)]
colors_advanced = [(0, 1, 1), (1, 0.5, 0), (0, 0, 0.9), (0, 0.5, 1)]

feature_colors = {
    feat: color for feat, color in zip(FEATURE_LABELS_SIMPLE, colors)
}

feature_colors.update({
  feat: color for feat, color in zip(FEATURE_LABELS_ADVANCED , colors_advanced)
})


class PharmacophoreViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer pharmacophore'
  _targets = [PharmacophoreChem]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    #Doking section
    form.addSection(label='Pharmacophore view')
    form.addParam('transparency', FloatParam, label='Transparency: ', default=0.5,
                  help='Transparency (alpha) for the pharmacophore spheres (1 means no transparency, 0 transparent)')
    form.addParam('displayPharmacophore', LabelParam, label='Display docking on set result: ',
                  help='Display all ligands docked in this set')

  def _getVisualizeDict(self):
    return {
      #Docking
      'displayPharmacophore': self._viewPharm,
    }


################# DOCKING VIEWS ###################
  def _viewPharm(self, paramName=None):
      pharm_centers, pharm_radii = self.getPharmPoints()
      pharmStr = self.buildPharmacophoreStr(pharm_centers, pharm_radii)
      pmlFile = self.writePharmacophoreFile(pharmStr)
      return [PyMolView(pmlFile, cwd=self.getPharmacophoreObject().getSetDir())]

#####################################################

  def getPharmPoints(self):
      inPharm = self.getPharmacophoreObject()
      pharm_centers, pharm_radii = {}, {}
      for featObj in inPharm:
          featType = featObj.getType()
          if featType in pharm_centers:
            pharm_centers[featType].append(featObj.getCoords())
            pharm_radii[featType].append(featObj.getRadius())
          else:
            pharm_centers[featType] = [featObj.getCoords()]
            pharm_radii[featType] = [featObj.getRadius()]
      return pharm_centers, pharm_radii

  def buildPharmacophoreStr(self, cenDic, radDic):
      spheresStr = ''
      for feature_type in cenDic:
          featSpheres = ''
          centers = cenDic[feature_type]
          radii = radDic[feature_type]
          if centers:
              for i, loc in enumerate(centers):
                  if feature_type in feature_colors:
                      feature_color = feature_colors[feature_type]
                      featSpheres += SPHERE.format(self.transparency.get(), *feature_color, *loc, radii[i])

              spheresStr += SPHERE_LIST.format(feature_type, featSpheres, feature_type, feature_type)

      protFile = self.getPharmacophoreObject().getProteinFile()
      if not protFile:
          protStr = ''
      else:
          protStr = 'load {}, receptor'.format(os.path.abspath(protFile))

      ligStr = self.getLigandsStr()

      pharStr = PML_PHARM.format(protStr, ligStr, spheresStr)
      return pharStr

  def writePharmacophoreFile(self, pharmStr):
      pmlFile = os.path.abspath(self.getPharmacophoreObject().getSetDir('pharmacophore.pml'))
      with open(pmlFile, 'w') as f:
          f.write(pharmStr)
      return pmlFile

  def getPharmacophoreObject(self):
      return self.protocol

  def getLigandsStr(self, disabled=True):
      return ''


class GeneratePharmacophoreViewer(PharmacophoreViewer):
    _label = 'Viewer generate pharmacophore'
    _targets = [ProtocolPharmacophoreGeneration]
    _environments = [pwviewer.DESKTOP_TKINTER]

    def getLigandsFiles(self):
      paramsFile = os.path.abspath(self.protocol._getExtraPath('inputParams.txt'))
      with open(paramsFile) as f:
        for line in f:
          if line.startswith('ligandFiles'):
            files = line.split()[1:]
      return files

    def getLigandsStr(self, disabled=True):
      ligStr = ''
      for file in self.getLigandsFiles():
        ligBase = getBaseName(file)
        ligStr += LOAD_LIGAND.format(file, ligBase)
        if disabled:
          ligStr += DISABLE_LIGAND.format(ligBase)
      return ligStr

    def getPharmacophoreObject(self):
      return self.protocol.outputPharmacophore


class FilterPharmacophoreViewer(PharmacophoreViewer):
  _label = 'Viewer pharmacophore filtering'
  _targets = [ProtocolPharmacophoreFiltering]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def getLigandsFiles(self):
    files = []
    for mol in self.protocol.outputSmallMolecules:
        files.append(os.path.abspath(mol.getPoseFile()))
    return files

  def getLigandsStr(self, disabled=True):
    ligStr = ''
    for file in self.getLigandsFiles():
      ligBase = getBaseName(file)
      ligStr += LOAD_LIGAND.format(file, ligBase)
      if disabled:
        ligStr += DISABLE_LIGAND.format(ligBase)
    return ligStr

  def getPharmacophoreObject(self):
    return self.protocol.inputPharmacophore.get()


