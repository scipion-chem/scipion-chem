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
from subprocess import Popen

from pwchem.objects import SetOfSmallMolecules
import pyworkflow.utils as pwutils
from pyworkflow.protocol.params import EnumParam, BooleanParam
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer

class SmallMoleculesViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer small molecules'
  _targets = [SetOfSmallMolecules]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)

  def _defineParams(self, form):
    form.addSection(label='Visualization of small molecules')
    form.addParam('displayPymolConformers', EnumParam,
                  choices=self.getChoicesConformers(), default=0,
                  label='Display small molecule conformers: ',
                  help='Docking results are grouped by their pocket, choose the one to visualize')
    form.addParam('displayPymolSingle', EnumParam,
                  choices=self.getChoicesSingle(), default=0,
                  label='Display single ligand: ',
                  help='Display this single ligand with the target')

  def getChoicesSingle(self):
    self.outputLigands = {}
    for mol in self.protocol:
        curMol = mol.clone()
        pName = curMol.getUniqueName()
        self.outputLigands[pName] = curMol

    outputLabels = list(self.outputLigands.keys())
    outputLabels.sort()
    return outputLabels

  def getChoicesConformers(self):
    self.outputLigandsC = {}
    for mol in self.protocol:
        curMol = mol.clone()
        pName = curMol.getMolBase()
        if pName in self.outputLigandsC:
          self.outputLigandsC[pName].append(curMol)
        else:
          self.outputLigandsC[pName] = [curMol]

    outputLabels = list(self.outputLigandsC.keys())
    outputLabels.sort()
    return outputLabels

  def _getVisualizeDict(self):
    return {
      'displayPymolSingle': self._viewSinglePymol,
      'displayPymolConformers': self._viewConformersPymol,
    }

  def _viewSinglePymol(self, e=None):
    ligandLabel = self.getEnumText('displayPymolSingle')
    mol = self.outputLigands[ligandLabel].clone()
    pmlFile = mol.getFileName()

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewConformersPymol(self, e=None):
    ligandLabel = self.getEnumText('displayPymolConformers')
    mols = self.outputLigandsC[ligandLabel]

    outDir = self.protocol.getSetDir()
    pmlFile = os.path.join(outDir, '{}.pml'.format(ligandLabel))
    with open(pmlFile, 'w') as f:
        for mol in mols:
            molFile, molName = os.path.abspath(mol.getFileName()), mol.getUniqueName()
            f.write('load {}, {}\nhide spheres, {}\nshow sticks, {}\n'.format(molFile, molName, molName, molName))

    pymolV = PyMolViewer(project=self.getProject())
    pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _validate(self):
    return []

