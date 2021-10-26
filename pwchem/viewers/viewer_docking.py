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

from pyworkflow.protocol.params import EnumParam, BooleanParam
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer
from pwchem.utils.utilsViewer import *

SINGLE, MOLECULE, POCKET = 'single', 'molecule', 'pocket'

class DockingViewer(pwviewer.ProtocolViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Docking viewer with Pymol'
    _targets = []
    _environments = [pwviewer.DESKTOP_TKINTER]

    def __init__(self, **args):
        pwviewer.ProtocolViewer.__init__(self, **args)
        self.singleLabels, self.singleLigandsDic = self.getChoices(type=SINGLE)
        self.moleculeLabels, self.moleculeLigandsDic = self.getChoices(type=MOLECULE)
        self.pocketLabels, self.pocketLigandsDic = self.getChoices(type=POCKET)

    def _defineParams(self, form):
        form.addSection(label='Visualize with Pymol')
        group = form.addGroup('Pocket ligands')
        group.addParam('displayPymolPocket', EnumParam,
                       choices=self.getChoices(type=POCKET)[0], default=0,
                       label='Display docking on pocket result: ',
                       help='Display all ligands docked in this pocket')

        group = form.addGroup('Each molecule')
        group.addParam('displayPymolMolecule', EnumParam,
                       choices=self.getChoices(type=MOLECULE)[0], default=0,
                       label='Display one ligand type: ',
                       help='Display all conformers and positions of this molecule')

        group = form.addGroup('Single Ligand')
        group.addParam('displayPymolSingle', EnumParam,
                       choices=self.getChoices(type=SINGLE)[0], default=0,
                       label='Display single ligand: ',
                       help='Display this single ligand with the target')

    def getChoices(self, type=POCKET, pymol=True):
        outputLigandsDic = {}
        for oAttr in self.protocol.iterOutputAttributes():
          if 'outputSmallMolecules' in oAttr[0]:
            if type == POCKET:
              oLabel = oAttr[0]
            molSet = getattr(self.protocol, oAttr[0])
            for mol in molSet:
              curMol = mol.clone()
              if type == SINGLE:
                oLabel = curMol.getUniqueName()
              elif type == MOLECULE:
                oLabel = curMol.getMolBase()
              if not oLabel in outputLigandsDic:
                outputLigandsDic[oLabel] = [curMol]
              else:
                outputLigandsDic[oLabel] += [curMol]

        outputLabels = list(outputLigandsDic.keys())
        outputLabels.sort()
        if type == POCKET and pymol and len(outputLabels) > 1:
          outputLabels = ['All'] + outputLabels
        return outputLabels, outputLigandsDic

    def _getVisualizeDict(self):
        return {'displayPymolSingle': self._viewSinglePymol,
                'displayPymolMolecule': self._viewMoleculePymol,
                'displayPymolPocket': self._viewPocketPymol,
                }

    def _viewSinglePymol(self, e=None):
        ligandLabel = self.getEnumText('displayPymolSingle')
        pmlsDir = getPmlsDir(self.protocol)
        pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

        mol = self.singleLigandsDic[ligandLabel][0].clone()
        writePmlFile(pmlFile, buildPMLDockingSingleStr(self, mol, ligandLabel, addTarget=True, disable=False))

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewMoleculePymol(self, e=None):
        ligandLabel = self.getEnumText('displayPymolMolecule')
        pmlsDir = getPmlsDir(self.protocol)
        pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

        mols = self.moleculeLigandsDic[ligandLabel]
        mols = sortMolsByUnique(mols)
        pmlStr, addTarget = '', True
        for mol in mols:
          pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
          addTarget = False
        writePmlFile(pmlFile, pmlStr)

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    def _viewPocketPymol(self, e=None):
        pmlsDir = getPmlsDir(self.protocol)
        ligandLabel = self.getEnumText('displayPymolPocket')
        if ligandLabel != 'All':
          pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))
          mols = self.pocketLigandsDic[ligandLabel]
          mols, pmlStr = sortMolsByUnique(mols), ''
          addTarget=True
          for mol in mols:
            pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
            addTarget = False

        else:
          pmlFile, pmlStr = os.path.join(pmlsDir, 'allDockedMolecules.pml'), ''
          addTarget = True
          for ligandLabel in self.pocketLabels:
            if ligandLabel != 'All':
              mols = self.pocketLigandsDic[ligandLabel]
              mols = sortMolsByUnique(mols)
              for mol in mols:
                pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
                addTarget = False
        writePmlFile(pmlFile, pmlStr)

        pymolV = PyMolViewer(project=self.getProject())
        pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))



