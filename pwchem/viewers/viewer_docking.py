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

from pyworkflow.protocol.params import EnumParam, LabelParam
import pyworkflow.viewer as pwviewer
from pwchem.utils.utilsViewer import *
from pwchem.utils import runOpenBabel, mergePDBs, clean_PDB
from pyworkflow.viewer import DESKTOP_TKINTER
from pwchem.protocols import ProtocolConsensusDocking
from pwchem import Plugin as pwchemPlugin
from pyworkflow.gui.dialog import showError
from pwchem.viewers import BioinformaticsDataViewer

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class DockingViewer(pwviewer.ProtocolViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Docking viewer with Pymol'
    _targets = []
    _environments = [pwviewer.DESKTOP_TKINTER]

    def __init__(self, **args):
        pwviewer.ProtocolViewer.__init__(self, **args)
        self.singleLabels, self.singleLigandsDic = self.getChoices(vType=SINGLE)
        self.moleculeLabels, self.moleculeLigandsDic = self.getChoices(vType=MOLECULE)
        self.pocketLabels, self.pocketLigandsDic = self.getChoices(vType=POCKET)
        self.setLabels, self.setLigandsDic = self.getChoices(vType=SET)

    def _defineParams(self, form):
        form.addSection(label='Visualize with Pymol')
        group = form.addGroup('Pocket ligands')
        group.addParam('displayPymolPocket', EnumParam,
                       choices=self.getChoices(vType=POCKET)[0], default=0,
                       label='Display docking on pocket result: ',
                       help='Display all ligands docked in this pocket')

        group = form.addGroup('Each molecule')
        group.addParam('displayPymolMolecule', EnumParam,
                       choices=self.getChoices(vType=MOLECULE)[0], default=0,
                       label='Display one ligand vType: ',
                       help='Display all conformers and positions of this molecule')

        group = form.addGroup('Single Ligand')
        group.addParam('displayPymolSingle', EnumParam,
                       choices=self.getChoices(vType=SINGLE)[0], default=0,
                       label='Display single ligand: ',
                       help='Display this single ligand with the target')

        form.addSection(label='Visualize with PLIP')
        form.addParam('displayPymolPLIP', EnumParam,
                       choices=self.getChoices(vType=SINGLE)[0], default=0,
                       label='Display ligand interactions: ',
                       help='Display this single ligand with the binding site and interactions')

        form.addSection(label='Table view')
        form.addParam('displayTable', EnumParam,
                      choices=self.getChoices(vType=SET)[0], default=0,
                      label='Display ligands set and attributes in table format: ',
                      help='Display the chosen ligands set in the set in table format with their respective attributes')


    def getChoices(self, vType=POCKET, pymol=True):
        outputLigandsDic = {}
        for oAttr in self.protocol.iterOutputAttributes():
          if 'outputSmallMolecules' in oAttr[0]:
            if vType == POCKET or vType == SET:
              oLabel = oAttr[0]
            molSet = getattr(self.protocol, oAttr[0])
            if vType == SET:
              outputLigandsDic[oLabel] = molSet
            else:
              for mol in molSet:
                curMol = mol.clone()
                if vType == SINGLE:
                  oLabel = curMol.getUniqueName()
                elif vType == MOLECULE:
                  oLabel = curMol.getMolBase()
                if not oLabel in outputLigandsDic:
                  outputLigandsDic[oLabel] = [curMol]
                else:
                  outputLigandsDic[oLabel] += [curMol]

        outputLabels = list(outputLigandsDic.keys())
        outputLabels.sort()
        if vType == POCKET and pymol and len(outputLabels) > 1:
          outputLabels = ['All'] + outputLabels
        return outputLabels, outputLigandsDic

    def _getVisualizeDict(self):
        return {'displayPymolSingle': self._viewSinglePymol,
                'displayPymolMolecule': self._viewMoleculePymol,
                'displayPymolPocket': self._viewPocketPymol,
                'displayPymolPLIP': self._viewPLIPPymol,
                'displayTable': self._viewSet,
                }

    def _viewSet(self, e=None):
        ligandLabel = self.getEnumText('displayTable')
        mols = self.setLigandsDic[ligandLabel]

        setV = BioinformaticsDataViewer(project=self.getProject())
        views = setV._visualize(mols)
        views[0].show()

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

    def _viewPLIPPymol(self, e=None):
      ligandLabel = self.getEnumText('displayPymolPLIP')
      pmlsDir = getPmlsDir(self.protocol)
      mol = self.singleLigandsDic[ligandLabel][0].clone()

      mergedPDB = self.createComplexPDB(self.protocol.getOriginalReceptorFile(), mol.getPoseFile(),
                                        os.path.join(pmlsDir, ligandLabel+'.pdb'))

      pwchemPlugin.runPLIP('-f {} -yt -o {}'.format(os.path.abspath(mergedPDB), ligandLabel),
                           cwd=os.path.abspath(pmlsDir))

      pmlFile = ''
      for file in os.listdir(os.path.abspath(os.path.join(pmlsDir, ligandLabel))):
          if file.endswith('.pse') and self.typicalLigNames(file):
              pmlFile = file

      if pmlFile != '':
          pmlFile = os.path.join(os.path.abspath(pmlsDir), ligandLabel, pmlFile)
          pymolV = PyMolViewer(project=self.getProject())
          pymolV.visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))
      else:
          showError('PLIP error', 'PLIP found no interactions in this docking position', self.getTkRoot())
          print('PLIP found no interactions in the docking position')

    def typicalLigNames(self, file):
        ligNames = ['UNK', 'UNL', 'LIG', 'X']
        for ln in ligNames:
            if ln in file:
                return True
        return False

    def createComplexPDB(self, receptorFile, molFile, outPath):
        outBase, ext = os.path.splitext(molFile)
        auxPath = '/tmp/{}.pdb'.format(os.path.basename(outBase))
        if not molFile.endswith('.pdb'):
            outFile = os.path.abspath(outBase + '.pdb')
            runOpenBabel(self.protocol, '-i{} {} -opdb -O {}'.format(ext[1:], molFile, outFile),
                         cwd=self.protocol._getExtraPath(), popen=True)
            molFile = outFile

        if not receptorFile.endswith('.pdb'):
            outBase, ext = os.path.splitext(receptorFile)
            outFile = os.path.abspath(outBase + '.pdb')
            runOpenBabel(self.protocol, '-i{} {} -opdb -O {}'.format(ext[1:], receptorFile, outFile),
                         cwd=self.protocol._getExtraPath(), popen=True)
            receptorFile = outFile

        mergePDBs(receptorFile, molFile, auxPath, hetatm2=True)
        clean_PDB(auxPath, outPath, waters=True, HETATM=False)
        return outPath



class ProtConsensusDockingViewer(DockingViewer):
    """ Visualize the output of protocol autodock """
    _label = 'Viewer autodock docking'
    _targets = [ProtocolConsensusDocking]
    _environments = [DESKTOP_TKINTER]

    def __init__(self, **args):
        DockingViewer.__init__(self, **args)

    def getChoices(self, vType=POCKET, pymol=True):
        outputLigandsDic, oLabel = {}, None
        for oAttr in self.protocol.iterOutputAttributes():
          if 'outputSmallMolecules' in oAttr[0]:
            if vType == POCKET:
              oLabel = oAttr[0]
            molSet = getattr(self.protocol, oAttr[0])
            for mol in molSet:
              curMol = mol.clone()
              if vType == SINGLE:
                oLabel = curMol.getUniqueName()
              elif vType == MOLECULE and oAttr[0] != 'outputSmallMoleculesAll':
                oLabel = curMol.getMolBase()

              if oLabel is not None:
                if not oLabel in outputLigandsDic:
                  outputLigandsDic[oLabel] = [curMol]
                else:
                  outputLigandsDic[oLabel] += [curMol]

        outputLabels = list(outputLigandsDic.keys())
        outputLabels.sort()
        return outputLabels, outputLigandsDic

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

