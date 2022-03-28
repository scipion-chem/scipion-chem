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

from pwchem.objects import SetOfSmallMolecules
from pyworkflow.protocol.params import EnumParam, LabelParam
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer, BioinformaticsDataViewer
from pwchem.utils.utilsViewer import sortMolsByUnique, buildPMLDockingSingleStr, writePmlFile, getPmlsDir
from pwchem.utils import runOpenBabel, mergePDBs, clean_PDB
from pyworkflow.gui.dialog import showError
from pwchem import Plugin as pwchemPlugin

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class SmallMoleculesViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer small molecules'
  _targets = [SetOfSmallMolecules]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)
    self.singleLabels, self.singleLigandsDic = self.getChoices(vType=SINGLE)
    self.moleculeLabels, self.moleculeLigandsDic = self.getChoices(vType=MOLECULE)
    self.pocketLabels, self.pocketLigandsDic = self.getChoices(vType=POCKET)
    self.setLabels, self.setLigandsDic = self.getChoices(vType=SET)

  def defineParamsTableProtocol(self, form):
      form.addParam('displayTableProtocol', EnumParam,
                    choices=self.getChoices(vType=SET)[0], default=0,
                    label='Display ligands set and attributes in table format: ',
                    help='Display the ligands set in the set in table format with their respective attributes')

  def defineParamsTableSet(self, form):
      form.addParam('displayTableSet', LabelParam,
                    label='Display ligands set and attributes in table format: ',
                    help='Display the ligands set in the set in table format with their respective attributes')

  def _defineParams(self, form):
    #Doking section
    if self.checkIfDocked():
        form.addSection(label='Docking view')
        group = form.addGroup('Each set')
        if self.checkIfProtocol():
            group.addParam('displayPymolPocket', EnumParam,
                           choices=self.getChoices(vType=POCKET)[0], default=0,
                           label='Display docking on set result: ',
                           help='Display all ligands docked in this set')
        else:
            group.addParam('displayPymolPocketSet', LabelParam,
                           label='Display docking on set result: ',
                           help='Display all ligands docked in this set')

        group = form.addGroup('Each molecule')
        group.addParam('displayPymolMoleculeDock', EnumParam,
                       choices=self.getChoices(vType=MOLECULE)[0], default=0,
                       label='Display one ligand type: ',
                       help='Display all conformers and positions of this molecule')

        group = form.addGroup('Each position')
        group.addParam('displayPymolSingleDock', EnumParam,
                       choices=self.getChoices(vType=SINGLE)[0], default=0,
                       label='Display single ligand: ',
                       help='Display this single ligand with the target')

        group = form.addGroup('Visualize with PLIP')
        group.addParam('displayPymolPLIP', EnumParam,
                       choices=self.getChoices(vType=SINGLE)[0], default=0,
                       label='Display ligand interactions: ',
                       help='Display this single ligand with the binding site and interactions')

    #Molecules section
    form.addSection(label='Small molecules view')
    form.addParam('displayPymolMolecules', EnumParam,
                  choices=self.getChoices(vType=MOLECULE)[0], default=0,
                  label='Display small molecule conformers: ',
                  help='Docking results are grouped by their pocket, choose the one to visualize')
    form.addParam('displayPymolSingle', EnumParam,
                  choices=self.getChoices(vType=SINGLE)[0], default=0,
                  label='Display single ligand: ',
                  help='Display this single ligand with the target')

    #Table section
    form.addSection(label='Table view')
    if not self.checkIfProtocol():
        self.defineParamsTableSet(form)
    else:
        self.defineParamsTableProtocol(form)


  def getChoices(self, vType=POCKET, pymol=True):
    outputLigandsDic = {}
    if issubclass(type(self.protocol), SetOfSmallMolecules):
        '''If the viewer has been called for a SetOfMolecules'''
        molSet = self.protocol
        oLabel = 'outputSmallMolecules'
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
    else:
        '''If the viewer has been called for a protocol with SetOfMolecules (can be several) as output'''
        for oAttr in self.protocol.iterOutputAttributes():
          if type(getattr(self.protocol, oAttr[0])) == SetOfSmallMolecules:
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
    return {
      #Docking
      'displayPymolPocketSet': self._viewPocketPymolDock,
      'displayPymolPocket': self._viewPocketPymolDock,
      'displayPymolMoleculeDock': self._viewMoleculePymolDock,
      'displayPymolSingleDock': self._viewSinglePymolDock,
      'displayPymolPLIP': self._viewPLIPPymol,

      #Table
      'displayTableSet': self._viewSet,

      #Ligand
      'displayPymolSingle': self._viewSinglePymol,
      'displayPymolMolecules': self._viewMoleculesPymol,
      'displayTableProtocol': self._viewSet,
    }


################# DOCKING VIEWS ###################

  def _viewPocketPymolDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    if self.checkIfProtocol():
        ligandLabel = self.getEnumText('displayPymolPocket')
    else:
        ligandLabel = 'All'

    if ligandLabel != 'All':
      pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))
      mols = self.pocketLigandsDic[ligandLabel]
      mols, pmlStr = sortMolsByUnique(mols), ''
      addTarget = True
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
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewMoleculePymolDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolMoleculeDock')
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

    mols = self.moleculeLigandsDic[ligandLabel]
    mols = sortMolsByUnique(mols)
    pmlStr, addTarget = '', True
    for mol in mols:
      pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
      addTarget = False
    writePmlFile(pmlFile, pmlStr)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewSinglePymolDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolSingleDock')
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

    mol = self.singleLigandsDic[ligandLabel][0].clone()
    writePmlFile(pmlFile, buildPMLDockingSingleStr(self, mol, ligandLabel, addTarget=True, disable=False))

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPLIPPymol(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolPLIP')
    mol = self.singleLigandsDic[ligandLabel][0].clone()

    mergedPDB = self.createComplexPDB(self.protocol.getOriginalReceptorFile(), mol.getPoseFile(),
                                      os.path.join(pmlsDir, ligandLabel + '.pdb'))

    pwchemPlugin.runPLIP('-f {} -yt -o {}'.format(os.path.abspath(mergedPDB), ligandLabel),
                         cwd=os.path.abspath(pmlsDir))

    pmlFile = ''
    for file in os.listdir(os.path.abspath(os.path.join(pmlsDir, ligandLabel))):
      if file.endswith('.pse') and ligandLabel.upper().replace('-', '_') in file:
        pmlFile = file

    if pmlFile != '':
      pmlFile = os.path.join(os.path.abspath(pmlsDir), ligandLabel, pmlFile)
      pymolV = PyMolViewer(project=self.getProject())
      return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))
    else:
      showError('PLIP error', 'PLIP found no interactions in this docking position', self.getTkRoot())
      print('PLIP found no interactions in the docking position')


###################   TABLE VIEW ####################

  def _viewSet(self, e=None):
    if self.checkIfProtocol():
        ligandLabel = self.getEnumText('displayTableProtocol')
        molSet = self.setLigandsDic[ligandLabel]
    else:
        molSet = self.protocol

    setV = BioinformaticsDataViewer(project=self.getProject())
    views = setV._visualize(molSet)
    return views


################### LIGANDS VIEWS #################

  def _viewSinglePymol(self, e=None):
    ligandLabel = self.getEnumText('displayPymolSingle')
    mol = self.singleLigandsDic[ligandLabel][0].clone()
    pmlFile = mol.getFileName()

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewMoleculesPymol(self, e=None):
    ligandLabel = self.getEnumText('displayPymolMolecules')
    mols = self.moleculeLigandsDic[ligandLabel]

    outDir = self.getPmlsDir()
    pmlFile = os.path.join(outDir, '{}.pml'.format(ligandLabel))
    with open(pmlFile, 'w') as f:
        for mol in mols:
            molFile, molName = os.path.abspath(mol.getFileName()), mol.getUniqueName()
            f.write('load {}, {}\nhide spheres, {}\nshow sticks, {}\n'.format(molFile, molName, molName, molName))

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))


#################### UTILS ###########################33

  def checkIfDocked(self):
      if not self.checkIfProtocol():
          molSet = self.protocol
      else:
          for oAttr in self.protocol.iterOutputAttributes():
              if type(getattr(self.protocol, oAttr[0])) == SetOfSmallMolecules:
                  molSet = getattr(self.protocol, oAttr[0])
      return molSet.isDocked()

  def checkIfProtocol(self):
      if issubclass(type(self.protocol), SetOfSmallMolecules):
          return False
      else:
          return True

  def typicalLigNames(self, file):
    ligNames = ['UNK', 'UNL', 'LIG', 'X', 'LG']
    for ln in ligNames:
      if ln in file:
        return True
    return False

  def getPmlsDir(self):
    if self.checkIfProtocol():
      pmlsDir = getPmlsDir(self.protocol)
    else:
      pmlsDir = os.path.join(self.protocol.getSetDir(), 'pmls')
      if not os.path.exists(pmlsDir):
          os.mkdir(pmlsDir)
    return pmlsDir

  def createComplexPDB(self, receptorFile, molFile, outPath):
    receptorFile, molFile = os.path.abspath(receptorFile), os.path.abspath(molFile)
    outBase, ext = os.path.splitext(molFile)
    auxPath = '/tmp/{}.pdb'.format(os.path.basename(outBase))
    if not molFile.endswith('.pdb'):
      outFile = os.path.abspath(outBase + '.pdb')
      runOpenBabel(self.protocol, '-i{} {} -opdb -O {}'.format(ext[1:], molFile, outFile),
                   popen=True)
      molFile = outFile

    if not receptorFile.endswith('.pdb'):
      outBase, ext = os.path.splitext(receptorFile)
      outFile = os.path.abspath(outBase + '.pdb')
      runOpenBabel(self.protocol, '-i{} {} -opdb -O {}'.format(ext[1:], receptorFile, outFile),
                   popen=True)
      receptorFile = outFile

    mergePDBs(receptorFile, molFile, auxPath, hetatm2=True)
    clean_PDB(auxPath, outPath, waters=True, HETATM=False)
    return outPath

  def _validate(self):
    return []

