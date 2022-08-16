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

from pyworkflow.protocol.params import EnumParam, LabelParam
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import showError

from pwchem.objects import SetOfSmallMolecules
from pwchem.viewers import PyMolViewer, BioinformaticsDataViewer
from pwchem.utils.utilsViewer import *
from pwchem.utils import runOpenBabel, mergePDBs, clean_PDB, natural_sort
from pwchem import Plugin as pwchemPlugin
from pwchem.protocols import ProtocolConsensusDocking, ProtocolLigandsFetching

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class SmallMoleculesViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer small molecules'
  _targets = [ProtocolConsensusDocking, ProtocolLigandsFetching, SetOfSmallMolecules]
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)
    self.singleLabels, self.singleLigandsDic = self.getChoices(vType=SINGLE)
    self.moleculeLabels, self.moleculeLigandsDic = self.getChoices(vType=MOLECULE)
    self.pocketLabels, self.pocketLigandsDic = self.getChoices(vType=POCKET)
    self.setLabels, self.setLigandsDic = self.getChoices(vType=SET)

  def defineParamsTable(self, form):
      form.addParam('displayTable', EnumParam,
                    choices=self.getChoices(vType=SET, pymol=False)[0], default=0,
                    label='Display ligands set and attributes in table format: ',
                    help='Display the ligands set in the set in table format with their respective attributes')


  def _defineParams(self, form):
    #Doking section
    if self.checkIfDocked():
        form.addSection(label='Docking view')
        group = form.addGroup('Each docked set')
        group.addParam('displayPymolSetDock', EnumParam,
                       choices=self.getChoices(vType=SET)[0], default=0,
                       label='Display docking on set result: ',
                       help='Display all ligands docked in this set')

        group = form.addGroup('Each docked ROI')
        group.addParam('displayPymolROIDock', EnumParam,
                       choices=self.getChoices(vType=POCKET)[0], default=0,
                       label='Display molecules in ROI: ',
                       help='Display all conformers and positions docked on this ROI')

        group = form.addGroup('Each docked molecule')
        group.addParam('displayPymolMoleculeDock', EnumParam,
                       choices=self.getChoices(vType=MOLECULE)[0], default=0,
                       label='Display one ligand type: ',
                       help='Display all conformers and positions of this molecule')

        group = form.addGroup('Each docked single')
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
    group = form.addGroup('Each set')
    group.addParam('displayPymolSet', EnumParam,
                   choices=self.getChoices(vType=SET)[0], default=0,
                   label='Display set: ',
                   help='Display all ligands in this set')

    group = form.addGroup('Each molecule')
    group.addParam('displayPymolMolecule', EnumParam,
                  choices=self.getChoices(vType=MOLECULE)[0], default=0,
                  label='Display small molecule conformers: ',
                  help='Docking results are grouped by their pocket, choose the one to visualize')

    group = form.addGroup('Each single')
    group.addParam('displayPymolSingle', EnumParam,
                  choices=self.getChoices(vType=SINGLE)[0], default=0,
                  label='Display single ligand: ',
                  help='Display this single ligand with the target')

    #Table section
    form.addSection(label='Table view')
    self.defineParamsTable(form)

  def updateLigandsDic(self, outputLigandsDic, molSet, vType, oLabelSet=None):
    if vType == SET:
        oLabel = oLabelSet
        outputLigandsDic[oLabel] = molSet
    else:
        for mol in molSet:
          curMol = mol.clone()
          if vType == SINGLE:
            oLabel = curMol.getUniqueName()
          elif vType == MOLECULE:
            oLabel = curMol.getMolName()
          elif vType == POCKET:
            oLabel = 'g_{}'.format(curMol.getGridId())

          if not oLabel in outputLigandsDic:
            outputLigandsDic[oLabel] = [curMol]
          else:
            outputLigandsDic[oLabel] += [curMol]
    return outputLigandsDic

  def getChoices(self, vType=SET, pymol=True):
    outputLigandsDic = {}
    if issubclass(type(self.protocol), SetOfSmallMolecules):
        '''If the viewer has been called for a SetOfMolecules'''
        molSet = self.protocol
        oLabel = 'outputSmallMolecules'
        outputLigandsDic = self.updateLigandsDic(outputLigandsDic, molSet, vType, oLabel)
    else:
        '''If the viewer has been called for a protocol with SetOfMolecules (can be several) as output'''
        for oAttr in self.protocol.iterOutputAttributes():
            if type(getattr(self.protocol, oAttr[0])) == SetOfSmallMolecules:
              molSet = getattr(self.protocol, oAttr[0])
              outputLigandsDic = self.updateLigandsDic(outputLigandsDic, molSet, vType, oAttr[0])

    outputLabels = list(outputLigandsDic.keys())
    outputLabels = natural_sort(outputLabels)
    if vType in SET and pymol and len(outputLabels) > 1:
        outputLabels = ['All'] + outputLabels
    return outputLabels, outputLigandsDic


  def _getVisualizeDict(self):
    return {
      #Docking
      'displayPymolSetDock': self._viewPymolSetDock,
      'displayPymolROIDock': self._viewPymolROIDock,
      'displayPymolMoleculeDock': self._viewPymolMoleculeDock,
      'displayPymolSingleDock': self._viewPymolSingleDock,
      'displayPymolPLIP': self._viewPymolPLIP,

      #Ligand
      'displayPymolSet': self._viewPymolSet,
      'displayPymolMolecule': self._viewPymolMolecule,
      'displayPymolSingle': self._viewPymolSingle,

      # Table
      'displayTable': self._viewSet,
    }


################# DOCKING VIEWS ###################

  def _viewPymolSetDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    if self.checkIfProtocol():
        ligandLabel = self.getEnumText('displayPymolSetDock')
    else:
        ligandLabel = 'All'

    if ligandLabel != 'All':
      pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))
      mols = sortMolsByUnique(self.setLigandsDic[ligandLabel])
      pmlStr = ''
      addTarget = True
      for mol in mols:
        pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
        addTarget = False

    else:
      pmlFile, pmlStr = os.path.join(pmlsDir, 'allSetDockedMolecules.pml'), ''
      addTarget = True
      for ligandLabel in self.setLabels:
        if ligandLabel != 'All':
          mols = sortMolsByUnique(self.setLigandsDic[ligandLabel])
          for mol in mols:
            pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
            addTarget = False
    writePmlFile(pmlFile, pmlStr)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolROIDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolROIDock')
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

    mols = sortMolsByUnique(self.pocketLigandsDic[ligandLabel])
    pmlStr, addTarget = '', True
    for mol in mols:
      pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
      addTarget = False
    writePmlFile(pmlFile, pmlStr)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolMoleculeDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolMoleculeDock')
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

    mols = sortMolsByUnique(self.moleculeLigandsDic[ligandLabel])
    pmlStr, addTarget = '', True
    for mol in mols:
      pmlStr += buildPMLDockingSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=addTarget)
      addTarget = False
    writePmlFile(pmlFile, pmlStr)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolSingleDock(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolSingleDock')
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))

    mol = self.singleLigandsDic[ligandLabel][0].clone()
    writePmlFile(pmlFile, buildPMLDockingSingleStr(self, mol, ligandLabel, addTarget=True, disable=False))

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolPLIP(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolPLIP')
    mol = self.singleLigandsDic[ligandLabel][0].clone()

    mergedPDB = self.createComplexPDB(self.protocol.getOriginalReceptorFile(), mol.getPoseFile(),
                                      os.path.join(pmlsDir, ligandLabel + '.pdb'))

    pwchemPlugin.runPLIP('-f {} -yt -o {}'.format(os.path.abspath(mergedPDB), ligandLabel),
                         cwd=os.path.abspath(pmlsDir))

    pmlFile, pmlFiles = '', []
    for file in os.listdir(os.path.abspath(os.path.join(pmlsDir, ligandLabel))):
        if file.endswith('.pse') and ligandLabel.upper().replace('-', '_') in file:
            pmlFiles.append(file)

    for file in pmlFiles:
        for ligName in ['UNK', 'UNL', 'LIG']: #typical ligand names
            if ligName in file:
                pmlFile = file
                break
    if pmlFile == '' and len(pmlFiles) > 0:
        pmlFile = pmlFiles[0]

    if pmlFile != '':
      pmlFile = os.path.join(os.path.abspath(pmlsDir), ligandLabel, pmlFile)
      pymolV = PyMolViewer(project=self.getProject())
      return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))
    else:
      showError('PLIP error', 'PLIP found no interactions in this docking position', self.getTkRoot())
      print('PLIP found no interactions in the docking position')

################### LIGANDS VIEWS #################

  def _viewPymolSet(self, e=None):
    pmlsDir = self.getPmlsDir()

    if self.checkIfProtocol():
        ligandLabel = self.getEnumText('displayPymolSet')
    else:
        ligandLabel = 'All'

    if ligandLabel != 'All':
      pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))
      mols = sortMolsByUnique(self.setLigandsDic[ligandLabel])
      pmlStr = ''
      for mol in mols:
        pmlStr += buildPMLFileNameSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=False)

    else:
      pmlFile, pmlStr = os.path.join(pmlsDir, 'allSetMolecules.pml'), ''
      for ligandLabel in self.setLabels:
        if ligandLabel != 'All':
          mols = sortMolsByUnique(self.setLigandsDic[ligandLabel])
          for mol in mols:
            pmlStr += buildPMLFileNameSingleStr(self, mol.clone(), mol.getUniqueName(), addTarget=False)
    writePmlFile(pmlFile, pmlStr)

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolSingle(self, e=None):
    ligandLabel = self.getEnumText('displayPymolSingle')
    mol = self.singleLigandsDic[ligandLabel][0].clone()
    pmlFile = mol.getFileName()

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def _viewPymolMolecule(self, e=None):
    ligandLabel = self.getEnumText('displayPymolMolecule')
    mols = sortMolsByUnique(self.moleculeLigandsDic[ligandLabel])

    outDir = self.getPmlsDir()
    pmlFile = os.path.join(outDir, '{}.pml'.format(ligandLabel))
    with open(pmlFile, 'w') as f:
        for mol in mols:
            molFile, molName = os.path.abspath(mol.getFileName()), mol.getUniqueName()
            f.write('load {}, {}\nhide spheres, {}\nshow sticks, {}\n'.format(molFile, molName, molName, molName))

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

    ###################   TABLE VIEW ####################

  def _viewSet(self, e=None):
    if self.checkIfProtocol():
      ligandLabel = self.getEnumText('displayTable')
      molSet = self.setLigandsDic[ligandLabel]
    else:
      molSet = self.protocol

    setV = BioinformaticsDataViewer(project=self.getProject())
    views = setV._visualize(molSet)
    return views


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


