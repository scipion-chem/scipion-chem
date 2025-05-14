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

from matplotlib import pyplot as plt
import numpy as np

from pyworkflow.protocol import params
import pyworkflow.viewer as pwviewer
from pyworkflow.gui.dialog import showError, askYesNoCancel

from pwem.viewers.viewer_chimera import Chimera, ChimeraView
from pwem.viewers.mdviewer.viewer import MDViewer

from pwchem.objects import SetOfSmallMolecules, SmallMoleculesLibrary
from pwchem.viewers import PyMolViewer, BioinformaticsDataViewer
from pwchem.utils.utilsViewer import *
from pwchem.utils import runOpenBabel, mergePDBs, cleanPDB, natural_sort
from pwchem import Plugin as pwchemPlugin
from pwchem.protocols import ProtocolConsensusDocking, ProtocolLigandsFetching

PYMOL, CHIMERAX, VIEWDOCKX = 0, 1, 2
SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

def chimeraInstalled():
  return Chimera.getHome() and os.path.exists(Chimera.getProgram())

CHIMERA_ERROR = 'Chimera program is not found were it was expected: \n\n{}\n\n' \
                'Either install ChimeraX in this path or install our ' \
                'scipion-em-chimera plugin'.format(Chimera.getProgram())

class SmallMoleculesViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer small molecules'
  _targets = [ProtocolConsensusDocking, ProtocolLigandsFetching, SetOfSmallMolecules]
  _environments = [pwviewer.DESKTOP_TKINTER]
  _viewerOptions = ['PyMol', 'ChimeraX', 'ViewDockX']

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)
    self.protocolObject = self.protocol
    self.singleLabels, self.singleLigandsDic = self.getChoices(vType=SINGLE)
    self.moleculeLabels, self.moleculeLigandsDic = self.getChoices(vType=MOLECULE)
    self.pocketLabels, self.pocketLigandsDic = self.getChoices(vType=POCKET)
    self.setLabels, self.setLigandsDic = self.getChoices(vType=SET)

  def getViewerOptions(self):
    if type(self.protocol).__name__ == 'ProtChemVinaDocking':
      return self._viewerOptions
    return self._viewerOptions[:-1]

  def defineParamsTable(self, form):
      form.addParam('displayTable', params.EnumParam,
                    choices=self.getChoices(vType=SET, pymol=False)[0], default=0,
                    label='Display ligands set and attributes in table format: ',
                    help='Display the ligands set in the set in table format with their respective attributes')

  def getOutputSetLabels(self):
    if not self.setLabels:
      for oAttr in self.protocolObject.iterOutputAttributes():
        if type(getattr(self.protocol, oAttr[0])) == SetOfSmallMolecules:
          self.setLabels.append(oAttr[0])
    return self.setLabels

  def _defineParams(self, form):
    #Doking section
    if self.checkIfDocked():
        form.addSection(label='Docking view')
        group = form.addGroup('Visualize docking poses with:')
        group.addParam('displaySoftware', params.EnumParam,
                       choices=self.getViewerOptions(), default=PYMOL,
                       label='Display docking poses with: ',
                       help='Display the selected group of docking poses with which software. '
                            'Available: PyMol, ChimeraX and viewDockX (ChimeraX with stored attributes)')

        group = form.addGroup('Display docking poses')
        group.addParam('displaySetDock', params.EnumParam,
                       choices=self.setLabels, default=0,
                       label='Display docking poses in set: ',
                       help='Display all ligands docked in this set')

        group.addParam('displayROIDock', params.EnumParam,
                       choices=self.pocketLabels, default=0,
                       label='Display docking poses in ROI: ',
                       help='Display all conformers and positions docked on this ROI')

        group.addParam('displayMoleculeDock', params.EnumParam,
                       choices=self.moleculeLabels, default=0,
                       label='Display docking poses of molecule: ',
                       help='Display all conformers and positions of this molecule')

        group.addParam('displaySingleDock', params.EnumParam,
                       choices=self.singleLabels, default=0,
                       label='Display docking pose of ligand: ',
                       help='Display this single ligand with the target')

        group = form.addGroup('Visualize with PLIP')
        group.addParam('displayPymolPLIP', params.EnumParam,
                       choices=self.singleLabels, default=0,
                       label='Display ligand interactions: ',
                       help='Display this single ligand with the binding site and interactions')

    #Molecules section
    form.addSection(label='Small molecules view')
    group = form.addGroup('Visualize molecules with:')
    group.addParam('displaySoftwareMol', params.EnumParam,
                   choices=self._viewerOptions[:2], default=PYMOL,
                   label='Display docking poses with: ',
                   help='Display the selected group of docking poses with which software. '
                        'Available: PyMol, ChimeraX and viewDockX (ChimeraX with stored attributes)')

    group = form.addGroup('Display molecules')
    group.addParam('displaySet', params.EnumParam,
                   choices=self.setLabels, default=0,
                   label='Display molecules in set: ',
                   help='Display all ligands in this set')

    group.addParam('displayMolecule', params.EnumParam,
                  choices=self.moleculeLabels, default=0,
                  label='Display molecule conformers: ',
                  help='Display all conformers saved for this molecule')

    group.addParam('displaySingle', params.EnumParam,
                  choices=self.singleLabels, default=0,
                  label='Display single molecule: ',
                  help='Display this single ligand with the target')

    #Table section
    form.addSection(label='Table view')
    self.defineParamsTable(form)

  def getChoices(self, vType=SET, pymol=True):
    '''Return based on the subgroups of the output set(s):
        outputLabels: list of labels with the choices of subgroups that can be made for selection vType
                      (single, molecule, pocket or set)
        outputLigandsDic: dictionary containing the items contained in each of the subgroups in the labels
                          the dictionary is  separated by each of the output sets:
                          {outputSetLabel1: {outputSubGroupLabel1: [indexes]}, }
    '''
    outputLigandsDic, outputLabels = {}, []
    if issubclass(type(self.protocol), SetOfSmallMolecules):
      '''If the viewer has been called for a SetOfMolecules'''
      molSet = self.protocol
      setLabel = 'outputSmallMolecules'
      outputLigandsDic[setLabel] = molSet.getGroupIndexes()[vType]
      outputLabels = [setLabel]

    else:
      '''If the viewer has been called for a protocol with SetOfMolecules (can be several) as output'''
      molSets = self.getOutputMolSets()
      for setLabel, molSet in molSets.items():
          molSet = getattr(self.protocol, setLabel)
          curLigDic = molSet.getGroupIndexes(setLabel)[vType]

          outputLabels += list(curLigDic.keys())
          outputLigandsDic[setLabel] = curLigDic

    outputLabels = list(set(outputLabels))
    outputLabels = natural_sort(outputLabels)
    if vType in SET and pymol and len(outputLabels) > 1:
        outputLabels = ['All'] + outputLabels
    return outputLabels, outputLigandsDic


  def _getVisualizeDict(self):
    return {
      #Docking
      'displaySetDock': self._viewSetDock,
      'displayROIDock': self._viewROIDock,
      'displayMoleculeDock': self._viewMoleculeDock,
      'displaySingleDock': self._viewSingleDock,
      'displayPymolPLIP': self._viewPymolPLIP,

      #Ligand
      'displaySet': self._viewSet,
      'displayMolecule': self._viewMolecule,
      'displaySingle': self._viewSingle,

      # Table
      'displayTable': self._viewTable,
    }


################# MAIN VIEWER FUNCTIONS ###################

  def getOutputMolSets(self):
    '''Return the sets of output molecules in the protocol'''
    outputMolSets = {}
    for oAttr in self.protocolObject.iterOutputAttributes():
      if type(getattr(self.protocolObject, oAttr[0])) == SetOfSmallMolecules:
        outputMolSets[oAttr[0]] = getattr(self.protocolObject, oAttr[0])
    return outputMolSets

  def index2MolDic(self, molSets, groupDic):
    '''From the sets of molecules selected in the output, return a dictionary with the index dictionaries
    of form {groupLabel: [mols]}
    '''
    molDic = {}
    for molSetLabel, molSet in molSets.items():
      indexDic = groupDic[molSetLabel]
      for groupLabel, indexes in indexDic.items():
        if groupLabel in molDic:
          molDic[groupLabel] += molSet.getMolsFromIds(indexes)
        else:
          molDic[groupLabel] = molSet.getMolsFromIds(indexes)
    return molDic

  def getGroupMols(self, groupDic, sLabel):
    '''Return the molecules determined by the sLabel in the dictionary of index (groupDic)
    '''
    molSets = self.getOutputMolSets()
    molDic = self.index2MolDic(molSets, groupDic)
    if sLabel != 'All':
        mols = sortMolsByUnique(molDic[sLabel])
    else:
        mols = []
        for sLab in self.setLabels:
          if sLab != 'All':
              mols += sortMolsByUnique(molDic[sLab])
    if len(mols) > 1000:
      res = askYesNoCancel('Slow display',
                  f'Trying to display a high number of molecules ({len(mols)}) might be slow, do you want to continue?',
                  self.getTkRoot())
      if res != 0:
        return []
    return mols

  def viewPymolMols(self, mols, ligandLabel, addTarget=True, disable=True, pose=True, e=None):
    pmlsDir = self.getPmlsDir()
    pmlFile = os.path.join(pmlsDir, '{}.pml'.format(ligandLabel))
    writePmlFile(pmlFile, buildPMLDockingGroupsStr(self, mols, addTarget=addTarget, disable=disable, pose=pose))

    pymolV = PyMolViewer(project=self.getProject())
    return pymolV._visualize(os.path.abspath(pmlFile), cwd=os.path.dirname(pmlFile))

  def getDockFile(self, mol):
    poseFile = os.path.abspath(mol.getPoseFile())
    molName, gridId = mol.getMolName(), mol.getGridId()

    dockDir = os.path.dirname(poseFile).replace('outputLigands', f'extra/pocket_{gridId}')
    molFile = os.path.join(dockDir, f'{molName}.pdbqt')
    if not os.path.exists(molFile):
      molFile = mol.getPoseFile()
    return molFile

  def getMolDic(self, mols, pose, vinaDock):
    molDic = {}
    for mol in mols:
      poseB, dockB = True, True
      if vinaDock:
        poseB, dockB = False, False

      molName = mol.getUniqueName(pose=poseB, dock=dockB)
      if vinaDock:
        molFile = self.getDockFile(mol)
      elif pose:
        molFile = mol.getPoseFile()
      else:
        molFile = mol.getFileName()

      molDic[molName] = molFile
    return molDic

  def writeChimeraScript(self, molDic, ligandLabel, addTarget=True, disable=True, vinaDock=False):
      pmlsDir = self.getPmlsDir()
      chimScript = os.path.join(pmlsDir, f'{ligandLabel}_chimeraX.py')
      with open(chimScript, "w") as f:
        f.write("from chimerax.core.commands import run\n")

        f.write(f"run(session, 'cd {os.getcwd()}')\n")
        f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates

        if addTarget:
          _inputStruct = os.path.abspath(self.protocol.getOriginalReceptorFile())
          f.write(f"run(session, 'open {_inputStruct}')\n")

        i = 2
        for molName, dockFile in molDic.items():
          f.write(f"run(session, 'open {dockFile} name {molName}')\n")
          if disable:
            f.write(f"run(session, 'hide #{i} models')\n")
          i += 1

        if vinaDock:
          f.write("run(session, 'viewdockx')\n")
      return chimScript

  def viewChimeraXMols(self, mols, ligandLabel, addTarget=True, disable=True,
                       pose=True, vinaDock=False, e=None):
    if not chimeraInstalled():
        print(CHIMERA_ERROR)
        return [self.warnMessage(CHIMERA_ERROR, 'Chimera not found')]
    else:
        molDic = self.getMolDic(mols, pose, vinaDock)
        chimScript = self.writeChimeraScript(molDic, ligandLabel, addTarget, disable, vinaDock)
        view = ChimeraView(chimScript)
        return [view]


  ################### DOCKING POSES VIEWS #################

  def _viewSetDock(self, e=None):
      if self.checkIfProtocol():
          ligandLabel = self.getEnumText('displaySetDock')
          sLabel = ligandLabel
      else:
          ligandLabel, sLabel = 'All', 'allSetDockedMolecules'

      mols = self.getGroupMols(self.setLigandsDic, ligandLabel)
      if len(mols) > 0:
        if self.displaySoftware.get() == PYMOL:
            return self.viewPymolMols(mols, sLabel)
        elif self.displaySoftware.get() == CHIMERAX:
            return self.viewChimeraXMols(mols, sLabel)
        elif self.displaySoftware.get() == VIEWDOCKX:
            return self.viewChimeraXMols(mols, sLabel, vinaDock=True)

  def _viewROIDock(self, e=None):
    ligandLabel = self.getEnumText('displayROIDock')

    mols = self.getGroupMols(self.pocketLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftware.get() == PYMOL:
          return self.viewPymolMols(mols, ligandLabel)
      elif self.displaySoftware.get() == CHIMERAX:
          return self.viewChimeraXMols(mols, ligandLabel)
      elif self.displaySoftware.get() == VIEWDOCKX:
        return self.viewChimeraXMols(mols, ligandLabel, vinaDock=True)


  def _viewMoleculeDock(self, e=None):
    ligandLabel = self.getEnumText('displayMoleculeDock')

    mols = self.getGroupMols(self.moleculeLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftware.get() == PYMOL:
          return self.viewPymolMols(mols, ligandLabel)
      elif self.displaySoftware.get() == CHIMERAX:
          return self.viewChimeraXMols(mols, ligandLabel)
      elif self.displaySoftware.get() == VIEWDOCKX:
        return self.viewChimeraXMols(mols, ligandLabel, vinaDock=True)

  def _viewSingleDock(self, e=None):
    ligandLabel = self.getEnumText('displaySingleDock')

    mols = self.getGroupMols(self.singleLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftware.get() == PYMOL:
          return self.viewPymolMols(mols, ligandLabel, disable=False)
      elif self.displaySoftware.get() == CHIMERAX:
          return self.viewChimeraXMols(mols, ligandLabel, disable=False)
      elif self.displaySoftware.get() == VIEWDOCKX:
        return self.viewChimeraXMols(mols, ligandLabel, vinaDock=True)

  def _viewPymolPLIP(self, e=None):
    pmlsDir = self.getPmlsDir()

    ligandLabel = self.getEnumText('displayPymolPLIP')
    mols = self.getGroupMols(self.singleLigandsDic, ligandLabel)
    if len(mols) > 0:
      mol = mols[0].clone()
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

################### MOLECULES VIEWS #################

  def _viewSet(self, e=None):
    if self.checkIfProtocol():
      ligandLabel = self.getEnumText('displaySet')
      sLabel = ligandLabel
    else:
      ligandLabel, sLabel = 'All', 'allSetMolecules'

    mols = self.getGroupMols(self.setLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftwareMol.get() == PYMOL:
        return self.viewPymolMols(mols, sLabel, pose=False, addTarget=False)
      elif self.displaySoftwareMol.get() == CHIMERAX:
        return self.viewChimeraXMols(mols, sLabel, pose=False, addTarget=False)

  def _viewMolecule(self, e=None):
    ligandLabel = self.getEnumText('displayMolecule')

    mols = self.getGroupMols(self.moleculeLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftwareMol.get() == PYMOL:
        return self.viewPymolMols(mols, ligandLabel, pose=False, addTarget=False)
      elif self.displaySoftwareMol.get() == CHIMERAX:
        return self.viewChimeraXMols(mols, ligandLabel, pose=False, addTarget=False)

  def _viewSingle(self, e=None):
    ligandLabel = self.getEnumText('displaySingle')

    mols = self.getGroupMols(self.singleLigandsDic, ligandLabel)
    if len(mols) > 0:
      if self.displaySoftwareMol.get() == PYMOL:
        return self.viewPymolMols(mols, ligandLabel, disable=False, pose=False, addTarget=False)
      elif self.displaySoftwareMol.get() == CHIMERAX:
        return self.viewChimeraXMols(mols, ligandLabel, disable=False, pose=False, addTarget=False)

    ###################   TABLE VIEW ####################

  def _viewTable(self, e=None):
    if self.checkIfProtocol():
      ligandLabel = self.getEnumText('displayTable')
      molSet = self.getOutputMolSets()[ligandLabel]
    else:
      molSet = self.protocol

    try:
      setV = MDViewer(project=self.getProject())
    except:
      setV = BioinformaticsDataViewer(project=self.getProject())
    views = setV._visualize(molSet)
    return views


#################### UTILS ###########################33

  def checkIfDocked(self):
      molSet = None
      if not self.checkIfProtocol():
          molSet = self.protocol
      else:
          for oAttr in self.protocol.iterOutputAttributes():
              if type(getattr(self.protocol, oAttr[0])) == SetOfSmallMolecules:
                  molSet = getattr(self.protocol, oAttr[0])

      if molSet is not None:
        return molSet.isDocked()
      else:
        return False

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
    cleanPDB(auxPath, outPath, waters=True, hetatm=False)
    return outPath

class CorrelationViewer(pwviewer.Viewer):
    from pwchem.protocols.VirtualDrugScreening.protocol_scores_correlation import ProtScoreCorrelation
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
      ProtScoreCorrelation,
    ]

    def __init__(self, **kwargs):
      pwviewer.Viewer.__init__(self, **kwargs)

    def _visualize(self, protocol, **kwargs):
      from PIL import Image
      image = Image.open(protocol.getOuputImgPath())
      image.show()


class SmallMoleculesLibraryViewer(pwviewer.ProtocolViewer):
    ''' Visualize the library characteristics'''
    _label = 'Viewer Small Molecules Library'
    _targets = [SmallMoleculesLibrary]
    _environments = [pwviewer.DESKTOP_TKINTER]


    def __init__(self, **args):
        super().__init__(**args)

    def _defineParams(self, form):
        form.addSection(label='Visualization of Small Molecules library')
        group = form.addGroup('Score 1')
        group.addParam('chooseHeader1', params.StringParam, default='', label='Library score 1: ',
                       help='Choose the library score to analyze')
        group.addParam('trueMin1', params.FloatParam, label='Minimum score value: ', default=0,
                       help='Minimum and maximum values for the selected score.'
                            'Use the wizard to get the actual minimum value or set yourself the histogram limit')
        group.addParam('trueMax1', params.FloatParam, label='Maximum score value: ', default=1,
                       help='Minimum and maximum values for the selected score. '
                            'Use the wizard to get the actual maximum value or set yourself the histogram limit')

        group = form.addGroup('Score 2')
        group.addParam('chooseHeader2', params.StringParam, default='', label='Library score 2: ',
                       help='Choose the library score 2 to analyze in a scatterplot')
        group.addParam('trueMin2', params.FloatParam, label='Minimum score value: ', default=0,
                       help='Minimum and maximum values for the selected score.'
                            'Use the wizard to get the actual minimum value or set yourself the histogram limit')
        group.addParam('trueMax2', params.FloatParam, label='Maximum score value: ', default=1,
                       help='Minimum and maximum values for the selected score. '
                            'Use the wizard to get the actual maximum value or set yourself the histogram limit')

        group = form.addGroup('Plots')
        group.addParam('recalc', params.BooleanParam, label='Recalculate: ', default='False',
                       help='Recalculate the graph instead of loading a previously saved')
        group.addParam('displayDistribution', params.LabelParam, label='Display histogram: ',
                       help='Display the distribution of the score 1 values as histogram')
        group.addParam('displayScatterplot', params.LabelParam, label='Display scatterplot: ',
                       help='Display the distribution of the score 1 and 2 values as scatterplot')
        group.addParam('displayHeatmap', params.LabelParam, label='Display heatmap: ',
                       help='Display the distribution of the score 1 and 2 values as heatmap')


    def _getVisualizeDict(self):
        return {
            'displayDistribution': self._showDistribution,
            'displayScatterplot': self._showScatterplot,
            'displayHeatmap': self._showHeatmap,
        }

    def getLibrary(self):
      return self.protocol

    def getMinValue(self, scoreIdx, doRound=True):
      minVal = 1e6
      inLib, head = self.getLibrary(), getattr(self, f'chooseHeader{scoreIdx}').get()
      colIdx = self.getHeaders().index(head)
      for batch in inLib.yieldLibraryValues([colIdx]):
        fBatch = [float(val[0]) for val in batch]
        minVal = min(minVal, min(fBatch))

      if doRound:
        minVal = round(minVal, 2)
      return minVal

    def getMaxValue(self, scoreIdx, doRound=True):
      maxVal = -1e6
      inLib, head = self.getLibrary(), getattr(self, f'chooseHeader{scoreIdx}').get()
      colIdx = self.getHeaders().index(head)
      for batch in inLib.yieldLibraryValues([colIdx]):
        fBatch = [float(val[0]) for val in batch]
        maxVal = max(maxVal, max(fBatch))

      if doRound:
        maxVal = round(maxVal, 2)
      return maxVal

    def getHeaders(self):
      return self.getLibrary().getHeaders()

    def _showDistribution(self, e=None):
      pngFile = self.getPngFile(plot='distribution')
      if not os.path.exists(pngFile) or self.recalc.get():
        inLib, head = self.getLibrary(), self.chooseHeader1.get()
        colIdx = self.getHeaders().index(head)

        iniMin, iniMax = self.trueMin1.get(), self.trueMax1.get()
        initialBins = 50

        binEdges = np.linspace(iniMin, iniMax, initialBins + 1)
        histCounts = np.zeros(initialBins)
        for batch in inLib.yieldLibraryValues([colIdx]):
          fBatch = [float(val[0]) for val in batch]
          histCounts += np.histogram(fBatch, bins=binEdges)[0]

          underflow = np.sum(np.array(fBatch) < iniMin)
          overflow = np.sum(np.array(fBatch) > iniMax)
          histCounts[0] += underflow
          histCounts[-1] += overflow

        plt.bar(binEdges[:-1], histCounts, width=np.diff(binEdges), align='edge', edgecolor='black')
        plt.title(f"Histogram of {head}")
        plt.xlabel(f"{head} Values")
        plt.ylabel("Frequency")
        plt.tight_layout()

        plt.savefig(pngFile)
      else:
        from PIL import Image
        image = Image.open(pngFile)
        image.show()

      plt.show()

    def _showScatterplot(self, e=None):
      pngFile = self.getPngFile(plot='scatterplot', nVars=2)
      if not os.path.exists(pngFile) or self.recalc.get():
        inLib, headers = self.getLibrary(), self.getHeaders()
        head1, head2 = self.chooseHeader1.get(), self.chooseHeader2.get()
        colIdx1, colIdx2 = headers.index(head1), headers.index(head2)

        for batch in inLib.yieldLibraryValues([colIdx1, colIdx2]):
          fBatch1, fBatch2 = [float(val[0]) for val in batch], [float(val[1]) for val in batch]
          plt.scatter(fBatch1, fBatch2, color='blue', s=2)

        plt.title(f"Scatterplot of {head1} vs {head2}")
        plt.xlabel(f"{head1} Values")
        plt.ylabel(f"{head2} Values")
        plt.xlim(self.trueMin1.get(), self.trueMax1.get())
        plt.ylim(self.trueMin2.get(), self.trueMax2.get())

        plt.legend()
        plt.grid(True)

        plt.savefig(pngFile)
      else:
        from PIL import Image
        image = Image.open(pngFile)
        image.show()

      plt.show()

    def _showHeatmap(self, e=None):
      pngFile = self.getPngFile(plot='heatmap', nVars=2)
      if not os.path.exists(pngFile) or self.recalc.get():
        xMin, xMax = self.trueMin1.get(), self.trueMax1.get()
        yMin, yMax = self.trueMin2.get(), self.trueMax2.get()
        nBinsX, nBinsY = 50, 50

        # Bin edges
        xEdges = np.linspace(xMin, xMax, nBinsX + 1)
        yEdges = np.linspace(yMin, yMax, nBinsY + 1)

        # Accumulator for 2D histogram
        hist = np.zeros((nBinsX, nBinsY))

        inLib, headers = self.getLibrary(), self.getHeaders()
        head1, head2 = self.chooseHeader1.get(), self.chooseHeader2.get()
        colIdx1, colIdx2 = headers.index(head1), headers.index(head2)
        # Incrementally update histogram
        for batch in inLib.yieldLibraryValues([colIdx1, colIdx2]):
          fBatch1, fBatch2 = [float(val[0]) for val in batch], [float(val[1]) for val in batch]
          h, _, _ = np.histogram2d(fBatch1, fBatch2, bins=[xEdges, yEdges])
          hist += h

        # Plot heatmap
        plt.figure(figsize=(8, 6))
        plt.imshow(hist.T, origin='lower', aspect='auto',
                   extent=[xMin, xMax, yMin, yMax],
                   cmap='plasma')
        plt.colorbar(label='Counts')
        plt.xlabel(f"{head1} Values")
        plt.ylabel(f"{head2} Values")
        plt.xlim(self.trueMin1.get(), self.trueMax1.get())
        plt.ylim(self.trueMin2.get(), self.trueMax2.get())
        plt.title(f"Heatplot of {head1} vs {head2}")
        plt.grid(False)
        plt.show()

    def getPngFile(self, nVars=1, plot='distribution'):
      oDir = os.path.abspath(os.path.dirname(self.getLibrary().getFileName()))
      hVars = f'{self.chooseHeader1.get()}' if nVars == 1 else f'{self.chooseHeader1.get()}_{self.chooseHeader2.get()}'
      return os.path.join(oDir, f'{hVars}_{plot}.png')
