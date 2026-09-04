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
import json
import pandas as pd
import matplotlib.pyplot as plt
from tkinter.messagebox import askokcancel

import os
from subprocess import Popen

import pyworkflow.viewer as pwviewer
from pyworkflow.protocol import params
import pyworkflow.utils as pwutils

import pwem.viewers.views as pwemViews
from pwem.viewers import Chimera, ChimeraView
import pwem.viewers.showj as showj
from pwem.viewers.mdviewer.viewer import MDViewer
from pwem.protocols import EMProtocol
from pwem.objects import SetOfSequences, AtomStruct, SetOfAtomStructs

from Bio.PDB import PDBParser, MMCIFParser, PPBuilder

import pwchem.objects
from pwchem import Plugin as pwchemPlugin
from pwchem.constants import *
from pwchem.utils import getBaseName


from pwchem.viewers.viewer_interaction import BaseInteractionViewer
from pwem.protocols import ProtSubSet

class PyMol:
  """ Help class to run PyMol and manage its environment. """

  @classmethod
  def getEnviron(cls):
    """ Return the proper environ to launch PyMol.
    PyMol_HOME variable is read from the ~/.config/scipion.conf file.
    """
    environ = pwutils.Environ(os.environ)
    environ.set('PATH', pwchemPlugin.getProgramHome(OPENBABEL_DIC, path='bin'),
                position=pwutils.Environ.BEGIN)
    return environ


class PyMolView(pwviewer.CommandView):
  """ View for calling an external command. """

  def __init__(self, pymolArgs, cwd, **kwargs):
    pwviewer.CommandView.__init__(self, [self.getPymolBin(), *pymolArgs.split()],
                                  cwd=cwd, **kwargs)
    
  def getPymolBin(self):
    return pwchemPlugin.getEnvPath(OPENBABEL_DIC, 'bin/pymol')

  def show(self):
    Popen(self._cmd, cwd=self._cwd, env=PyMol.getEnviron())


class PyMolViewer(pwviewer.Viewer):
  """ Wrapper to visualize pml objects with PyMol viewer. """
  _environments = [pwviewer.DESKTOP_TKINTER]

  def __init__(self, **args):
    pwviewer.Viewer.__init__(self, **args)

  def _visualize(self, pymolFile, cwd=None, **args):
    view = PyMolView(pymolFile, cwd)
    return [view]

class VmdViewPopen(pwviewer.CommandView):
  def __init__(self, vmdArgs, **kwargs):
    pwviewer.CommandView.__init__(self, 'vmd ' + vmdArgs, **kwargs)

  def show(self):
      fullProgram = '%s && %s' % (pwchemPlugin.getEnvActivationCommand(VMD_DIC), self._cmd)
      Popen(fullProgram, cwd=self._cwd, shell=True)

import numpy as np
def plotLocalizationHistogramFromDataFrame(df):
    localizationCols = [
        'Cytoplasm',
        'Nucleus',
        'Extracellular',
        'Cell membrane',
        'Mitochondrion',
        'Plastid',
        'Endoplasmic reticulum',
        'Lysosome/Vacuole',
        'Golgi apparatus',
        'Peroxisome',
        'Peripheral',
        'Transmembrane',
        'Lipid anchor',
        'Soluble'
    ]

    # Keep only columns that are actually present
    localizationCols = [
        col for col in localizationCols
        if col in df.columns
    ]

    proteinIds = df['Protein_ID'].astype(str).tolist()

    nProteins = len(proteinIds)
    nLocalizations = len(localizationCols)

    x = np.arange(nLocalizations)

    # Width of each individual bar
    width = 0.8 / nProteins

    fig, ax = plt.subplots(figsize=(14, 7))

    for i, proteinId in enumerate(proteinIds):
        values = df.loc[
            df['Protein_ID'].astype(str) == proteinId,
            localizationCols
        ].iloc[0].astype(float).values

        offset = (i - (nProteins - 1) / 2) * width

        ax.bar(
            x + offset,
            values,
            width,
            label=proteinId
        )

    ax.set_xlabel('Localization')
    ax.set_ylabel('Probability')
    ax.set_title('DeepLoc localization probabilities')

    ax.set_xticks(x)
    ax.set_xticklabels(
        localizationCols,
        rotation=45,
        ha='right'
    )

    ax.set_ylim(0, 1)
    ax.legend(title='Protein ID')

    plt.tight_layout()
    plt.show()

def plotLocalizationHistogram(csvFile):
    df = pd.read_csv(csvFile)
    plotLocalizationHistogramFromDataFrame(df)


class AtomStructViewer(pwviewer.ProtocolViewer):
    _label = 'Viewer AtomStruct'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [AtomStruct]
    _viewerOptions = ['PyMol', 'ChimeraX']

    def _defineParams(self, form):
      form.addSection(label='Visualization of AtomStruct')
      group = form.addGroup('AtomStruct General Viewer')
      group.addParam('displaySoftware', params.EnumParam,
                     choices=self._viewerOptions, default=0,
                     label='Display AtomStruct with: ',
                     help='Display the AtomStruct object with which software.\nAvailable: PyMol, ChimeraX')
      obj = self.getAtomStruct()

      if hasattr(obj, '_localizationPerc'):
          form.addSection(label='DeepLoc localization')
          form.addParam(
              'viewLocalization',
              params.LabelParam,
              label='Display localization probabilities: ',
              help='Display the DeepLoc predicted localization probabilities.'
          )
          form.addParam(
              'viewSequence',
              params.LabelParam,
              label='Display residue importance in Chimera: ',
              help='Display the DeepLoc predicted residue importance.'
          )

    def _getVisualizeDict(self):
        visDic = {
            'displaySoftware': self._viewAtomStruct,
        }

        obj = self.getAtomStruct()

        if hasattr(obj, '_localizationPerc'):
            visDic['viewLocalization'] = self._showLocalization
            visDic['viewSequence'] = self._showResidueImportance

        return visDic

    def _showLocalization(self, e=None):
        obj = self.getAtomStruct()

        localizationPerc = getattr(obj, '_localizationPerc', None)

        if localizationPerc is None:
            return

        if hasattr(localizationPerc, 'get'):
            localizationPerc = localizationPerc.get()

        if not os.path.exists(localizationPerc):
            raise FileNotFoundError(
                f"Localization probability file not found: {localizationPerc}"
            )

        plotLocalizationHistogram(localizationPerc)

    def _showResidueImportance(self, e=None):
        obj = self.getAtomStruct()

        if not hasattr(obj, '_residuePredictions'):
            raise ValueError(
                "No DeepLoc residue predictions available "
                "for this structure."
            )

        cxcFile = self.buildDeepLocCXCFile(obj)

        return [
            ChimeraView(cxcFile)
        ]

    def _viewAtomStruct(self, e=None):
      if self.displaySoftware.get() == 0:
        pymolViewer = AtomStructPymolViewer(project=self.getProject())
        return pymolViewer._visualize(self.getAtomStruct())
      elif self.displaySoftware.get() == 1:
        return self._viewChimera(self.getAtomStruct())

    def getAtomStruct(self):
      obj = self.protocol
      # If the input is a protocol (Analyze results was used), extract the AtomStruct obj
      if issubclass(type(obj), EMProtocol):
        for output in self.protocol.iterOutputAttributes(outputClass=AtomStruct):
          obj = output[1]
      return obj

    def _viewChimera(self, obj):
      fnCmd = os.path.abspath(self._getPath("chimera_output.cxc"))
      with open(fnCmd, 'w') as f:
        f.write('cd %s\n' % os.getcwd())
        f.write("cofr 0,0,0\n")  # set center of coordinates
        f.write("style stick\n")

        _inputVol = obj.getVolume()
        if _inputVol is not None:
          volID = 1
          dim, sampling = _inputVol.getDim()[0], _inputVol.getSamplingRate()

          f.write("open %s\n" % _inputVol.getFileName())
          x, y, z = _inputVol.getOrigin(force=True).getShifts()
          f.write("volume #%d style surface voxelSize %f\nvolume #%d origin "
                  "%0.2f,%0.2f,%0.2f\n"
                  % (volID, sampling, volID, x, y, z))
        else:
          dim, sampling = 150., 1.

        bildFileName = self._getPath("axis_output.bild")
        Chimera.createCoordinateAxisFile(dim, bildFileName=bildFileName, sampling=sampling)
        f.write("open %s\n" % bildFileName)

        f.write("open %s\n" % obj.getFileName())

        view = ChimeraView(fnCmd)
        return [view]


    def _getDeepLocResidueData(self, obj):
        csvFile = obj._residuePredictions.get()

        if not csvFile or not os.path.exists(csvFile):
            raise FileNotFoundError(
                f"DeepLoc residue prediction file not found: {csvFile}"
            )

        df = pd.read_csv(csvFile, sep=None, engine='python')
        if 'AA' not in df.columns:
            raise ValueError(
                f"Invalid DeepLoc residue file: {csvFile}"
            )
        scoreColumns = [
            col for col in df.columns
            if col != 'AA'
        ]

        if not scoreColumns:
            raise ValueError(
                f"No residue prediction columns found in {csvFile}"
            )

        scoreColumn = scoreColumns[0]
        fileName = obj.getFileName()

        if fileName.lower().endswith(('.cif', '.mmcif')):
            parser = MMCIFParser(QUIET=True)
        else:
            parser = PDBParser(QUIET=True)

        structure = parser.get_structure("protein", fileName)
        ppb = PPBuilder()
        residues = []
        structureSequence = []

        for peptide in ppb.build_peptides(structure):
            sequence = str(peptide.get_sequence())

            for residue, aa in zip(peptide, sequence):
                residues.append(residue)
                structureSequence.append(aa)

        deepLocSequence = ''.join(
            df['AA'].astype(str).tolist()
        )

        structureSequence = ''.join(structureSequence)

        if len(deepLocSequence) != len(structureSequence):
            raise ValueError(
                f"DeepLoc/structure residue count mismatch for "
                f"{obj.getFileName()}: "
                f"DeepLoc={len(deepLocSequence)}, "
                f"structure={len(structureSequence)}"
            )

        if deepLocSequence != structureSequence:
            raise ValueError(
                f"DeepLoc sequence does not match the structure "
                f"sequence for {obj.getFileName()}"
            )

        residueData = []

        for residue, score in zip(
                residues,
                df[scoreColumn]
        ):
            if pd.isna(score):
                continue

            chain = residue.get_parent()
            model = chain.get_parent()

            residueData.append({
                'model': model.id,
                'chain': chain.id,
                'residueId': residue.id,
                'score': float(score)
            })

        return scoreColumn, residueData

    def _writeDeepLocAttributeFile(self, obj, scoreColumn, residueData):
        attrFile = self._getPath(
            f"deeploc_{scoreColumn}.defattr"
        )

        with open(attrFile, 'w') as f:
            f.write("attribute: deeploc_score\n")
            f.write("recipient: residues\n\n")

            for data in residueData:
                chain = data['chain']
                _, resnum, insertionCode = data['residueId']

                if insertionCode and insertionCode != ' ':
                    residueSpec = f"/{chain}:{resnum}{insertionCode}"
                else:
                    residueSpec = f"/{chain}:{resnum}"

                f.write(
                    f"\t{residueSpec}\t{data['score']}\n"
                )

        return attrFile

    def buildDeepLocCXCFile(self, obj):
        scoreColumn, residueData = self._getDeepLocResidueData(obj)

        if not residueData:
            raise ValueError(
                "No valid DeepLoc residue predictions were found."
            )

        attrFile = self._writeDeepLocAttributeFile(
            obj,
            scoreColumn,
            residueData
        )

        cxcFile = self._getPath(
            f"deeploc_{scoreColumn}.cxc"
        )

        structureFile = os.path.abspath(
            obj.getFileName()
        )

        scores = [
            data['score']
            for data in residueData
        ]

        minScore = min(scores)
        maxScore = max(scores)

        with open(cxcFile, 'w') as f:
            f.write(
                f"open {structureFile}\n"
            )
            f.write(
                f"defattr {os.path.abspath(attrFile)}\n"
            )
            f.write(
                "cartoon\n"
            )
            f.write(
                "color byattribute r:deeploc_score "
                f"palette white:yellow:red "
                f"range {minScore},{maxScore} "
                "target r\n"
            )
            f.write(
                "lighting default\n"
            )

            f.write(
                "view\n"
            )

        return cxcFile


class AtomStructPymolViewer(PyMolViewer):
    _label = 'Pymol viewer AtomStruct'
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = []

    def _visualize(self, obj, **args):
      pymolV = PyMolViewer(project=self.getProject())
      return pymolV._visualize(obj.getFileName())



class SetOfAtomStructViewer(AtomStructViewer, BaseInteractionViewer):
  _label = 'Viewer Set Of AtomStruct'
  _targets = [SetOfAtomStructs]

  def __init__(self, **kwargs):
    pwviewer.ProtocolViewer.__init__(self, **kwargs)
    self.protocolObject = self.protocol

  def _defineParams(self, form):
    form.addSection(label='Structures view')
    form.addParam('displaySoftware', params.EnumParam, choices=self._viewerOptions, default=0,
                   label='Display atom structures with: ',
                   help='Display the atom structures with this software')
    form.addParam('alignStructures', params.BooleanParam, default=True,
                  label='Align the structures: ',
                  help='Whether to align the structures in the viewer')

    form.addSection(label='Table view')
    form.addParam('displayTable', params.LabelParam,
                  label='Display Atom Struct set in table format: ',
                  help='Display the Atom Struct set in the set in table format with their respective attributes')
    structs = self.getStructSet()
    if hasattr(structs, '_getData'):
        data = structs._getData()
        if data:
            BaseInteractionViewer._defineInteractionParams(self, form=form, data=data)

    hasLocalization = False
    if hasattr(structs, '_localizationPerc'):
        hasLocalization = True

    if hasLocalization:
        form.addSection(label='DeepLoc localization')
        form.addParam(
            'viewLocalization',
            params.LabelParam,
            label='Display localization probabilities: ',
            help='Display the DeepLoc predicted localization probabilities.'
        )
        form.addParam(
            'viewSequence',
            params.LabelParam,
            label='Display residue importance in Chimera: ',
            help='Display the DeepLoc predicted residue importance.'
        )

  def _getVisualizeDict(self):
      d = {
          'displaySoftware': self._viewSetStructure,
          'displayTable': self._viewTable,
      }

      structs = self.getStructSet()

      if hasattr(structs, '_localizationPerc'):
          d['viewLocalization'] = self._showLocalization
          d['viewSequence'] = self._showResidueImportance

      d.update(BaseInteractionViewer._getVisualizeDict(self))

      return d

  def _showLocalization(self, e=None):
      structs = self.getStructSet()

      localizationPerc = getattr(
          structs,
          '_localizationPerc',
          None
      )

      if hasattr(localizationPerc, 'get'):
          localizationPerc = localizationPerc.get()

      if not localizationPerc:
          raise FileNotFoundError(
              "No DeepLoc localization probability file was found."
          )

      if not os.path.exists(localizationPerc):
          raise FileNotFoundError(
              f"DeepLoc localization file not found: "
              f"{localizationPerc}"
          )

      df = pd.read_csv(localizationPerc)

      plotLocalizationHistogramFromDataFrame(df)

  def _showResidueImportance(self, e=None):
      structs = self.getAtomStructs()

      cxcFile = self.buildDeepLocSetCXCFile(structs)

      return [
          ChimeraView(cxcFile)
      ]

  def _viewSetStructure(self, e=None):
    if self.displaySoftware.get() == 0:
      pymolViewer = PyMolViewer(project=self.getProject())
      pmlFile = self.buildPMLFile(align=self.alignStructures.get())
      return pymolViewer._visualize(pmlFile)
    elif self.displaySoftware.get() == 1:
      view = ChimeraView(self.buildCXCFile(align=self.alignStructures.get()))
      return [view]

  def _viewTable(self, e=None):
    ASs = self.getAtomStructs()

    setV = MDViewer(project=self.getProject())
    views = setV._visualize(ASs)
    return views

  def getAtomStructs(self):
      obj = self.protocol
      if issubclass(type(obj), EMProtocol):
        for output in self.protocol.iterOutputAttributes(outputClass=SetOfAtomStructs):
          obj = output[1]
      return obj

  def buildPMLFile(self, align=False):
      ASs = self.getAtomStructs()
      oDir = os.path.dirname(ASs.getFileName())
      oFile = os.path.join(oDir, 'visualization.pml')

      fAS = ASs.getFirstItem()
      fFile, fName = fAS.getFileName(), getBaseName(fAS.getFileName())
      oStr = f'load {fFile}, {fName}\n\n'
      for struct in ASs:
          sName = getBaseName(struct.getFileName())
          oStr += f'load {struct.getFileName()}, {sName}\n'
          if align:
            oStr += f'align {sName}, {fName}\n'
          oStr += '\n'
      oStr += 'zoom'

      with open(oFile, 'w') as f:
        f.write(oStr)
      return oFile

  def buildCXCFile(self, align=False):
    ASs = self.getAtomStructs()
    oDir = os.path.dirname(ASs.getFileName())
    oFile = os.path.join(oDir, 'visualization.cxc')

    oStr = ''
    for struct in ASs:
      oStr += f'open {os.path.abspath(struct.getFileName())}\n\n'

    if align:
      for i in range(len(ASs)-1):
        oStr += f'matchmaker #{i+2} to #1\n'
    oStr += 'view'

    with open(oFile, 'w') as f:
      f.write(oStr)
    return oFile

  def buildDeepLocSetCXCFile(self, structSet):
      cxcFile = self._getPath(
          "deeploc_residue_importance.cxc"
      )

      with open(cxcFile, 'w') as f:

          for struct in structSet:

              if not hasattr(struct, '_residuePredictions'):
                  continue

              try:
                  scoreColumn, residueData = (
                      self._getDeepLocResidueData(struct)
                  )
              except Exception as exc:
                  self.warning(
                      f"Could not process DeepLoc residue "
                      f"predictions for "
                      f"{struct.getFileName()}: {exc}"
                  )
                  continue

              if not residueData:
                  continue

              attrFile = self._writeDeepLocAttributeFile(
                  struct,
                  scoreColumn,
                  residueData
              )

              structureFile = os.path.abspath(
                  struct.getFileName()
              )

              scores = [
                  data['score']
                  for data in residueData
              ]

              minScore = min(scores)
              maxScore = max(scores)

              f.write(
                  f"open {structureFile}\n"
              )
              f.write(
                  f"defattr {os.path.abspath(attrFile)}\n"
              )

              f.write("cartoon\n")

              f.write(
                  "color byattribute r:deeploc_score "
                  f"palette white:yellow:red "
                  f"range {minScore},{maxScore} "
                  "target r\n"
              )

          f.write("view\n")

      return cxcFile

  def _getEntityNames(self, data):
      protNames = sorted(data.keys())
      molNames = sorted(next(iter(data.values())).keys())
      scoreTypes = sorted(next(iter(next(iter(data.values())).values())).keys())
      return protNames, molNames, scoreTypes

  def _getLabels(self):
      return "Protein", "Ligand", "Affinity"

  def getStructSet(self):
      if hasattr(self.protocol, 'outputAtomStructs'):
          return self.protocol.outputAtomStructs
      return self.protocol

  def _generateProts(self, paramName=None):
      f1 = self.getEnumText('chooseEnt1')
      f2 = self.getEnumText('chooseEnt2')
      fScore = self.getEnumText('chooseScore')

      data = self._getData()

      _, e1, _, _ = self._getFilteredData(data, f1, f2, fScore)

      objIds = []
      structSet = self.getStructSet()

      for obj in structSet:
          if os.path.splitext(os.path.basename(obj.getFileName()))[0] in e1:
              objIds.append(str(obj.getObjId()))

      if not objIds:
          return

      if askokcancel("Generate proteins subset",
                     f"Generate subset with {len(objIds)} proteins?"):
          project = self.getProject()
          prot = project.newProtocol(
              ProtSubSet,
              inputFullSet=structSet,
              selectIds=True,
              range=','.join(objIds)
          )

          project.launchProtocol(prot, wait=True)

  def _getMolSet(self):
      structs = self.getStructSet()
      return structs.getInteractMols()

  def getMolecules(self, jsonPath):
      with open(jsonPath, "r") as f:
          data = json.load(f)

      molecules = set()

      for modelData in data.values():
          molecules.update(modelData.keys())

      return sorted(molecules)

  def getInteractionSet(self):
      return self.getStructSet()


class SetOfDatabaseIDView(pwemViews.ObjectView):
    """ Customized ObjectView for SetOfDatabaseID. """
    def __init__(self, project, inputid, path, other='',
                 viewParams={}, **kwargs):
        defaultViewParams = {showj.MODE: 'metadata',
                             showj.RENDER: '_PDBLigandImage'}
        defaultViewParams.update(viewParams)
        pwemViews.ObjectView.__init__(self, project, inputid, path, other,
                                  defaultViewParams, **kwargs)

class BioinformaticsDataViewer(pwviewer.Viewer):
    """ Wrapper to visualize different type of objects
    with the Xmipp program xmipp_showj
    """
    _environments = [pwviewer.DESKTOP_TKINTER]
    _targets = [
        SetOfSequences,
        pwchem.objects.SetOfDatabaseID,
        pwchem.objects.SetOfSmallMolecules,
        pwchem.objects.SetOfBindingSites,
    ]

    def __init__(self, **kwargs):
        pwviewer.Viewer.__init__(self, **kwargs)
        self._views = []

    def _getObjView(self, obj, fn, viewParams={}):
        return pwemViews.ObjectView(
            self._project, obj.strId(), fn, viewParams=viewParams)

    def _visualize(self, obj, **kwargs):
        views = []
        cls = type(obj)
        objectViews = (SetOfSequences, pwchem.objects.SetOfSmallMolecules, 
                       pwchem.objects.SetOfStructROIs, pwchem.objects.SetOfSequenceROIs)

        if issubclass(cls, pwchem.objects.SetOfDatabaseID):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, pwchem.objects.SetOfBindingSites):
            views.append(SetOfDatabaseIDView(self._project, obj.strId(), obj.getFileName()))
        elif issubclass(cls, objectViews):
            views.append(pwemViews.ObjectView(self._project, obj.strId(), obj.getFileName()))

        return views
