import os

from pwchem.objects import SetOfSmallMolecules
from pyworkflow.protocol.params import EnumParam, LabelParam
import pyworkflow.viewer as pwviewer
from pwchem.viewers import PyMolViewer, BioinformaticsDataViewer
from pwchem.utils.utilsViewer import *
from pwchem.utils import runOpenBabel, mergePDBs, clean_PDB
from pyworkflow.gui.dialog import showError
from pwchem import Plugin as pwchemPlugin

SINGLE, MOLECULE, POCKET, SET = 'single', 'molecule', 'pocket', 'set'

class SmilesViewer(pwviewer.ProtocolViewer):
  _label = 'Viewer SMILES'
  _targets = [SetOfSmallMolecules]
  _environments = [pwviewer.DESKTOP_TKINTER]
