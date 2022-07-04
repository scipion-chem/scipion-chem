from pyworkflow.gui import ListTreeProviderString, dialog
from pyworkflow.object import String
from pyworkflow.wizard import Wizard
from pwem.wizards import GetStructureChainsWizard, EmWizard
from pwem.convert import AtomicStructHandler
import pyworkflow.object as pwobj
import pyworkflow.wizard as pwizard
import os, requests
import pwchem.protocols as emprot

class SelectFileNameWizard(pwizard.Wizard):
  """Base wizard for selecting an attribute from those contained in the items of a given input
  inputParam: Name of the input parameter where the items are stored
  outputParam:
  """
  _targets, _inputs, _outputs = [], {}, {}

  def addTarget(self, protocol, targets, inputs, outputs):
      self._targets += [(protocol, targets)]
      self._inputs.update({protocol: inputs})
      self._outputs.update({protocol: outputs})

  def getInputOutput(self, form):
      '''Retrieving input and output corresponding to the protocol where the wizard is used'''
      outParam = ''
      for target in self._targets:
          if form.protocol.__class__ == target[0]:
              inParam, outParam = self._inputs[target[0]], self._outputs[target[0]]
      return inParam, outParam

  def getInputFiles(self, form, inputParam):

      inputPointer = getattr(form.protocol, inputParam)


      #fileNames = []

      #for mol in inputPointer:
          #mol = mol.get()
          #fileNames.append(os.path.abspath(mol.getFileName()))

      return inputPointer.get()

  def show(self, form, *params):
    inputParam, outputParam = self.getInputOutput(form)
    filesList = self.getInputFiles(form, inputParam[0])
    finalFilesList = []

    for i in filesList:
      finalFilesList.append(pwobj.String(i))

    provider = ListTreeProviderString(finalFilesList)

    dlg = dialog.ListDialog(form.root, "Filter set", provider,
                                "Select one of the molecules")


    form.setVar(outputParam[0], dlg.values[0].get())




#Defining target for the SelectAttributeWizard
SelectFileNameWizard().addTarget(protocol=emprot.ProtocolFingerprintFiltering,
                                  targets=['inputReferenceMolecule'],
                                  inputs=['inputSmallMolecules'],
                                  outputs=['inputReferenceMolecule'])

SelectFileNameWizard().addTarget(protocol=emprot.ProtocolShapeDistancesFiltering,
                                  targets=['inputReferenceMolecule'],
                                  inputs=['inputSmallMolecules'],
                                  outputs=['inputReferenceMolecule'])
