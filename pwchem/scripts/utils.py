# Utils for the scripts
import os
from rdkit import Chem

def parseParams(paramsFile, listParams=[], sep=':'):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split(sep)
      if key in listParams:
        paramsDic[key] = value.strip().split()
      else:
        paramsDic[key] = value.strip()
  return paramsDic

def parseMoleculeFile(molFile):
  if molFile.endswith('.mol2'):
    mol = Chem.MolFromMol2File(molFile)
  elif molFile.endswith('.mol'):
    mol = Chem.MolFromMolFile(molFile)
  elif molFile.endswith('.pdb'):
    mol = Chem.MolFromPDBFile(molFile)
  elif molFile.endswith('.smi'):
    f = open(molFile, "r")
    firstline = next(f)
    mol = Chem.MolFromSmiles(str(firstline))
  elif molFile.endswith('.sdf'):
    suppl = Chem.SDMolSupplier(molFile)
    for mol in suppl:
      break
  else:
    mol = Chem.MolFromSmiles(molFile)

  return mol

def getMolFilesDic(molFiles):
  molsDict = {}
  for molFile in molFiles:
    m = parseMoleculeFile(molFile)
    molsDict[m] = molFile

  mols = list(molsDict.keys())
  return molsDict, mols

def writeMol(mol, outFile, cid=-1, setName=False):
  w = Chem.SDWriter(outFile)
  molName = os.path.split(os.path.splitext(outFile)[0])[-1]
  if setName:
    mol.SetProp('_Name', molName)
  w.write(mol, cid)
  w.close()

def getBaseName(file):
  return os.path.splitext(os.path.basename(file.strip()))[0]
