# Utils for the scripts
import os

def getBaseName(file):
  return os.path.splitext(os.path.basename(file.strip()))[0]

def parseParams(paramsFile, listParams=[], evalParams=[], sep=':'):
  paramsDic = {}
  with open(paramsFile) as f:
    for line in f:
      key, value = line.strip().split(sep)
      if key in listParams:
        paramsDic[key.strip()] = value.strip().split()
      elif key in evalParams:
        paramsDic[key] = eval(value.strip())
      else:
        paramsDic[key.strip()] = value.strip()
  return paramsDic


# RDKIT utils
def parseMoleculeFile(molFile):
  from rdkit import Chem
  if molFile.endswith('.mol2'):
    mol = Chem.MolFromMol2File(molFile)
  elif molFile.endswith('.mol'):
    mol = Chem.MolFromMolFile(molFile)
  elif molFile.endswith('.pdb'):
    mol = Chem.MolFromPDBFile(molFile)
  elif molFile.endswith('.smi'):
    with open(molFile, "r") as f:
      line = f.readline()
      if line.startswith('SMILES'):
        line = f.readline()
      mol = Chem.MolFromSmiles(line)
  elif molFile.endswith('.sdf'):
    suppl = Chem.SDMolSupplier(molFile)
    for mol in suppl:
      break
  else:
    mol = Chem.MolFromSmiles(molFile)

  return mol

def writeMol(mol, outFile, cid=-1, setName=False):
  from rdkit import Chem
  w = Chem.SDWriter(outFile)
  molName = os.path.split(os.path.splitext(outFile)[0])[-1]
  if setName:
    mol.SetProp('_Name', molName)
  w.write(mol, cid)
  w.close()

def getMolFilesDic(molFiles):
  molsDict = {}
  for molFile in molFiles:
    m = parseMoleculeFile(molFile)
    molsDict[m] = molFile

  mols = list(molsDict.keys())
  return molsDict, mols
