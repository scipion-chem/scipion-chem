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

def typeParamsDic(paramsDic, typeDic):
  '''Return a paramsDic with the values typed as determined in the typeDic.
  paramsDic: {parkey: parValue}
  typeDic: {parTypeFunc: [parKeys]}

  returns: {parkey: parTypeFunc(parValue)}
  '''
  for typeFunc, parKeys in typeDic.items():
    for parKey in parKeys:
      if parKey in paramsDic:
        paramsDic[parKey] = typeFunc(paramsDic[parKey])
  return paramsDic

# RDKIT utils
def parseMoleculeFile(molFile, sanitize=True):
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
    suppl = Chem.SDMolSupplier(molFile, sanitize=sanitize)
    for mol in suppl:
      break
  else:
    mol = Chem.MolFromSmiles(molFile)

  return mol

def writeMol(mol, outFile, cid=-1, setName=False):
  from rdkit import Chem
  try:
    w = Chem.SDWriter(outFile)
    molName = os.path.split(os.path.splitext(outFile)[0])[-1]
    if setName:
      mol.SetProp('_Name', molName)
    w.write(mol, cid)
    w.close()
    return mol
  except:
    return False

def getMolFilesDic(molFiles, sanitize=True):
  molsDict = {}
  for molFile in molFiles:
    m = parseMoleculeFile(molFile, sanitize=sanitize)
    if m:
      molsDict[m] = molFile

  mols = list(molsDict.keys())
  return molsDict, mols
