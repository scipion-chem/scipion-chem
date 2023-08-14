import sys
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

###Funcion para procesar los archivos de input

def preprocessLigands(ligandsFiles):
    molDic = {}
    for molFile in ligandsFiles:
        if molFile.endswith('.mol2'):
            m = Chem.MolFromMol2File(molFile)
            molDic[m] = molFile

        elif molFile.endswith('.mol'):
            m = Chem.MolFromMolFile(molFile)
            molDic[m] = molFile

        elif molFile.endswith('.pdb'):
            m = Chem.MolFromPDBFile(molFile)
            molDic[m] = molFile

        elif molFile.endswith('.smi'):
            f = open(molFile, "r")
            firstline = next(f)
            m = Chem.MolFromSmiles(str(firstline))
            molDic[m] = molFile

        elif molFile.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(molFile)
            for mol in suppl:
                molDic[mol] = molFile

    return molDic


def parseParams(paramsFile):
    paramsDic = {}
    with open(paramsFile) as f:
        for line in f:
            key, value = line.strip().split('::')
            if key == 'ligandFiles':
                paramsDic[key] = value.strip().split()
            elif key == 'filters':
                paramsDic[key] = eval(value.strip())
            else:
                paramsDic[key] = value.strip()
    return paramsDic


######################################################
def atomTypeFilt(files, minNum, atomType):
    listPass = []
    molDic = preprocessLigands(files)
    for mol, name in molDic.items():
        atomCount = {}
        for atom in mol.GetAtoms():
            aSymbol = atom.GetSymbol()
            if aSymbol in atomCount:
                atomCount[aSymbol] += 1
            else:
                atomCount[aSymbol] = 1

        if atomType in atomCount and atomCount[atomType] >= minNum:
            listPass.append(name)
    return listPass

def atomNumFilt(files, minAtoms):
    listPass = []
    molDic = preprocessLigands(files)
    for mol, name in molDic.items():
        nAtoms = len(mol.GetAtoms())
        if nAtoms >= minAtoms:
            listPass.append(name)
    return set(listPass)

def cycleNumFilt(files, minCycles):
    listPass = []
    molDic = preprocessLigands(files)
    for mol, name in molDic.items():
        infoCycles = mol.GetRingInfo()
        nCycles = infoCycles.NumRings()
        if nCycles >= minCycles:
            listPass.append(name)
    return set(listPass)

def performFilter(files, fType, filter):
    if fType == 'typeAtom':
        listPass = atomTypeFilt(files, minNum=filter[1], atomType=filter[2])
    elif fType == 'numAtoms':
        listPass = atomNumFilt(files, minAtoms=filter[1])
    elif fType == 'numCycles':
        listPass = cycleNumFilt(files, minCycles=filter[1])
    return set(listPass)


#################################################################################################################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    files = paramsDic['ligandFiles']
    filDic = paramsDic['filters']

    molDic = preprocessLigands(files)
    listPass = set(molDic.values())

    for fType in filDic:
        for filter in filDic[fType]:
            newListPass = performFilter(files, fType, filter)
            if filter[0] == 'Remove':
                listPass = listPass.difference(newListPass)
            else:
                listPass = listPass.intersection(newListPass)

    with open(paramsDic['outputPath'], 'w') as f:
        f.write("#The following molecules have been passed the defined filters:" + "\n")
        for molName in listPass:
            f.write(molName + "\n")

#################################################################################################################

