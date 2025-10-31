import sys

from utils import getMolFilesDic, parseParams

######################################################
def atomTypeFilt(files, minNum, atomType):
    listPass = []
    molDic, _ = getMolFilesDic(files)
    for mol, name in molDic.items():
        if mol:
            atomCount = {}
            for atom in mol.GetAtoms():
                aSymbol = atom.GetSymbol()
                if aSymbol in atomCount:
                    atomCount[aSymbol] += 1
                else:
                    atomCount[aSymbol] = 1

            if atomType in atomCount and atomCount[atomType] >= int(minNum):
                listPass.append(name)
    return listPass

def atomNumFilt(files, minAtoms):
    listPass = []
    molDic, _ = getMolFilesDic(files)
    for mol, name in molDic.items():
        if mol:
            nAtoms = len(mol.GetAtoms())
            if nAtoms >= int(minAtoms):
                listPass.append(name)
    return set(listPass)

def cycleNumFilt(files, minCycles):
    listPass = []
    molDic, _ = getMolFilesDic(files)
    for mol, name in molDic.items():
        infoCycles = mol.GetRingInfo()
        nCycles = infoCycles.NumRings()
        if nCycles >= int(minCycles):
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
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'], evalParams=['filters'], sep='::')
    files = paramsDic['ligandFiles']
    filDic = paramsDic['filters']

    molDic, _ = getMolFilesDic(files)
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

