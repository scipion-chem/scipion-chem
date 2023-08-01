from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys, os

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
    mols_dict = {}
    for molFile in molFiles:
        m = parseMoleculeFile(molFile)
        mols_dict[m] = molFile

    mols = list(mols_dict.keys())
    return mols_dict, mols


def parseParams(paramsFile):
    paramsDic = {}
    with open(paramsFile) as f:
        for line in f:
            key, value = line.strip().split(':')
            if key == 'ligandFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic

#################################################################################################################
def conformer_generation(mol, outBase, ffMethod, restrainMethod, numConf, rmsThres):
    param = getattr(rdDistGeom, restrainMethod)()
    param.pruneRmsThresh = float(rmsThres)
    cids = rdDistGeom.EmbedMultipleConfs(mol, int(numConf), param)
    if len(cids) > 0:
        AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant=ffMethod)
        rdMolAlign.AlignMolConformers(mol)
        return mol
    else:
        return False

def filterMolsSize(mols, size):
    nMols = []
    for mol in mols:
        if mol.GetNumAtoms() > int(size):
            nMols.append(mol)
    return nMols

def embedAndOptimize(mol):
    embedMess = AllChem.EmbedMolecule(mol, randomSeed=44)
    if embedMess == 0:
        AllChem.MMFFOptimizeMolecule(mol)
    else:
        embedMess = AllChem.EmbedMolecule(mol, randomSeed=44, useRandomCoords=True)
        if embedMess == 0:
            AllChem.MMFFOptimizeMolecule(mol)
        else:
            return False
    return mol

def writeMol(mol, outFile, cid=-1):
    w = Chem.SDWriter(outFile)
    w.write(mol, cid)
    w.close()

###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']
    ffMethod = paramsDic['ffMethod']

    failedMols = []
    mols_dict, _ = getMolFilesDic(molFiles)
    for inmol in mols_dict:
        molBase = os.path.splitext(os.path.basename(mols_dict[inmol]))[0]
        outBase = os.path.join(outDir, molBase)
        if 'numAtoms' in paramsDic:
            inMols = []
            frags = Chem.GetMolFrags(inmol, asMols=True)
            frags = filterMolsSize(frags, paramsDic['numAtoms'])
        else:
            frags = [inmol]

        for i, mol in enumerate(frags):
            if len(frags) > 1:
                outBasef = outBase + 'f{}'.format(i)
            else:
                outBasef = outBase

            if paramsDic['doHydrogens']:
                mol = Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)

            if paramsDic['doGasteiger']:
                AllChem.ComputeGasteigerCharges(mol)

            if 'numConf' in paramsDic:
                mol = conformer_generation(mol, outBasef, ffMethod, paramsDic['restrainMethod'],
                                                 paramsDic['numConf'], paramsDic['rmsThres'])
                if mol:
                    for cid, cMol in enumerate(mol.GetConformers()):
                        outFile = outBasef + '-{}.sdf'.format(cid+1)
                        writeMol(mol, outFile, cid=cid)
                else:
                    failedMols.append(outBase)
            else:
                mol = embedAndOptimize(mol)
                if mol:
                    outFile = outBasef + '.sdf'
                    writeMol(mol, outFile)
                else:
                    failedMols.append(outBase)


    paramsIdx = sys.argv[1].split('_')[-1].split('.')[0]
    if len(failedMols) > 0:
        with open(os.path.join(outDir, f'failedPreparations_{paramsIdx}.txt'), 'w') as f:
            for oBase in failedMols:
                f.write(oBase.split('/')[-1] + '\n')
