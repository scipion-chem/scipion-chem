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

    mp = AllChem.MMFFGetMoleculeProperties(mol, mmffVariant=ffMethod)
    AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant=ffMethod)

    for cid in cids:
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        rdMolAlign.AlignMolConformers(mol)

        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))

    return mol

def filterMolsSize(mols, size):
    nMols = []
    for mol in mols:
        if mol.GetNumAtoms() > int(size):
            nMols.append(mol)
    return nMols

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

            try:
                if paramsDic['doHydrogens']:
                    mol = Chem.AddHs(Chem.RemoveHs(mol), addCoords=True)

                if paramsDic['doGasteiger']:
                    AllChem.ComputeGasteigerCharges(mol)

                AllChem.MMFFOptimizeMoleculeConfs(mol, numThreads=0, mmffVariant=ffMethod)

                if 'numConf' in paramsDic:
                    mol = conformer_generation(mol, outBasef, ffMethod, paramsDic['restrainMethod'],
                                                     paramsDic['numConf'], paramsDic['rmsThres'])
                    for cid, cMol in enumerate(mol.GetConformers()):
                        outFile = outBasef + '-{}.sdf'.format(cid+1)
                        writeMol(mol, outFile, cid=cid)
                else:
                    outFile = outBasef + '.sdf'
                    writeMol(mol, outFile)
            except:
                failedMols.append(outBase)

    with open(os.path.join(outDir, 'failedPreparations.txt'), 'w') as f:
        for oBase in failedMols:
            f.write(oBase + '\n')

