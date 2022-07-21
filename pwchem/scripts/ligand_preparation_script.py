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

    res, confRes = [], {}
    for cid in cids:

        outFile = outBase + '-{}.sdf'.format(cid)
        ff = AllChem.MMFFGetMoleculeForceField(mol, mp, confId=cid)
        e = ff.CalcEnergy()
        res.append((cid, e))
        confRes[outFile] = e
        sorted_res = sorted(res, key=lambda x: x[1])
        rdMolAlign.AlignMolConformers(mol)

        w = Chem.SDWriter(outFile)
        mol.SetProp('CID', str(cid))
        mol.SetProp('Energy', str(e))
        w.write(mol, confId=cid)
        w.close()
    return confRes


###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']

    mols_dict, _ = getMolFilesDic(molFiles)
    for mol in mols_dict:
        molBase = os.path.splitext(os.path.basename(mols_dict[mol]))[0]
        if paramsDic['doHydrogens']:
            mol = Chem.AddHs(Chem.RemoveHs(mol))


        if paramsDic['doGasteiger']:
            AllChem.ComputeGasteigerCharges(mol)

        outBase = os.path.join(outDir, molBase)
        confFiles = conformer_generation(mol, outBase, paramsDic['ffMethod'], paramsDic['restrainMethod'],
                                         paramsDic['numConf'], paramsDic['rmsThres'])

