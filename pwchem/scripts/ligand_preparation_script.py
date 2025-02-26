from rdkit import Chem, RDConfig
from rdkit.Chem import AllChem, rdMolAlign, rdShapeHelpers, rdDistGeom
import sys, os

from utils import getMolFilesDic, parseParams, writeMol

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

###################################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'])
    molFiles = paramsDic['ligandFiles']
    outDir = paramsDic['outputDir']

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
                mol = conformer_generation(mol, outBasef, paramsDic['ffMethod'], paramsDic['restrainMethod'],
                                                 paramsDic['numConf'], paramsDic['rmsThres'])
                if mol:
                    setMolName = not mol.HasProp('_Name') or not mol.GetProp('_Name')
                    for cid, cMol in enumerate(mol.GetConformers()):
                        outFile = outBasef + '-{}.sdf'.format(cid+1)
                        writeMol(mol, outFile, cid=cid, setName=setMolName)
                else:
                    failedMols.append(outBase)
            else:
                mol = embedAndOptimize(mol)
                if mol:
                    setMolName = not mol.HasProp('_Name') or not mol.GetProp('_Name')
                    outFile = outBasef + '.sdf'
                    writeMol(mol, outFile, setName=setMolName)
                else:
                    failedMols.append(outBase)


    paramsIdx = sys.argv[1].split('_')[-1].split('.')[0]
    if len(failedMols) > 0:
        with open(os.path.join(outDir, f'failedPreparations_{paramsIdx}.txt'), 'w') as f:
            for oBase in failedMols:
                f.write(oBase.split('/')[-1] + '\n')
