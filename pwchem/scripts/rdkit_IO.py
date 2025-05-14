# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *
# * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
# *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************

'''Script to convert molecule files using the rdkit Chem module in the rdkit-env.
Mainly used for parsinf mae files (which openbabel is not able to read)'''

import sys, os, argparse, shutil, gzip, threading
from rdkit import Chem
from rdkit.Chem import AllChem

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'utils'))
from scriptUtils import *

def Mol2MolSupplier(file=None,sanitize=True):
    mols=[]
    with open(file, 'r') as f:
        line =f.readline()
        while not f.tell() == os.fstat(f.fileno()).st_size:
            if line.startswith("@<TRIPOS>MOLECULE"):
                mol = []
                mol.append(line)
                line = f.readline()
                while not line.startswith("@<TRIPOS>MOLECULE"):
                    mol.append(line)
                    line = f.readline()
                    if f.tell() == os.fstat(f.fileno()).st_size:
                        mol.append(line)
                        break
                mol[-1] = mol[-1].rstrip() # removes blank line at file end
                block = "".join(mol)
                m=Chem.MolFromMol2Block(block, sanitize=sanitize)
            if m:
                mols.append(m)
            else:
                print('Error parsing mol2 file {} with RDKit. Mol2 format is sometimes tricky to parse, use Corina '
                      'typing for RDKit'.format(file))
    return mols

def decompressFile(inFile):
    if inFile.endswith('.gz'):
        newInputFile = inputFile.replace('.gz', '')
    elif inFile.endswith('gz'):
        newInputFile = inputFile.replace('gz', '')
    else:
        print('Decompress failed for file {}'.format(inFile))

    with gzip.open(inputFile) as fIn:
        with open(newInputFile, 'wb') as f:
            shutil.copyfileobj(fIn, f)

    return newInputFile

def getMolsFromFile(inFile, ext=None, nameKey=None):
    '''Parse molecules stored in a file based on extension'''
    if not ext:
        ext = os.path.splitext(inFile)[1][1:]

    if not nameKey:
        nameKey = '_Name'

    if ext == 'maegz' or ext == 'gz':
        inFile = decompressFile(inFile)

    mols = []
    if ext == 'mae':
        mols = list(Chem.MaeMolSupplier(inFile))

    elif ext == 'smi' or ext == 'smiles':
        with open(inFile) as f:
            for line in f:
                smi, name = line.strip().split()
                mol = Chem.MolFromSmiles(smi)
                if mol:
                    mol.SetProp(nameKey, name)
                    mols.append(mol)

    elif ext == 'mol2':
        mols = Mol2MolSupplier(inFile)

    elif ext == 'sdf' or ext == 'sd':
        mols = list(Chem.SDMolSupplier(inFile))

    elif ext == 'pdb':
        for pdbBlock in divideMultiPDB(inFile):
            mols.append(Chem.MolFromPDBBlock(pdbBlock))

    else:
        print('Unrecognized format {} for file {}'.format(ext, inFile))

    return mols, nameKey

def getMolName(mol, nameKey, outBase):
    if mol.HasProp(nameKey):
        molName = mol.GetProp(nameKey)
    else:
        molName = '{}_{}'.format(outBase, i + 1)
    return molName

def make3DCoords(mols, mols3dLists, it, errBase):
    '''Optimize the 3D coordinates of a rdkit molecule'''
    for i, mol in enumerate(mols):
        mol2 = Chem.AddHs(mol, addCoords=True)
        AllChem.EmbedMolecule(mol2)
        try:
            AllChem.MMFFOptimizeMolecule(mol2)
        except:
            if mol.HasProp('_Name'):
                print('Could not optimize 3D structure of molecule: {}'.format(mol.GetProp('_Name')))
                errFile = '{}_{}.txt'.format(errBase, it)
                mode = 'a' if os.path.exists(errFile) else 'w'
                with open(errFile, mode) as f:
                    f.write(mol.GetProp('_Name') + '\n')

        if len(mol.GetAtoms()) != len(mol2.GetAtoms()):
            mol2 = Chem.RemoveHs(mol2)

        mols3dLists[it].append(mol2)
    return mols3dLists[it]

if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputFilename <mol(s)File> -of/--outputFormat <outputFormat> 
    -o/--outputName [<outputName>] [<outputDirectory>] 
    The script will parse the input molFile (which can have one or several molecules) and write the molecules
    in the specified output format.
    If an outputName is specified, the output molecules will be written in a single file with that name. Else,
    each the name of the output file will tried to be parsed from the molecule name (if not found, just numbering)
    If the output directory is not specified, molecule file(s) will be saved in the input file directory
    '''
    parser = argparse.ArgumentParser(description='Handles the IO for molecule files using openbabel')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input molecule file')
    parser.add_argument('-of', '--outputFormat', type=str, required=False, default='sdf', help='Output format')
    parser.add_argument('-o', '--outputName', type=str, required=False, help='Output name')
    parser.add_argument('-ob', '--outputBase', type=str, required=False, help='Output basename for multiple outputs')
    parser.add_argument('-od', '--outputDir', type=str, required=False, help='Output directory')
    parser.add_argument('--make3D', default=False, action='store_true', help='Optimize 3D coordinates')
    parser.add_argument('--overWrite', default=False, action='store_true', help='Overwrite output')
    parser.add_argument('--nameKey', default='', type=str, required=False, help='molecule name key in file')
    parser.add_argument('-nt', '--nthreads', default=1, type=int, required=False, help='Number of threads')

    args = parser.parse_args()
    inputFile, outFormat = args.inputFilename, args.outputFormat
    outFormat = outFormat if not outFormat.startswith('.') else outFormat[1:]
    inFormat = os.path.splitext(inputFile)[1][1:]

    if args.outputName:
        singleOutFile, outName = True, os.path.splitext(args.outputName)[0]
    else:
        singleOutFile, outName = False, None

    if args.outputDir:
        outDir = args.outputDir
    else:
        outDir = os.path.dirname(inputFile)

    overW = args.overWrite
    make3d = args.make3D
    nameKey = args.nameKey
    nt = args.nthreads

    mols, nameKey = getMolsFromFile(inputFile, nameKey=nameKey)
    if len(mols) > 0:
        if make3d:
            mols = performBatchThreading(make3DCoords, mols, nt, cloneItem=False,
                                         errBase=os.path.join(outDir, 'errors3D'))

        if outFormat == 'smi' or outFormat == 'smiles':
            writter, ext = Chem.SmilesWriter, 'smi'
        elif outFormat == 'pdb':
            writter, ext = Chem.PDBWriter, 'pdb'
        else:
            writter, ext = Chem.SDWriter, 'sdf'

        if singleOutFile:
            outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(outName, ext)))
            with writter(outFile) as f:
                for mol in mols:
                    if mol:
                        f.write(mol)

        else:
            outBase = args.outputBase if args.outputBase else 'molecule'
            for i, mol in enumerate(mols):
                if mol:
                    if mol.HasProp(nameKey):
                        molName = mol.GetProp(nameKey)
                    else:
                        molName = '{}_{}'.format(outBase, i+1)
                    molName = molName.replace('/', '-')
                    outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(molName, ext)))
                    with writter(outFile) as f:
                        f.write(mol)



