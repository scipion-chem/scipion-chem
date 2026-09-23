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

import sys, os, argparse, shutil, gzip, threading, fnmatch, tempfile
from rdkit import Chem
from rdkit.Chem import AllChem
import csv

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
    '''Decompress inFile into a temporary directory and return the path of the decompressed copy.'''
    baseName = os.path.basename(inFile)
    if baseName.endswith('.gz'):
        baseName = baseName[:-len('.gz')]
    elif baseName.endswith('gz'):
        baseName = baseName[:-len('gz')]
    else:
        print('Decompress failed for file {}'.format(inFile))
        return inFile

    newInputFile = os.path.join(tempfile.mkdtemp(prefix='rdkitIO_'), baseName)
    with gzip.open(inFile) as fIn:
        with open(newInputFile, 'wb') as f:
            shutil.copyfileobj(fIn, f)

    return newInputFile

def readMaeMols(inFile, keepHs=False):
    '''Parse a Maestro file. With keepHs the structure is preserved as it is in the file: (RDKit
    removes every hydrogen by defaul). Sanitization is tried first and dropped only if it fails'''
    if not keepHs:
        return list(Chem.MaeMolSupplier(inFile))

    try:
        return list(Chem.MaeMolSupplier(inFile, removeHs=False))
    except Exception as e:
        print('Sanitization failed for {} ({}), parsing it without sanitizing'.format(inFile, e))
        return list(Chem.MaeMolSupplier(inFile, sanitize=False, removeHs=False))

def readSmilesMols(inFile, nameKey):
    '''Parse a smiles file, taking the name from the second column when it is there'''
    mols = []
    with open(inFile) as f:
        for line in f:
            if not line.strip():
                continue

            values = line.strip().split('\t')
            if len(values) >= 2:
                smi, name = values[0], values[1]
            else:
                smi, name = values[0], os.path.basename(os.path.splitext(inFile)[0])

            mol = Chem.MolFromSmiles(smi)
            if mol:
                mol.SetProp(nameKey, name)
                mols.append(mol)
    return mols

def getMolsFromFile(inFile, ext=None, nameKey=None, keepHs=False):
    '''Parse molecules stored in a file based on extension'''
    if not ext:
        ext = os.path.splitext(inFile)[1][1:]

    if not nameKey:
        nameKey = '_Name'

    if ext == 'maegz' or ext == 'gz':
        inFile = decompressFile(inFile)
        ext = os.path.splitext(inFile)[1][1:]

    mols = []
    if ext == 'mae':
        mols = list(Chem.MaeMolSupplier(inFile))

    elif ext == 'smi' or ext == 'smiles':
        with open(inFile) as f:
            for line in f:
                if not line.strip():
                    continue

                values = line.split()
                if len(values) >= 2:
                    smi, name = values[0], values[1]
                else:
                    smi, name = values[0], os.path.basename(os.path.splitext(inFile)[0])
                
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
        return [], nameKey

    return readers[ext](), nameKey

def getMolName(mol, nameKey, outBase, idx=0):
    '''Name of the molecule, falling back to outBase plus its position (idx was undefined before)'''
    if mol.HasProp(nameKey):
        return mol.GetProp(nameKey)

    return '{}_{}'.format(outBase, idx + 1)

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

def loadInputFiles(inputFile, nameKey, keepHs=False):
    ext = os.path.splitext(inputFile)[1].lower()

    if ext == ".txt":
        mols = []
        with open(inputFile) as f:
            files = [line.strip() for line in f if line.strip()]

        for fpath in files:
            m, _ = getMolsFromFile(fpath, nameKey=nameKey, keepHs=keepHs)
            mols.extend(m)

        return mols, "_Name"
    return getMolsFromFile(inputFile, nameKey=nameKey, keepHs=keepHs)

def globInputFiles(inputDir, pattern):
    '''Files directly inside inputDir whose name matches pattern.'''
    inputDir = os.path.realpath(inputDir)
    if not os.path.isdir(inputDir):
        print('No such input directory {}'.format(inputDir))
        return []

    inFiles = []
    for name in sorted(os.listdir(inputDir)):
        if name.startswith('.') or not fnmatch.fnmatch(name, pattern):
            continue

        inFile = os.path.join(inputDir, name)
        if os.path.isfile(inFile):
            inFiles.append(inFile)

    return inFiles

def writeSmilesCsv(mols, outDir, outName, outBase, nameKey):
    '''Write the molecules as a csv of smiles and names'''
    outFile = os.path.abspath(os.path.join(outDir, '{}.csv'.format(outName or 'molecules')))
    allRows = []
    for i, mol in enumerate(mols):
        if not mol:
            continue

        try:
            Chem.SanitizeMol(mol)
            mol = Chem.RemoveHs(mol)
        except Exception as e:
            print('Could not sanitize molecule {}: {}'.format(i + 1, e))
            continue

        smiles = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
        allRows.append({"smiles": smiles, "name": getMolName(mol, nameKey, outBase, i)})

    with open(outFile, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["smiles", "name"])
        writer.writeheader()
        writer.writerows(allRows)

    print("SMILES CSV saved to:", outFile)

def getMolWriter(outFormat):
    '''Writer class and file extension for the requested output format'''
    if outFormat in ('smi', 'smiles'):
        return Chem.SmilesWriter, 'smi'

    if outFormat == 'pdb':
        return Chem.PDBWriter, 'pdb'

    return Chem.SDWriter, 'sdf'

def writeSingleFile(mols, writter, outFile):
    '''Write every molecule into one file'''
    with writter(outFile) as f:
        for mol in mols:
            if mol:
                f.write(mol)

def writeFilePerMol(mols, writter, ext, outDir, outBase, nameKey):
    '''Write one file per molecule, named after the molecule'''
    for i, mol in enumerate(mols):
        if not mol:
            continue

        molName = getMolName(mol, nameKey, outBase, i).replace('/', '-').replace(' ', '_')
        outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(molName, ext)))
        with writter(outFile) as f:
            f.write(mol)

def convertFile(inputFile, outFormat, outDir, singleOutFile, outName, outBase,
                make3d, nameKey, nt, keepHs):
    '''Convert one molecule file.'''
    mols, nameKey = loadInputFiles(inputFile, nameKey, keepHs)
    if not mols:
        return

    if make3d:
        mols = performBatchThreading(make3DCoords, mols, nt, cloneItem=False,
                                     errBase=os.path.join(outDir, 'errors3D'))

    if outFormat in ["smiles_csv", "csv", "smi_csv"]:
        writeSmilesCsv(mols, outDir, outName, outBase, nameKey)
        return

    writter, ext = getMolWriter(outFormat)
    if singleOutFile:
        outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(outName.replace(' ', '_'), ext)))
        writeSingleFile(mols, writter, outFile)
    else:
        writeFilePerMol(mols, writter, ext, outDir, outBase, nameKey)

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
    parser.add_argument('--multiFiles', default=False, action='store_true', help='Multiple files to convert')
    parser.add_argument('-iD', '--inputDir', type=str, help='Input molecule files directory if multiFiles')
    parser.add_argument('-pat', '--pattern', type=str, required=False, default='',
                        help='Input molecule files pattern if multiFiles')
    parser.add_argument('-i', '--inputFilename', default='', type=str, help='Input molecule file')
    parser.add_argument('-of', '--outputFormat', type=str, required=False, default='sdf', help='Output format')
    parser.add_argument('-o', '--outputName', type=str, required=False, help='Output name')
    parser.add_argument('-ob', '--outputBase', type=str, required=False, help='Output basename for multiple outputs')
    parser.add_argument('-od', '--outputDir', type=str, required=False, help='Output directory')
    parser.add_argument('--make3D', default=False, action='store_true', help='Optimize 3D coordinates')
    parser.add_argument('--keepHs', default=False, action='store_true',
                        help='Keep the hydrogens present in the file (RDKit removes them by default)')
    parser.add_argument('--overWrite', default=False, action='store_true', help='Overwrite output')
    parser.add_argument('--nameKey', default='', type=str, required=False, help='molecule name key in file')
    parser.add_argument('-nt', '--nthreads', default=1, type=int, required=False, help='Number of threads')

    args = parser.parse_args()
    inputFile, outFormat = args.inputFilename, args.outputFormat
    outFormat = outFormat if not outFormat.startswith('.') else outFormat[1:]

    if args.outputName:
        singleOutFile, outName = True, os.path.splitext(args.outputName)[0]
    else:
        singleOutFile, outName = False, None

    if args.outputDir:
        outDir = args.outputDir
    elif args.multiFiles:
        outDir = args.inputDir
    else:
        outDir = os.path.dirname(inputFile)

    make3d = args.make3D
    nameKey = args.nameKey
    nt = args.nthreads

    outBase = args.outputBase if args.outputBase else 'molecule'
    if args.multiFiles:
        # One output per input file, named after it, exactly like obabel_IO.py --multiFiles
        for inFile in globInputFiles(args.inputDir, args.pattern):
            convertFile(inFile, outFormat, outDir, True, os.path.splitext(os.path.basename(inFile))[0],
                        outBase, make3d, nameKey, nt, args.keepHs)
    else:
        convertFile(inputFile, outFormat, outDir, singleOutFile, outName, outBase,
                    make3d, nameKey, nt, args.keepHs)
