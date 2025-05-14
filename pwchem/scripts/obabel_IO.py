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

'''Script to convert molecule files using the openbabel pybel module in the plip-env'''

import sys, os, argparse, glob

from openbabel import pybel

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'utils'))
from scriptUtils import *

def getPDBMols(inFile):
    mols = []
    pdbStrs = divideMultiPDB(inFile)
    for i, pdbStr in enumerate(pdbStrs):
        tFile = 'molecule_{}'.format(i)
        with open(tFile, 'w') as f:
            f.write(pdbStr)
        mols += list(pybel.readfile('pdb', tFile))
        os.remove(tFile)
    return mols

def parseMols(inFile):
    inFormat = os.path.splitext(inFile)[1][1:]
    if inFormat == 'pdb':
        mols = getPDBMols(inFile)
    else:
        mols = pybel.readfile(inFormat, inFile)
    return mols

def make3DCoords(mols, mols3dLists, it):
    '''Optimize the 3D coordinates of a rdkit molecule'''
    for i, mol in enumerate(mols):
        mol.make3D()
        mols3dLists[it].append(mol)
    return mols3dLists[it]

def oBabelConversion(inputFile, outFormat, singleOutFile, outDir, outName=None, outBase=None,
                     overW=True, make3d=False, nameKey=None, nt=1):
    mols = list(parseMols(inputFile))
    inFormat = os.path.splitext(inputFile)[1][1:]
    if inFormat in ['smi', 'smiles'] or make3d:
        mols = performBatchThreading(make3DCoords, mols, nt, cloneItem=False)

    if singleOutFile:
        outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(outName, outFormat)))
        outObj = pybel.Outputfile(outFormat, outFile, overwrite=overW)

        for i, mol in enumerate(mols):
            if outName:
                if len(mols) == 1:
                    mol.title = outName
                else:
                    mol.title = f'{outName}_{i}'
            outObj.write(mol)
        outObj.close()
        return outFile

    else:
        names = []
        for i, mol in enumerate(mols):
            if nameKey and nameKey in mol.data:
                molName = mol.data[nameKey]
            elif mol.title and not mol.title in names and not outBase:
                molName = os.path.splitext(os.path.basename(mol.title))[0]
            else:
                if not outBase:
                    outBase = 'molecule'
                molName = '{}_{}'.format(outBase, i + 1)
            names.append(molName)
            mol.write(outFormat, os.path.abspath(os.path.join(outDir, '{}.{}'.format(molName, outFormat))),
                      overwrite=overW)

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
    parser.add_argument('-of', '--outputFormat', type=str, default='mol2', help='Output format')
    parser.add_argument('-o', '--outputName', type=str, required=False, help='Output name')
    parser.add_argument('-ob', '--outputBase', type=str, required=False, help='Output basename for multiple outputs')
    parser.add_argument('-od', '--outputDir', type=str, required=False, help='Output directory')
    parser.add_argument('--make3D', default=False, action='store_true', help='Optimize 3D coordinates')
    parser.add_argument('--overWrite', default=True, action='store_true', help='Overwrite output')
    parser.add_argument('--nameKey', type=str, required=False, help='molecule name key in file')
    parser.add_argument('-nt', '--nthreads', default=1, type=int, required=False, help='Number of threads')

    args = parser.parse_args()

    inputFile = args.inputFilename
    if args.outputDir:
        outDir = args.outputDir
    else:
        outDir = os.path.dirname(inputFile)

    outFormat = args.outputFormat
    outFormat = outFormat if not outFormat.startswith('.') else outFormat[1:]

    outBase = args.outputBase
    overW = args.overWrite
    make3d = args.make3D
    nameKey = args.nameKey
    nt = args.nthreads

    if not args.multiFiles:

        if args.outputName:
            singleOutFile, outName = True, os.path.splitext(args.outputName)[0]
        else:
            singleOutFile, outName = False, None

        oBabelConversion(inputFile, outFormat, singleOutFile, outDir, outName, outBase, overW, make3d, nameKey, nt)

    else:
        inputDir, pattern = args.inputDir, args.pattern
        pattern = os.path.join(inputDir, pattern)
        for inFile in glob.glob(pattern):
            outName = os.path.splitext(os.path.basename(inFile))[0]
            oBabelConversion(inFile, outFormat, True, outDir, outName, outBase, overW, make3d, nameKey, nt)