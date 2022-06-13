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

import sys, os, argparse

from openbabel import pybel

def oBabelConversion(inputFile, outFormat, singleOutFile, outDir, outName=None, outBase=None, overW=False):
    inFormat = os.path.splitext(inputFile)[1][1:]

    if singleOutFile:
        outObj = pybel.Outputfile(outFormat, os.path.abspath(os.path.join(outDir, '{}.{}'.format(outName, outFormat))),
                                  overwrite=overW)
        for mol in pybel.readfile(inFormat, inputFile):
            if inFormat in ['smi', 'smiles']:
                mol.make3D()
            outObj.write(mol)
        outObj.close()

    else:
        if inFormat == 'pdb':
            print('Warning: openbabel is not able to parse combo pdb files')

        if not outBase:
            outBase = 'molecule'

        names = []
        for i, mol in enumerate(pybel.readfile(inFormat, inputFile)):
            if inFormat in ['smi', 'smiles']:
                mol.make3D()
            if mol.title and not mol.title in names and not outBase:
                molName = os.path.splitext(os.path.basename(mol.title))[0]
                names.append(molName)
            else:
                molName = '{}_{}'.format(outBase, i + 1)
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
    parser.add_argument('-i', '--inputFilename', type=str, help='Input molecule file')
    parser.add_argument('-of', '--outputFormat', type=str, default='mol2', help='Output format')
    parser.add_argument('-o', '--outputName', type=str, required=False, help='Output name')
    parser.add_argument('-ob', '--outputBase', type=str, required=False, help='Output basename for multiple outputs')
    parser.add_argument('-od', '--outputDir', type=str, required=False, help='Output directory')
    parser.add_argument('--overWrite', default=False, action='store_true', help='Overwrite output')

    args = parser.parse_args()
    inputFile, outFormat = args.inputFilename, args.outputFormat

    if args.outputName:
        singleOutFile, outName = True, os.path.splitext(args.outputName)[0]
    else:
        singleOutFile, outName = False, None

    if args.outputDir:
        outDir = args.outputDir
    else:
        outDir = os.path.dirname(inputFile)

    outBase = args.outputBase
    overW = args.overWrite

    oBabelConversion(inputFile, outFormat, singleOutFile, outDir, outName, outBase, overW)
