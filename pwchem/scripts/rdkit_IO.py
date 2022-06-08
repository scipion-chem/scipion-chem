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

import sys, os, argparse, shutil, gzip
from rdkit import Chem

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
    parser.add_argument('-of', '--outputFormat', type=str, required=False, default='mol2', help='Output format')
    parser.add_argument('-o', '--outputName', type=str, required=False, help='Output name')
    parser.add_argument('-ob', '--outputBase', type=str, required=False, help='Output basename for multiple outputs')
    parser.add_argument('-od', '--outputDir', type=str, required=False, help='Output directory')
    parser.add_argument('--overWrite', default=False, action='store_true', help='Overwrite output')

    args = parser.parse_args()
    inputFile, outFormat = args.inputFilename, args.outputFormat
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

    if inFormat == 'maegz':
        newInputFile = inputFile.replace('.maegz', '.mae')
        with gzip.open(inputFile) as fIn:
            with open(newInputFile, 'wb') as f:
                shutil.copyfileobj(fIn, f)
        inputFile, inFormat = newInputFile, 'mae'

    if inFormat == 'mae':
        if singleOutFile:
            outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(outName, 'sdf')))
            with Chem.SDWriter(outFile) as f:
                for mol in Chem.MaeMolSupplier(inputFile):
                    f.write(mol)

        else:
            if args.outputBase:
                outBase = args.outputBase
            else:
                outBase = 'molecule'

            for i, mol in enumerate(Chem.MaeMolSupplier(inputFile)):
                molName = '{}_{}'.format(outBase, i+1)
                outFile = os.path.abspath(os.path.join(outDir, '{}.{}'.format(molName, 'sdf')))
                with Chem.SDWriter(outFile) as f:
                    f.write(mol)
    else:
        print('Script currently only prepared to convert from mae / maegz files to sdf using rdkit')


