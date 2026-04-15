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

'''Script to calculate the RMSD of a set of molecules pair to pair using SPYRMSD.
It iterates over the molecule files and calculates the RMSD to each of the downstream files.
The list, which represents the flatten upper part of a distance matrix is printed'''

import sys, os, argparse, shutil, gzip, threading

from spyrmsd import io
from spyrmsd.rmsd import rmsdwrapper

sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'utils'))
from scriptUtils import *

def runRMSDCalculation(ref, mols, symmetry, center, minimize, strip):
    RMSDlist = rmsdwrapper(ref, mols,
                           symmetry=symmetry, center=center, minimize=minimize, strip=strip, cache=False)
    return RMSDlist

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
    parser.add_argument('molecules',  type=str, nargs="+", help="Input molecule file(s)")
    parser.add_argument("-m", "--minimize", action="store_true", help="Minimize (fit)")
    parser.add_argument("-c", "--center", action="store_true", help="Center molecules at origin")
    parser.add_argument("--hydrogens", action="store_true", help="Keep hydrogen atoms")
    parser.add_argument("-n", "--nosymm", action="store_false", help="No graph isomorphism")

    args = parser.parse_args()

    try:
        mols = [mol for molfile in args.molecules for mol in io.loadallmols(molfile)]
    except OSError:
        print("ERROR: Molecule file(s) not found.", file=sys.stderr)
        exit(-1)

    rmsds = []
    for i in range(len(mols)-1):
        ref, targets = mols[i], mols[i + 1:]
        try:
            RMSDlist = runRMSDCalculation(ref, targets, symmetry=args.nosymm, center=args.center, minimize=args.minimize,
                                          strip=not args.hydrogens)
        except:
            try:
                RMSDlist = runRMSDCalculation(ref, targets, symmetry=False, center=args.center, minimize=args.minimize,
                                              strip=not args.hydrogens)
            except:
                RMSDlist = ['1000'] * len(targets)


        rmsds += RMSDlist

    print(rmsds)
