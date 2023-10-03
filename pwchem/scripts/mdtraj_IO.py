# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors: Daniel Del Hoyo (ddelhoyo@cnb.csic.es)
# *          James M. Krieger (jmkrieger@cnb.csic.es)
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
# *  e-mail address 'ddelhoyo@cnb.csic.es'
# *
# **************************************************************************

import argparse, os
import mdtraj

def convertSystem(sysFile, outFile):
    system = mdtraj.load(sysFile, top=sysFile)
    system.save(outFile)

def convertTraj(trjFile, sysFile, outFile):
    traj = mdtraj.load(trjFile, top=sysFile)
    traj.save(outFile)

if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputFilename <sysFile> -o/--outputName <outputName> 
    -t/--inputTraj [<trjFile>]
    If only an MD system file is provided, then the script will convert the system to the required format,
    based on outputName.
    If a trajectory file is provided too, then the script will convert the trajectory to the required format,
    based on outputName, using the system file to load it.
    '''
    parser = argparse.ArgumentParser(description='Handles the IO for molecule files using openbabel')
    parser.add_argument('-s', '--inputSystem', type=str, help='Input system file')
    parser.add_argument('-o', '--outputName', type=str, help='Output name')
    parser.add_argument('-t', '--inputTraj', type=str, help='Input trajectory file', 
                        required=False, default='')

    args = parser.parse_args()
    sysFile, outFile, trjFile = args.inputSystem, args.outputName, args.inputTraj

    if trjFile == '':
        convertSystem(sysFile, outFile)
    else:
        convertTraj(trjFile, sysFile, outFile)

