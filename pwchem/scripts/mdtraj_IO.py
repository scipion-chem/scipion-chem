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
import parmed


def convertTopology(topFile, outFile, gromacs_topdir=None):
    """Convert topology file using ParmEd"""
    print(f"Converting topology: {topFile} -> {outFile}")

    if topFile.endswith('.top') and gromacs_topdir:
        parmed.gromacs.GROMACS_TOPDIR = gromacs_topdir

    topology = parmed.load_file(topFile)
    topology.save(outFile)

def convertSystem(sysFile, outFile):
    """Convert system file using MDTraj"""
    print(f"Converting system: {sysFile} -> {outFile}")
    system = mdtraj.load(sysFile, top=sysFile)
    system.save(outFile)

def convertTraj(trjFile, sysFile, outFile):
    """Convert trajectory file using MDTraj"""
    print(f"Converting trajectory: {trjFile} -> {outFile}")
    traj = mdtraj.load(trjFile, top=sysFile)
    traj.save(outFile)

if __name__ == "__main__":
    '''
    Use: python <scriptName> -s <sysFile> [-top <topFile>] [-t <trjFile>] 
         [-os <outputSystem>] [-otop <outputTopology>] [-otj <outputTraj>]
    '''
    parser = argparse.ArgumentParser(
        description='Convert MD system, topology, and trajectory files')
    parser.add_argument('-s', '--inputSystem', type=str,
                        help='Input system file (e.g., .pdb, .gro)')
    parser.add_argument('-top', '--inputTopology', type=str, default='',
                        help='Input topology file (e.g., .top, .prmtop, .psf). If not provided, uses -s for topology operations')
    parser.add_argument('-t', '--inputTraj', type=str, default='',
                        help='Input trajectory file (e.g., .xtc, .dcd, .trr)')

    # Output files (only run conversion if provided)
    parser.add_argument('-os', '--outputSystem', type=str, default='',
                        help='Output system filename (triggers system conversion)')
    parser.add_argument('-otop', '--outputTopology', type=str, default='',
                        help='Output topology filename (triggers topology conversion)')
    parser.add_argument('-otj', '--outputTraj', type=str, default='',
                        help='Output trajectory filename (triggers trajectory conversion)')

    # GROMACS topology directory for .top files
    parser.add_argument('-gtopdir', '--gromacs_topdir', type=str, default=None,
                        help='GROMACS topology directory (needed for loading .top files)')

    args = parser.parse_args()

    if args.outputSystem:
        convertSystem(args.inputSystem, args.outputSystem)

    if args.outputTopology:
        convertTopology(args.inputTopology, args.outputTopology, args.gromacs_topdir)

    if args.outputTraj:
        convertTraj(args.inputTraj, args.inputSystem, args.outputTraj)

