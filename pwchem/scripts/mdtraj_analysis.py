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
# *  e-mail address 'ddelhoyo@cnb.csic.es'
# *
# **************************************************************************

import argparse
import mdtraj
import matplotlib.pyplot as plt

def getAtomSelection(selStr, ligName='LIG'):
    if selStr == 'CA':
        aSel = 'protein and name CA'
    elif selStr == 'Ligand':
        aSel = f'resname {ligName}'
    else:
        aSel = selStr.lower()
    return aSel

def getSelectionIndex(topo, selStr):
    aSel = getAtomSelection(selStr)
    index = set(topo.topology.select(aSel))
    if not index and selStr == 'Ligand':
        aSel = getAtomSelection(selStr, ligName='UNK')
        index = set(topo.topology.select(aSel))
    return index


if __name__ == "__main__":
    '''Use: python <scriptName> -i/--inputTopology <topFile> -t/--inputTraj [<trjFile>]
    -o/--outputName <outputName> 
    
    Performs several analysis on a MD trajectory using MDTraj
    '''
    parser = argparse.ArgumentParser(description='Handles the IO for molecule files using openbabel')
    parser.add_argument('-i', '--inputFilename', type=str, help='Input system file (PDB or similar), reference')
    parser.add_argument('-o', '--outputName', type=str, help='Output name')
    parser.add_argument('-t', '--inputTraj', type=str, help='Input trajectory file')

    parser.add_argument('-rmsd', default=False, action='store_true',
                        help='Plots the RMSD of the system through the trajectory, with respect to the original system')
    parser.add_argument('-rmsf', default=False, action='store_true',
                        help='Plots the RMSF of the system through the trajectory, with respect to the original system')

    parser.add_argument('-sa', '--selectAtoms', type=str, help='Atoms to select for the analysis')
    parser.add_argument('-ha', '--heavyAtoms', default=False, action='store_true', help='Analysis only on heavy atoms')


    args = parser.parse_args()
    inpFile, outFile, trjFile = args.inputFilename, args.outputName, args.inputTraj

    topo = mdtraj.load(inpFile)
    trajectory = mdtraj.load(trjFile, top=topo)

    # Atom selection
    index = getSelectionIndex(topo, args.selectAtoms)
    if index:
        if args.heavyAtoms:
            index = index.intersection(set([atom.index for atom in topo.topology.atoms if atom.element.symbol != 'H']))

        plt.figure()
        if args.rmsd:
            rmsd = mdtraj.rmsd(trajectory, topo, 0, atom_indices=list(index))
            plt.plot(trajectory.time, rmsd, 'r', label='RMSD')
            plt.title('RMSD')
            plt.xlabel('Simulation time (ps)')
            plt.ylabel('RMSD (nm)')

        elif args.rmsf:
            rmsf = mdtraj.rmsf(trajectory, topo, 0, atom_indices=list(index))
            plt.plot(rmsf, 'r', label='RMSF')
            plt.title('RMSF')
            plt.xlabel('Protein sequence')
            plt.ylabel('RMSF (nm)')

        plt.legend()
        plt.show()

    else:
        print('No atoms found for the atom selection')
