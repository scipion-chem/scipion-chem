# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
# *  e-mail address 'scipion@cnb.csic.es'
# *
# **************************************************************************

import sys
from rdkit import Chem
from rdkit.Chem import Draw

if __name__ == "__main__":
    if len(sys.argv)==1:
        print("Usage: python3 rdkitUtils.py [options]")
        print("   draw <smallMoleculeFile> <pngFile>: Small molecules .smi, .sdf, .mae, .mol2. .pdb")
    elif sys.argv[1]=="draw":
        fnIn = sys.argv[2]
        fnOut = sys.argv[3]
        if fnIn.endswith('.smi'):
            smile=open(fnIn).readlines()[0]
            smile=smile.split()[0]
            Draw.MolToFile(Chem.MolFromSmiles(smile),fnOut)
        elif fnIn.endswith('.sdf'):
            supplier = Chem.rdmolfiles.SDMolSupplier(fnIn)
            for molecule in supplier:
                Draw.MolToFile(molecule, fnOut)
        elif fnIn.endswith('.mae') or fnIn.endswith('.maegz'):
            supplier = Chem.rdmolfiles.MaeMolSupplier(fnIn)
            for molecule in supplier:
                Draw.MolToFile(molecule, fnOut)
        elif fnIn.endswith('.mol2'):
            supplier = Chem.rdmolfiles.MolFromMol2Block(fnIn)
            for molecule in supplier:
                Draw.MolToFile(molecule, fnOut)
        elif fnIn.endswith('.pdb'):
            molecule = Chem.rdmolfiles.MolFromPDBFile(fnIn)
            Draw.MolToFile(molecule, fnOut)
