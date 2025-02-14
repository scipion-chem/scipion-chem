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

import sys, os
from rdkit import Chem
from rdkit.Chem import Draw

def parseMoleculeFile(molFile):
    if molFile.endswith('.mol2'):
        mol = Chem.MolFromMol2File(molFile)
    elif molFile.endswith('.mol'):
        mol = Chem.MolFromMolFile(molFile)
    elif molFile.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(molFile)
    elif molFile.endswith('.smi'):
        f = open(molFile, "r")
        firstline = next(f)
        mol = Chem.MolFromSmiles(str(firstline))
    elif molFile.endswith('.sdf'):
        suppl = Chem.SDMolSupplier(molFile)
        for mol in suppl:
            break
    else:
        mol = Chem.MolFromSmiles(molFile)

    return mol

def getMolFilesDic(molFiles):
    molsDict = {}
    for molFile in molFiles:
        m = parseMoleculeFile(molFile)
        molsDict[m] = molFile

    mols = list(molsDict.keys())
    return molsDict, mols

def writeMol(mol, outFile, cid=-1, setName=False):
    w = Chem.SDWriter(outFile)
    molName = os.path.split(os.path.splitext(outFile)[0])[-1]
    if setName:
        mol.SetProp('_Name', molName)
    w.write(mol, cid)
    w.close()

def drawMolecule(mol, fnOut):
    Draw.MolToFile(mol, fnOut)

if __name__ == "__main__":
    if len(sys.argv)==1:
        print("Usage: python3 rdkitUtils.py [options]")
        print("   draw <smallMoleculeFile> <pngFile>: Small molecules .smi, .sdf, .mae, .mol2. .pdb")
    elif sys.argv[1]=="draw":
        fnIn = sys.argv[2]
        fnOut = sys.argv[3]

        mol = parseMoleculeFile(fnIn)
        drawMolecule(mol)

