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

'''Calculates the SASA of an atomic structure file. Biopython 1.79 needed (plip-env)'''

import sys, os, shutil

from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.SASA import ShrakeRupley


if __name__ == "__main__":
    '''Use: python <scriptName> <structFile> <outFile>'''
    structFile, outFile = sys.argv[1], sys.argv[2]

    if structFile.endswith('.pdb') or structFile.endswith('.ent'):
        p = PDBParser(QUIET=1)
    elif structFile.endswith('.cif'):
        p = MMCIFParser(QUIET=1)
    struct = p.get_structure("SASAstruct", structFile)

    sr = ShrakeRupley()
    sr.compute(struct, level="R")

    with open(outFile, 'w') as f:
        for model in struct:
            modelID = model.get_id()
            for chain in model:
                chainID = chain.get_id()
                for residue in chain:
                    resId = residue.get_id()[1]
                    f.write('{}:{}\t{}\n'.format(chainID, resId, residue.sasa))

