# coding: latin-1
# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano
# *
# * [1] uam, madrid, Spain
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

#Constant dictionaries
MGL_DIC =       {'name': 'mgltools',    'version': '1.5.7',         'home': 'MGL_HOME'}
JCHEM_DIC =     {'name': 'jchempaint',  'version': '3.2.0',         'home': 'JCHEM_HOME'}
OPENBABEL_DIC = {'name': 'openbabel',   'version': '2.2',           'home': 'OPENBABEL_HOME'}
ALIVIEW_DIC =   {'name': 'aliview',     'version': '1.28',          'home': 'ALIVIEW_HOME'}
SHAPEIT_DIC =   {'name': 'shape-it',    'version': '2.0.0',         'home': 'SHAPEIT_HOME'}
VMD_DIC =       {'name': 'vmd',         'version': '1.9.3',         'home': 'VMD_CHEM_HOME'}
RDKIT_DIC =     {'name': 'rdkit',       'version': '2022.09.1', 'home': 'RDKIT_HOME'}
BIOCONDA_DIC =  {'name': 'bioconda',    'version': '1.0'}
MDTRAJ_DIC =    {'name': 'mdtraj',      'version': '1.9.8',         'home': 'MDTRAJ_HOME'}
DEAP_DIC =      {'name': 'deap',        'version': '1.4',           'home': 'DEAP_HOME'}
RANX_DIC =     {'name': 'ranx',      'version': '0.3.20',         'home': 'RANKX_HOME'}

#Autoligand
POCKET_ATTRIBUTES_MAPPING = {'Pocket Score': 'score', 'Drug Score': 'druggability', 'nPoints': 'nPoints',
                      'Total Volume': 'volume', 'Total Energy per Vol': 'energy', 'class': 'class',
                      'contactAtoms': 'contactAtoms', 'contactResidues': 'contactResidues'}
#Fpocket
POCKET_ATTRIBUTES_MAPPING.update({'Pocket Score': 'score', 'Drug Score': 'druggability', 'Number of alpha spheres': 'nPoints',
                      'Pocket volume (Monte Carlo)': 'volume', 'class': 'class',
                      'contactAtoms': 'contactAtoms', 'contactResidues': 'contactResidues'})
#P2Rank
POCKET_ATTRIBUTES_MAPPING.update({'score': 'score', 'sas_points': 'nPoints', 'class': 'class',
                      'surf_atom_ids': 'contactAtoms', 'residue_ids': 'contactResidues'})
#Sitemap
POCKET_ATTRIBUTES_MAPPING.update({'SiteScore': 'score', 'Dscore': 'druggability', 'size': 'nPoints',
                      'volume': 'volume', 'exposure': 'exposure', 'enclosure': 'enclosure',
                      'contact': 'contact', 'phobic': 'hidrophobic', 'philic': 'hidrophilic',
                      'balance': 'balance', 'don/acc': 'don/acc', 'class': 'class',
                      'contactAtoms': 'contactAtoms', 'contactResidues': 'contactResidues'})

# Pharmacophore
FEATURE_LABELS_SIMPLE = ["Donor", "Acceptor", "Hydrophobe", "Aromatic"]
FEATURE_LABELS_ADVANCED = ["LumpedHydrophobe", "PosIonizable", "NegIonizable", "ZnBinder"]

MAX_MOLS_SET = 'MAX_MOLS_SET'
WARNLIBBIG = f'WARNING: you are about to split an immense library of molecules, ' \
             'which can cause severe disk IO traffic and storage use.\n' \
             f'You can update this value setting the {MAX_MOLS_SET} scipion variable in your scipion.conf file'

PML_STR = '''from pymol import cmd,stored
load {}
#select pockets, resn STP
stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")	#read info about residues STP

aux = {}
aux.sort()

lastSTP=max(list(map(int, stored.list)))	#get the index of the last residu
stored.list = list(map(str, stored.list))
hide lines, resn STP

#show spheres, resn STP
for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(aux[my_index-1]), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.color(my_index+1,"pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.3","pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(aux[my_index-1]))'''

PML_SURF_STR = '''from pymol import cmd,stored
load {}, protein
create ligands, protein and organic
select xlig, protein and organic
delete xlig

hide everything, all

color white, elem c
color bluewhite, protein
show surface, protein

show sticks, ligands
set stick_color, magenta

{}
zoom visible
'''

PML_SURF_EACH = '''set_color pcol{} = {}
select surf_pocket{}, protein and id {} 
set surface_color,  pcol{}, surf_pocket{}\n'''

TCL_STR = '''proc highlighting { colorId representation id selection } {
   puts "highlighting $id"
   mol representation $representation
   mol material "Diffuse" 
    mol color $colorId
   mol selection $selection
   mol addrep $id
}

set id [mol new %s type pdb]
mol delrep top $id
highlighting Name "Lines" $id "protein"
highlighting Name "Licorice" $id "not protein and not resname STP"
highlighting Element "NewCartoon" $id "protein"
set id [mol new %s type pqr]
                        mol selection "all" 
                         mol material "Glass3" 
                         mol delrep top $id 
                         mol representation "QuickSurf 0.3" 
                         mol color ResId $id 
                         mol addrep $id 
highlighting Index "Points 1" $id "resname STP"
display rendermode GLSL'''#.format(proteinHETATMFile, proteinPQRFile)

FUNCTION_GRID_BOX = '''from pymol.cgo import *
from pymol import cmd
from random import randint

#############################################################################
#
# Acknowledgement:
# This script was written based on the drawgridbox by Cunliang Geng
#
#############################################################################
python
def drawgridbox(center, size, gridName=None, nx=10, ny=10, nz=10, padding=0.0, lw=2.0, r=1.0, g=1.0, b=1.0):
    """
    DESCRIPTION
        Given a center and a side size, draw a grid box around it.

    USAGE:
        drawgridbox [center, size, [nx, [ny, [nz, [padding, [lw, [r, [g, b]]]]]]]]

    PARAMETERS:
        center,       center of the grid box

        size,         size of the sides for each direction x, y, z

        nx,           number of grids on axis X
                      defaults to 10

        ny,           number of grids on axis Y
                      defaults to 10

        nz,           number of grids on axis Z
                      defaults to 10

        padding,      defaults to 0

        lw,           line width
                      defaults to 2.0

        r,            red color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        g,            green color component, valid range is [0.0, 1.0]
                      defaults to 1.0

        b,            blue color component, valid range is [0.0, 1.0]
                      defaults to 1.0

    RETURNS
        string, the name of the CGO box

    NOTES
        * This function creates a randomly named CGO grid box. The user can
        specify the number of grids on X/Y/Z axis, the width of the lines,
        the padding and also the color.
    """

    minX, minY, minZ = center[0] - size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2
    maxX, maxY, maxZ = center[0] + size[0]/2, center[1] + size[1]/2, center[2] + size[2]/2

    print("Box dimensions (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))

    minX, minY, minZ = minX - float(padding), minY - float(padding), minZ - float(padding)
    maxX, maxY, maxZ = maxX + float(padding), maxY + float(padding), maxZ + float(padding)
    nX, nY, nZ =int(nx), int(ny), int(nz)
    dX = (maxX-minX)/nX
    dY = (maxY-minY)/nY
    dZ = (maxZ-minZ)/nZ

    if padding != 0:
        print("Box dimensions + padding (%.2f, %.2f, %.2f)" % (maxX-minX, maxY-minY, maxZ-minZ))
    gridbox = [
        LINEWIDTH, float(lw),
        BEGIN, LINES,
        COLOR, float(r), float(g), float(b),
        ]
    for i in range(nX):
        for j in range(nY):
            for k in range(nZ):
                dots= [
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+k*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+k*dZ,
                    VERTEX, minX+i*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+i*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+j*dY, minZ+(k+1)*dZ,
                    VERTEX, minX+(i+1)*dX, minY+(j+1)*dY, minZ+(k+1)*dZ,
                ]
                gridbox += dots
    gridbox.append(END)

    if gridName == None:
      boxName = "gridbox_" + str(randint(0,10000))
      while boxName in cmd.get_names():
          boxName = "gridbox_" + str(randint(0,10000))
    else:
      boxName = gridName

    cmd.load_cgo(gridbox,boxName)
    return boxName
python end

cmd.extend ("drawgridbox", drawgridbox)
'''

FUNCTION_BOUNDING_BOX = '''
from pymol.cgo import *
from pymol import cmd
from random import randint

#############################################################################
#
# Acknowledgement:
# This script was written based on the drawBoundingBox by Jason Vertrees
#
#############################################################################python
def drawBoundingBox(center, size, gridName=None, padding=0.0, linewidth=2.0, r=1.0, g=1.0, b=1.0):
    """                                                                  
    DESCRIPTION                                                          
            Given selection, draw the bounding box around it.          

    USAGE:
            drawBoundingBox [selection, [padding, [linewidth, [r, [g, b]]]]]

    PARAMETERS:
            selection,              the selection to enboxen.  :-)
                                    defaults to (all)

            padding,                defaults to 0

            linewidth,              width of box lines
                                    defaults to 2.0

            r,                      red color component, valid range is [0.0, 1.0]
                                    defaults to 1.0                               

            g,                      green color component, valid range is [0.0, 1.0]
                                    defaults to 1.0                                 

            b,                      blue color component, valid range is [0.0, 1.0]
                                    defaults to 1.0                                

    RETURNS
            string, the name of the CGO box

    NOTES
            * This function creates a randomly named CGO box that minimally spans the protein. The
            user can specify the width of the lines, the padding and also the color.                            
    """

    minX, minY, minZ = center[0] - size[0]/2, center[1] - size[1]/2, center[2] - size[2]/2
    maxX, maxY, maxZ = center[0] + size[0]/2, center[1] + size[1]/2, center[2] + size[2]/2

    if padding != 0.0:
        minX, minY, minZ = minX - float(padding), minY - float(padding), minZ - float(padding)
        maxX, maxY, maxZ = maxX + float(padding), maxY + float(padding), maxZ + float(padding)

    boundingBox = [
    LINEWIDTH, float(linewidth),

    BEGIN, LINES,
    COLOR, float(r), float(g), float(b),

    VERTEX, minX, minY, minZ,  
    VERTEX, minX, minY, maxZ,  
    VERTEX, minX, maxY, minZ,  
    VERTEX, minX, maxY, maxZ,  
    VERTEX, maxX, minY, minZ,  
    VERTEX, maxX, minY, maxZ,  
    VERTEX, maxX, maxY, minZ,  
    VERTEX, maxX, maxY, maxZ,  
    VERTEX, minX, minY, minZ,  
    VERTEX, maxX, minY, minZ,  
    VERTEX, minX, maxY, minZ,  
    VERTEX, maxX, maxY, minZ,  
    VERTEX, minX, maxY, maxZ,  
    VERTEX, maxX, maxY, maxZ,  
    VERTEX, minX, minY, maxZ,  
    VERTEX, maxX, minY, maxZ,  
    VERTEX, minX, minY, minZ,  
    VERTEX, minX, maxY, minZ,  
    VERTEX, maxX, minY, minZ,  
    VERTEX, maxX, maxY, minZ,  
    VERTEX, minX, minY, maxZ,  
    VERTEX, minX, maxY, maxZ,  
    VERTEX, maxX, minY, maxZ,  
    VERTEX, maxX, maxY, maxZ,  

    END
    ]

    if gridName == None:
        boxName = "box_" + str(randint(0,10000))
        while boxName in cmd.get_names():
            boxName = "box_" + str(randint(0,10000))
    else:
        boxName = gridName

    cmd.load_cgo(boundingBox,boxName)
    return boxName
python end

'''

PML_BBOX_STR = '''load {} 
set cartoon_color, white

python
{}
zoom
python end
'''

PML_BBOX_STR_POCK = '''load {} 
set cartoon_color, white

from pymol import cmd,stored
load {}
#select pockets, resn STP
stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")	#read info about residues STP

aux = {}
aux.sort()

stored.list = list(map(str, stored.list))
#print(stored.list)
lastSTP=stored.list[-1]	#get the index of the last residu
hide lines, resn STP

#show spheres, resn STP
for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(aux[my_index-1]), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.color(my_index+1,"pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.3","pocket"+str(aux[my_index-1]))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(aux[my_index-1]))

python
{}
zoom
python end
'''

PML_BBOX_STR_POCKSURF = '''
from pymol import cmd,stored
load {}, protein
create protein and organic
hide everything, all

color white, elem c
color bluewhite, protein
show surface, protein

{}
zoom visible

python
{}
zoom
python end
'''


PML_BBOX_STR_EACH = '''rgb = {}
boxName = drawBoundingBox(center={}, size={}, gridName="{}", r=rgb[0], g=rgb[1], b=rgb[2])
'''


PML_PHARM = '''from pymol.cgo import *
from pymol import cmd
{}
{}    #ligands

python
{}
python end
'''

SPHERE = '''\tALPHA,   {},
\tCOLOR,    {},    {},    {},
\tSPHERE,   {},   {},   {},  {},\n
'''

LOAD_LIGAND = '''load {}, {}\n'''
DISABLE_LIGAND = '''disable {}\n'''

SPHERE_LIST = '''spherelist_{} = [
    {}
    ]

cmd.load_cgo(spherelist_{}, '{}')\n'''

elements_mass = {'H' : 1.008,'HE' : 4.003, 'LI' : 6.941, 'BE' : 9.012,\
                 'B' : 10.811, 'C' : 12.011, 'N' : 14.007, 'O' : 15.999,\
                 'F' : 18.998, 'NE' : 20.180, 'NA' : 22.990, 'MG' : 24.305,\
                 'AL' : 26.982, 'SI' : 28.086, 'P' : 30.974, 'S' : 32.066,\
                 'CL' : 35.453, 'AR' : 39.948, 'K' : 39.098, 'CA' : 40.078,\
                 'SC' : 44.956, 'TI' : 47.867, 'V' : 50.942, 'CR' : 51.996,\
                 'MN' : 54.938, 'FE' : 55.845, 'CO' : 58.933, 'NI' : 58.693,\
                 'CU' : 63.546, 'ZN' : 65.38, 'GA' : 69.723, 'GE' : 72.631,\
                 'AS' : 74.922, 'SE' : 78.971, 'BR' : 79.904, 'KR' : 84.798,\
                 'RB' : 84.468, 'SR' : 87.62, 'Y' : 88.906, 'ZR' : 91.224,\
                 'NB' : 92.906, 'MO' : 95.95, 'TC' : 98.907, 'RU' : 101.07,\
                 'RH' : 102.906, 'PD' : 106.42, 'AG' : 107.868, 'CD' : 112.414,\
                 'IN' : 114.818, 'SN' : 118.711, 'SB' : 121.760, 'TE' : 126.7,\
                 'I' : 126.904, 'XE' : 131.294, 'CS' : 132.905, 'BA' : 137.328,\
                 'LA' : 138.905, 'CE' : 140.116, 'PR' : 140.908, 'ND' : 144.243,\
                 'PM' : 144.913, 'SM' : 150.36, 'EU' : 151.964, 'GD' : 157.25,\
                 'TB' : 158.925, 'DY': 162.500, 'HO' : 164.930, 'ER' : 167.259,\
                 'TM' : 168.934, 'YB' : 173.055, 'LU' : 174.967, 'HF' : 178.49,\
                 'TA' : 180.948, 'W' : 183.84, 'RE' : 186.207, 'OS' : 190.23,\
                 'IR' : 192.217, 'PT' : 195.085, 'AU' : 196.967, 'HG' : 200.592,\
                 'TL' : 204.383, 'PB' : 207.2, 'BI' : 208.980, 'PO' : 208.982,\
                 'AT' : 209.987, 'RN' : 222.081, 'FR' : 223.020, 'RA' : 226.025,\
                 'AC' : 227.028, 'TH' : 232.038, 'PA' : 231.036, 'U' : 238.029,\
                 'NP' : 237, 'PU' : 244, 'AM' : 243, 'CM' : 247, 'BK' : 247,\
                 'CT' : 251, 'ES' : 252, 'FM' : 257, 'MD' : 258, 'NO' : 259,\
                 'LR' : 262, 'RF' : 261, 'DB' : 262, 'SG' : 266, 'BH' : 264,\
                 'HS' : 269, 'MT' : 268, 'DS' : 271, 'RG' : 272, 'CN' : 285,\
                 'NH' : 284, 'FL' : 289, 'MC' : 288, 'LV' : 292, 'TS' : 294,\
                 'OG' : 294}

############### DB IDS #####################

DB_IDS = [{'UniProt': ['A1A4S6', 'A2RUC4']},                      # UniProt
          {'PDB': ['5ni1', '4erf'], 'ChEMBL': ['CHEMBL5061', 'CHEMBL3881']},    # Targets
          {'ChEMBL': ['CHEMBL25', 'CHEMBL353472'], 'BindingDB': ['429291', '4444'],
           'ZINC': ['ZINC480', 'ZINC1019'], 'PubChem': ['98514', '55748', '8739']}]    # Ligands

#  MOLECULAR DYNAMICS

TCL_MD_STR = '''
mol addrep 0
mol new {%s} type {%s} first 0 last -1 step 1 waitfor 1
mol addfile {%s} type {%s} first 0 last -1 step 1 waitfor 1 0

mol color Name
mol representation NewCartoon 0.300000 10.000000 4.100000 0
mol selection protein
mol material Opaque
mol modrep 0 0

mol addrep 0
mol color Name
mol representation Points 1.000000
mol selection hetero within 3 of protein
mol material Opaque
mol modrep 1 0
'''

TCL_MD_LIG_STR = '''
mol addrep 0
mol modstyle 2 0 Licorice 0.300000 12.000000 12.000000
mol modselect 2 0 resname {}
'''

PML_MD_STR = '''load {}
load_traj {}
hide everything, not br. all within 3 of (byres polymer & name CA)
set movie_fps, 15
'''

NORM_STRATEGY = ["None", "min-max", "min-max-inverted", "max", "sum", "rank", "borda"]
SCORE_BASED_METHODS = ["med", "anz"]
RANK_BASED_METHODS = ["isr", "log_isr", "logn_isr", "rrf", "rbc"]

