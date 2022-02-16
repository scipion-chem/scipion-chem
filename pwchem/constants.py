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

PYMOL_HOME = 'PYMOL_HOME'

PML_STR = '''from pymol import cmd,stored
load {}
#select pockets, resn STP
stored.list=[]
cmd.iterate("(resn STP)","stored.list.append(resi)")	#read info about residues STP

aux = list(map(int, stored.list))
aux.sort()
stored.list = list(map(str, aux))
#print(stored.list)
lastSTP=stored.list[-1]	#get the index of the last residu
hide lines, resn STP

#show spheres, resn STP
for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.color(my_index+1,"pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.3","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))'''

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

aux = list(map(int, stored.list))
aux.sort()
stored.list = list(map(str, aux))
#print(stored.list)
lastSTP=stored.list[-1]	#get the index of the last residu
hide lines, resn STP

#show spheres, resn STP
for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.color(my_index+1,"pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.3","pocket"+str(my_index))
for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))

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
