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
for my_index in range(1,int(lastSTP)+1): cmd.color(my_index,"pocket"+str(my_index))
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
