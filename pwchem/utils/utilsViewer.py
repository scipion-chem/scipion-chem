# **************************************************************************
# *
# * Authors:     Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

import os
from ..utils import addToDic, getBaseName, natural_sort

def buildPMLDockingSingleStr(viewer, molFile, molName, addTarget=True, disable=True):
    pmlStr = ''
    if addTarget:
        pmlStr = 'load {}\n'.format(os.path.abspath(viewer.protocol.getOriginalReceptorFile()))

    disableStr = ''
    if disable:
        disableStr = '\ndisable {}'.format(molName)

    pdbFile = os.path.abspath(molFile)
    pmlStr += 'load {}, {}{}\nhide spheres, {}\nshow sticks, {}\n'. \
        format(pdbFile, molName, disableStr, molName, molName)
    return pmlStr

def buildSchShowStr(name, isReceptor=True):
    schStrBase = ', mimic=1, object_props=*, atom_props=*'
    schStr = f'{schStrBase}\nhide lines, {name}'
    if isReceptor:
        schStr += f'\nshow cartoon, {name}'
    else:
        schStr += f'\nshow sticks, {name}'
    return schStr

def buildPMLDockingGroupsStr(viewer, mols, addTarget=True, pose=True, disable=True):
    if addTarget and pose:
        genRecFile = os.path.abspath(viewer.protocol.getOriginalReceptorFile())

    # Build complex groups, in case mols have different receptor files
    cGroups, uNames = {}, {}
    for mol in mols:
        molFile = mol.getPoseFile() if pose else mol.getFileName()
        if addTarget and pose:
            if mol.getProteinFile():
                recFile = mol.getProteinFile()
            else:
                recFile = genRecFile

            cGroups = addToDic(cGroups, recFile, molFile)
        else:
            cGroups = addToDic(cGroups, 'all', molFile)
        uNames[molFile] = mol.getUniqueName()

    pmlStr, ci = '', 1

    for recFile, molFiles in cGroups.items():
        gNames = []
        if recFile != 'all':
            gNames += [f'{getBaseName(recFile)}_{ci}']
            schStr = '' if '.mae' not in recFile else buildSchShowStr(gNames[-1], True)
            pmlStr += f'load {os.path.abspath(recFile)}, {gNames[-1]}{schStr}\n'

        molFiles = natural_sort(set(molFiles))
        for mi, molFile in enumerate(molFiles):
            gNames.append(uNames[molFile])
            schStr = '' if '.mae' not in molFile else buildSchShowStr(gNames[-1], False)
            disableStr = f'\ndisable {gNames[-1]}' if (disable and mi != 0) else ''
            pmlStr += f'load {os.path.abspath(molFile)}, {gNames[-1]}{schStr}{disableStr}\n'

        if recFile != 'all':
            complexName = f'{gNames[0]}_{ci}' if len(gNames) > 2 else f'{gNames[-1]}_{ci}'
            disableStr = f'\ndisable {complexName}' if (disable and ci != 1) else ''
            pmlStr += f'group {complexName}, {" ".join(gNames)} add{disableStr}\n'
            ci += 1
    return pmlStr


def writePmlFile(pmlFile, pmlStr):
    with open(pmlFile, 'w') as f:
        f.write(pmlStr)
        f.write('zoom')
    return pmlFile, pmlStr

def sortMolsByUnique(mols):
    uniques, newMols = [], []
    for mol in mols:
        uniq = mol.clone().getUniqueName()
        if not uniq in uniques:
            uniques.append(uniq)
            newMols.append(mol.clone())

    zipped_lists = sorted(zip(uniques, newMols))
    uniques, mols = zip(*zipped_lists)
    return mols

def getPmlsDir(protocol):
    pmlsDir = protocol._getExtraPath('pmls')
    if not os.path.exists(pmlsDir):
        os.mkdir(pmlsDir)
    return pmlsDir