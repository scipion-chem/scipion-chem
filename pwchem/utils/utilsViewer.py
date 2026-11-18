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
    # Build complex groups, in case mols have different receptor files
    cGroups, uNames = {}, {}
    for mol in mols:
        molFile = mol.getPoseFile() if pose else mol.getFileName()
        recFile = viewer.getReceptorFile(mol) if addTarget and pose else None
        cGroups = addToDic(cGroups, os.path.abspath(recFile) if recFile else 'all', molFile)
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


def buildPMLCovalentStr(covDic, disable=True):
    '''PyMOL script for covalently docked poses.

    :param covDic: {objectName: covalent complex file}
    :param disable: hide every object but the first, so poses do not overlap
    '''
    pmlStr = 'set cartoon_transparency, 0.65\n'
    for i, (objName, covFile) in enumerate(covDic.items()):
        pmlStr += f'load {os.path.abspath(covFile)}, {objName}\n'
        pmlStr += f'hide everything, {objName}\n'
        pmlStr += f'show cartoon, {objName} and polymer\n'
        pmlStr += f'show sticks, {objName} and not polymer and not solvent\n'
        pmlStr += f'show sticks, byres ({objName} and polymer within 5 of ' \
                  f'({objName} and resn LIG))\n'
        pmlStr += f'color grey70, {objName} and polymer and elem C\n'
        pmlStr += f'color cyan, {objName} and resn LIG and elem C\n'
        if disable and i != 0:
            pmlStr += f'disable {objName}\n'

    if covDic:
        pmlStr += f'orient {next(iter(covDic))} and resn LIG\n'
    return pmlStr

def buildChimeraDockingGroupsStr(viewer, mols, addTarget=True, pose=True, disable=True,
                                 vinaDock=False):
    """ChimeraX counterpart of buildPMLDockingGroupsStr: each receptor is shown with the poses
    docked into it, grouped together as a pymol group is.
    """
    cGroups, uNames = {}, {}
    for mol in mols:
        if vinaDock:
            molFile, molName = viewer.getDockFile(mol), mol.getUniqueName(pose=False, dock=False)
        else:
            molFile, molName = mol.getPoseFile() if pose else mol.getFileName(), mol.getUniqueName()

        recFile = viewer.getReceptorFile(mol) if addTarget and pose else None
        cGroups = addToDic(cGroups, os.path.abspath(recFile) if recFile else 'all', molFile)
        uNames[molFile] = molName

    chimStr, mi, cGroupLims = '', 1, []
    for recFile, molFiles in cGroups.items():
        firstId = mi
        if recFile != 'all':
            chimStr += f"run(session, 'open {recFile}')\n"
            mi += 1

        molFiles = natural_sort(set(molFiles))
        for i, molFile in enumerate(molFiles):
            chimStr += f"run(session, 'open {os.path.abspath(molFile)} name {uNames[molFile]}')\n"
            if disable and i != 0:
                chimStr += f"run(session, 'hide #{mi} models')\n"
            mi += 1

        if recFile != 'all':
            ci = len(cGroupLims) + 1
            complexName = f'{getBaseName(recFile)}_{ci}' if len(molFiles) > 1 \
                else f'{uNames[molFiles[0]]}_{ci}'
            cGroupLims.append((firstId, mi - 1, complexName))

    for ci, (firstId, lastId, complexName) in enumerate(cGroupLims):
        if firstId < lastId:
            chimStr += f"""run(session, 'rename #{firstId}-{lastId} "{complexName}" id #{mi + ci}')\n"""
            if disable and ci != 0:
                chimStr += f"run(session, 'hide #{mi + ci} models')\n"
    return chimStr


def buildChimeraCovalentStr(covDic, disable=True):
    """ChimeraX counterpart of buildPMLCovalentStr.

    :param covDic: {objectName: covalent complex file}
    :param disable: hide every complex but the first, so poses do not overlap
    """
    chimStr = ''
    for i, (objName, covFile) in enumerate(covDic.items()):
        chimStr += f"run(session, 'open {os.path.abspath(covFile)} name {objName}')\n"
        if disable and i != 0:
            chimStr += f"run(session, 'hide #{i + 1} models')\n"
    return chimStr


def writeChimeraFile(chimFile, chimStr, vinaDock=False):
    with open(chimFile, 'w') as f:
        f.write('from chimerax.core.commands import run\n')
        f.write(f"run(session, 'cd {os.getcwd()}')\n")
        f.write("run(session, 'cofr 0,0,0')\n")  # set center of coordinates
        f.write(chimStr)
        if vinaDock:
            f.write("run(session, 'viewdockx')\n")
    return chimFile, chimStr

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