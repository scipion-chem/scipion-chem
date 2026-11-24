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
import xml.etree.ElementTree as ET

from ..utils import addToDic, getBaseName, natural_sort

# PLIP interaction containers in its XML report: {container: (type drawn, ligand atom field, distance field)}
PLIP_INTERACTIONS = {
    'hydrophobic_interactions': ('hydrophobic', 'ligcarbonidx', 'dist'),
    'hydrogen_bonds': ('hbond', ('donoridx', 'acceptoridx'), 'dist_h-a'),
    'water_bridges': ('waterbridge', ('donor_idx', 'acceptor_idx'), 'dist_a-w'),
    'salt_bridges': ('saltbridge', 'lig_idx_list', 'dist'),
    'pi_stacks': ('pistacking', 'lig_idx_list', 'centdist'),
    'pi_cation_interactions': ('pication', 'lig_idx_list', 'dist'),
    'halogen_bonds': ('halogen', 'don_idx', 'dist'),
}

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

def getChimeraDockingGroups(viewer, mols, addTarget, pose, vinaDock):
    """Poses grouped by the receptor they were docked into, so that those that share it are shown together:
    ({receptor file: [pose files]}, {pose file: name of the pose}). Without receptor they are all in 'all'."""
    cGroups, uNames = {}, {}
    for mol in mols:
        if vinaDock:
            molFile, molName = viewer.getDockFile(mol), mol.getUniqueName(pose=False, dock=False)
        else:
            molFile, molName = mol.getPoseFile() if pose else mol.getFileName(), mol.getUniqueName()

        recFile = viewer.getReceptorFile(mol) if addTarget and pose else None
        cGroups = addToDic(cGroups, os.path.abspath(recFile) if recFile else 'all', molFile)
        uNames[molFile] = molName
    return cGroups, uNames

def getChimeraGroupOpenStr(recFile, molFiles, uNames, mi, disable):
    """Commands opening a receptor, if there is one, and the poses docked into it.
    return: the commands, the number of the next model and the pose files in the order they were opened
    """
    chimStr = ''
    if recFile != 'all':
        chimStr += f"run(session, 'open {recFile}')\n"
        mi += 1

    molFiles = natural_sort(set(molFiles))
    for i, molFile in enumerate(molFiles):
        chimStr += f"run(session, 'open {os.path.abspath(molFile)} name {uNames[molFile]}')\n"
        if disable and i != 0:
            chimStr += f"run(session, 'hide #{mi} models')\n"
        mi += 1
    return chimStr, mi, molFiles

def getChimeraGroupName(recFile, molFiles, uNames, ci):
    """Name of the group of a receptor and its poses: after the pose when it is the only one, else after the receptor"""
    baseName = getBaseName(recFile) if len(molFiles) > 1 else uNames[molFiles[0]]
    return f'{baseName}_{ci}'

def getChimeraGroupsRenameStr(cGroupLims, mi, disable):
    """Commands grouping the models of each receptor. The groups are hidden but the first."""
    chimStr = ''
    for ci, (firstId, lastId, groupName) in enumerate(cGroupLims):
        if firstId < lastId:
            chimStr += f"""run(session, 'rename #{firstId}-{lastId} "{groupName}" id #{mi + ci}')\n"""
            if disable and ci != 0:
                chimStr += f"run(session, 'hide #{mi + ci} models')\n"
    return chimStr

def buildChimeraDockingGroupsStr(viewer, mols, addTarget=True, pose=True, disable=True,
                                 vinaDock=False):
    """ChimeraX counterpart of buildPMLDockingGroupsStr: each receptor is shown with the poses
    docked into it, grouped together as a pymol group is.
    """
    cGroups, uNames = getChimeraDockingGroups(viewer, mols, addTarget, pose, vinaDock)

    chimStr, mi, cGroupLims = '', 1, []
    for recFile, molFiles in cGroups.items():
        firstId = mi
        groupStr, mi, molFiles = getChimeraGroupOpenStr(recFile, molFiles, uNames, mi, disable)
        chimStr += groupStr
        if recFile != 'all':
            groupName = getChimeraGroupName(recFile, molFiles, uNames, len(cGroupLims) + 1)
            cGroupLims.append((firstId, mi - 1, groupName))

    return chimStr + getChimeraGroupsRenameStr(cGroupLims, mi, disable)


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

def parsePlipXML(xmlFile):
    """Interactions PLIP found for a pose, in the shape drawInteractions.py expects.

    :return: {(residueId, ligandAtomIds): {interactionType: distance}}
    """
    intDic = {}
    for site in ET.parse(xmlFile).getroot().iter('bindingsite'):
        for container in site.find('interactions'):
            if container.tag not in PLIP_INTERACTIONS:
                continue
            intType, idxField, distField = PLIP_INTERACTIONS[container.tag]

            for contact in container:
                ligIdxs = getPlipLigandIdxs(contact, idxField)
                if not ligIdxs:
                    continue

                resId = '{}_{}'.format(getPlipText(contact, 'restype'), getPlipText(contact, 'resnr'))
                pairId = (resId, '_'.join(ligIdxs))
                intDic.setdefault(pairId, {})[intType] = float(getPlipText(contact, distField, 0))

    return intDic

def getPlipText(element, tag, default=None):
    found = element.find(tag)
    return found.text if found is not None and found.text else default

def getPlipLigandIdxs(contact, idxField):
    """Atoms of the ligand taking part in an interaction: a list (salt bridges, rings), a single atom or,
    for the donor/acceptor pairs, whichever of the two is not the protein"""
    if isinstance(idxField, tuple):
        donField, accField = idxField
        idxField = accField if getPlipText(contact, 'protisdon', 'False') == 'True' else donField
        return [getPlipText(contact, idxField)] if getPlipText(contact, idxField) else []

    elif idxField.endswith('_list'):
        found = contact.find(idxField)
        return [idx.text for idx in found] if found is not None else []

    idx = getPlipText(contact, idxField)
    return [idx] if idx else []

def getPlipLigandName(xmlFile):
    """Residue name PLIP gave to the ligand of the complex"""
    hetId = ET.parse(xmlFile).getroot().find('.//identifiers/hetid')
    return hetId.text if hetId is not None else 'UNL'

def writePlipLigandFile(plipDir, xmlFile, outFile, complexFile):
    """Ligand alone, with the atom numbers the report refers to: those of the PDB PLIP fixed, or those
    of the complex it was given when it found nothing to fix and wrote no fixed copy"""
    ligName = getPlipLigandName(xmlFile)
    fixedFiles = [f for f in os.listdir(plipDir) if f.startswith('plipfixed.') and f.endswith('.pdb')]
    source = os.path.join(plipDir, fixedFiles[0]) if fixedFiles else complexFile

    with open(outFile, 'w') as fOut, open(source) as fIn:
        for line in fIn:
            #mgltools (sometimes) puts H in separate chains, so they would be drawn as loose atoms -> remove them
            if line.startswith(('ATOM', 'HETATM')) and line[17:20].strip() == ligName and line[76:78].strip() != 'H':
                fOut.write(line)
    return outFile

def writeDrawInteractionsParams(paramsFile, molFile, intDic, poseFile=None, valueUnit=None):
    """Parameters file read by drawInteractions.py. The pose file (with the bond orders a PDB has
    lost) is optional, but the molecule is drawn wrong without it"""
    with open(paramsFile, 'w') as f:
        f.write(f'molFile :: {os.path.abspath(molFile)}\n')
        f.write(f'intDic :: {intDic}\n')
        if poseFile:
            f.write(f'poseFile :: {os.path.abspath(poseFile)}\n')
        if valueUnit:
            f.write(f'valueUnit :: {valueUnit}\n')
    return paramsFile

def sortMolsByUnique(mols):
    if not mols:
        return []

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