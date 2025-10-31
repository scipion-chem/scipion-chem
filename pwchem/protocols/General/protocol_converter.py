# -*- coding: utf-8 -*-
# **************************************************************************
# *
# * Authors:  Alberto Manuel Parra Pérez (amparraperez@gmail.com)
# *           James Krieger (jamesmkrieger@gmail.com)
# *
# * Biocomputing Unit, CNB-CSIC
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

# General imports
import os, shutil, parmed, importlib, mdtraj

# Scipion chem imports
from pyworkflow.protocol.params import PointerParam, EnumParam, BooleanParam, LEVEL_ADVANCED
from pyworkflow.utils import Message
from pwem.objects.data import AtomStruct
from pwem.protocols import EMProtocol
from pwem.convert.atom_struct import toCIF, toPdb

# Plugin imports
from ... import Plugin
from ...objects import SetOfSmallMolecules, SmallMolecule, MDSystem
from ...utils import getBaseName
from ...constants import RDKIT_DIC, OPENBABEL_DIC, MDTRAJ_DIC

RDKIT, OBABEL = 'RDKit', 'OpenBabel'
extDic = {'PDB': '.pdb', 'cif': '.cif', 'Mol2': '.mol2', 'SDF': '.sdf', 'Smiles': '.smi'}

class ConvertStructures(EMProtocol):
    """
    Convert a set of input ligands or a protein structure to a specific file format
    """

    _label = 'Convert structure format'
    _program = ""

    def _defineParams(self, form):
        # Input type condition
        inputTypeCondition = 'isinstance(inputObject, '

        form.addSection(label=Message.LABEL_INPUT)
        form.addParam('inputObject', PointerParam, pointerClass="SetOfSmallMolecules,AtomStruct,MDSystem",
                      label='Input object:', allowsNull=False, important=True,
                      help="Input object. Can be of type Small molecules, Atomic structure, and MDSystem.\n\n"
                        "Small molecules: The allowed formats are pdb, cif, mol2, sdf and smi.\n"
                        "Atomic structure: The allowed format are pdb or cif.\n"
                        "MDSystem: Any system with a system file and trajectory file that mdtraj can read.")

        group = form.addGroup('Output', condition='inputObject')
        group.addParam('outputFormatSmall', EnumParam, default=2, condition=f'{inputTypeCondition}SetOfSmallMolecules)',
                       choices=['PDB', 'Mol2', 'SDF', 'Smiles'],
                       label='Output format: ',
                       help="Output format for the converted molecules")
        group.addParam('usePose', BooleanParam, default=False,
                       label='Use the docked ligands: ', expertLevel=LEVEL_ADVANCED,
                       help='Use the docked ligand files for preparation.')

        group.addParam('useManager', EnumParam, default=0, label='Convert using: ',
                      condition=f'{inputTypeCondition}SetOfSmallMolecules)', choices=[RDKIT, OBABEL],
                      help='Whether to convert the input molecules using RDKit or OpenBabel')


        group.addParam('outputFormatTarget', EnumParam, default=0,
                      condition=f'{inputTypeCondition}AtomStruct)', choices=['PDB', 'cif'], label='Output format')

        group.addParam('convSysFile', BooleanParam, default=False, label='Convert coordinates file: ',
                       condition=f'{inputTypeCondition}MDSystem)', help="Convert coordinates file from the MDSystem")
        group.addParam('outputSysFormat', EnumParam, default=0, label='System coordinates output format: ',
                       condition=f'{inputTypeCondition}MDSystem) and convSysFile', choices=['PDB', 'GRO'],
                       help="Output format for the coordinates of the system.")
        
        group.addParam('convTopFile', BooleanParam, default=True, label='Convert topology file: ',
                       condition=f'{inputTypeCondition}MDSystem)', help="Convert topology file from the MDSystem")
        group.addParam('outputTopFormat', EnumParam, default=0, label='Topology output format: ',
                       condition=f'{inputTypeCondition}MDSystem) and convTopFile', choices=['PSF', 'TOP', 'PRMTOP'],
                       help="Output format for the topology.")

        group.addParam('convTrjFile', BooleanParam, default=True, label='Convert trajectory file: ',
                       condition=f'{inputTypeCondition}MDSystem)', help="Convert coordinates file from the MDSystem")
        group.addParam('outputTrjFormat', EnumParam, default=0, label='Trajectory output format: ',
                       condition=f'{inputTypeCondition}MDSystem) and convTrjFile', choices=['DCD', 'GRO', 'NETCDF', 'PDB', 'TRR', 'XTC'],
                       help="Output format for the trajectory.")

    # --------------------------- Steps functions --------------------

    def _insertAllSteps(self):
        self._insertFunctionStep('convertStep')

    def convertStep(self):

        if isinstance(self.inputObject.get(), SetOfSmallMolecules):
            outputSmallMolecules = SetOfSmallMolecules().create(outputPath=self._getPath(), suffix='SmallMols')

            self.convErrors = []  # Save the file paths that could not be transformed
            for mol in self.inputObject.get():
                fnSmall = self.getMolFile(mol)
                fnRoot = os.path.splitext(os.path.split(fnSmall)[1])[0]

                outFormat = extDic[self.getEnumText('outputFormatSmall')]
                outDir = os.path.abspath(self._getExtraPath())
                fnOut = os.path.join(outDir, fnRoot + outFormat)

                args = ' -i "{}" -of {} -o {} --outputDir {}'.format(fnSmall, outFormat, fnOut, outDir)
                if self.getEnumText('useManager') == OBABEL:
                    Plugin.runScript(self, 'obabel_IO.py', args, env=OPENBABEL_DIC, cwd=outDir)
                else:
                    Plugin.runScript(self, 'rdkit_IO.py', args, env=RDKIT_DIC, cwd=outDir)

                if os.path.exists(fnOut):
                    smallMolecule = SmallMolecule(smallMolFilename=fnOut, molName='guess')
                    outputSmallMolecules.append(smallMolecule)
                else:
                    self.convErrors.append(fnRoot)

            if len(outputSmallMolecules) > 0:
                outputSmallMolecules.updateMolClass()
                self._defineOutputs(outputSmallMolecules=outputSmallMolecules)
                self._defineSourceRelation(self.inputObject, outputSmallMolecules)

            if len(self.convErrors) > 0:
                print("The following entries could not be converted: %s" % ','.join(self.convErrors))

        elif isinstance(self.inputObject.get(), AtomStruct):
            fnStructure = os.path.abspath(self.inputObject.get().getFileName())
            fnRoot = os.path.splitext(os.path.split(fnStructure)[1])[0]

            outFormat = extDic[self.getEnumText('outputFormatTarget')]
            outDir = os.path.abspath(self._getExtraPath())
            fnOut = os.path.join(outDir, fnRoot + outFormat)

            if outFormat == '.pdb':
                convFn = toPdb(fnStructure, fnOut)
            elif outFormat == '.cif':
                convFn = toCIF(fnStructure, fnOut)

            if convFn == fnStructure:
                shutil.copy(convFn, fnOut)

            if os.path.exists(fnOut):
                target = AtomStruct(filename=fnOut)
                self._defineOutputs(outputStructure=target)
                self._defineSourceRelation(self.inputObject, target)

        elif isinstance(self.inputObject.get(), MDSystem):
            inSystem = self.inputObject.get()
            sysFile = inSystem.getSystemFile()
            topFile = inSystem.getTopologyFile() # may be None
            topPath = os.path.split(topFile)[0]

            # Load system with ParmEd: topology + coordinates if available
            # Set GROMACS_TOPDIR if using a .top file
            if importlib.util.find_spec('gromacs'):
                from gromacs import Plugin as gromacsPlugin
                from gromacs.constants import GROMACS_DIC
                parmed.gromacs.GROMACS_TOPDIR = gromacsPlugin._getLocation(
                    GROMACS_DIC, marker='GROMACS_INSTALLED') + '/share/top'

            parmedSystem = parmed.load_file(topFile)
            parmedSystem.coordinates = parmed.load_file(sysFile).coordinates

            mdtrajTop = mdtraj.load(sysFile, top=sysFile).top
            resseqs = [r.resSeq for r in mdtrajTop.residues]
            assignChainsResseqs(parmedSystem, topPath, resseqs)

            if self.convSysFile.get():
                outDir = os.path.abspath(self._getExtraPath())
                fnRoot = os.path.splitext(os.path.split(sysFile)[1])[0]
                outFormat = self.getEnumText('outputSysFormat').lower()
                fnOut = os.path.join(outDir, fnRoot + '.' + outFormat)

                if outFormat == 'pdb':
                    parmedSystem.save(fnOut, renumber=False)
                else:
                    parmedSystem.save(fnOut)
                sysFile = fnOut

            outSystem = MDSystem(filename=sysFile)
            outSystem.setSystemFile(sysFile)
            
            if inSystem.hasTopology():
                topFile = inSystem.getTopologyFile()
                
                if self.convTopFile.get():
                    if topFile.endswith('.top') and importlib.util.find_spec('gromacs'):
                        from gromacs import Plugin as gromacsPlugin
                        from gromacs.constants import GROMACS_DIC
                        parmed.gromacs.GROMACS_TOPDIR = gromacsPlugin._getLocation(GROMACS_DIC, marker='GROMACS_INSTALLED') + '/share/top'

                    top = parmed.load_file(topFile)
                    topFile = self._getExtraPath('{}.{}'.format(getBaseName(topFile),
                                                                self.getEnumText('outputTopFormat').lower()))
                    top.save(topFile)
                outSystem.setTopologyFile(topFile)
            
            if inSystem.hasTrajectory():
                trjFile = inSystem.getTrajectoryFile()

                if self.convTrjFile.get():
                    outDir = os.path.abspath(self._getExtraPath())
                    fnRoot = os.path.splitext(os.path.split(sysFile)[1])[0]
                    outFormat = self.getEnumText('outputTrjFormat').lower()
                    fnOut = os.path.join(outDir, fnRoot + '.' + outFormat)

                    args = ' -s {} -o {} -t {}'.format(os.path.abspath(sysFile), fnOut, os.path.abspath(trjFile))
                    Plugin.runScript(self, 'mdtraj_IO.py', args, env=MDTRAJ_DIC, cwd=outDir)
                    trjFile = fnOut

                outSystem.setTrajectoryFile(trjFile)

            outSystem.setForceField(inSystem.getForceField())
            outSystem.setWaterForceField(inSystem.getWaterForceField())

            self._defineOutputs(outputSystem=outSystem)
            self._defineSourceRelation(self.inputObject, outSystem)

    def inputArg(self, fn):  # Input format file (fn)

        if fn.endswith('.pdb'):  # Protein Data Bank
            args = "-ipdb"
        elif fn.endswith('.cif'):  # cif (crystallography information)
            args = "-icif"

        elif fn.endswith('.sdf'):  # MDL MOL FORMAT
            args = "-isdf"
        elif fn.endswith('.sd'):  # MDL MOL FORMAT
            args = "-isd"

        elif fn.endswith('.mol2'):  # Sybyl Mol2 format (3D)
            args = "-imol2"
        elif fn.endswith('.smi') or fn.endswith('.smiles'):  # Smiles format (2D)
            args = "-ismi"
        else:
            error = " Input format was not recognize (The allowed format are pdb, cif, mol2, sdf and smi)"
            raise Exception(error)
        return args + " %s" % os.path.abspath(fn)

    def getMolFile(self, mol):
      if self.usePose.get():
        molFile = mol.getPoseFile()
      else:
        molFile = mol.getFileName()
      return os.path.abspath(molFile)

    # --------------------------- Summary functions --------------------
    def _summary(self):
        summary=[]

        if isinstance(self.inputObject.get(), SetOfSmallMolecules):
            if self.outputFormatSmall.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatSmall.get() == 2:
                summary.append('Converted to Mol2')
            elif self.outputFormatSmall.get() == 3:
                summary.append('Converted to SDF')
            elif self.outputFormatSmall.get() == 4:
                summary.append('Converted to Smiles')

        elif isinstance(self.inputObject.get(), AtomStruct):
            if self.outputFormatTarget.get() == 0:
                summary.append('Converted to PDB')
            elif self.outputFormatTarget.get() == 1:
                summary.append('Converted to Cif')
            elif self.outputFormatTarget.get() == 2:
                summary.append('Converted to Mol2')
        return summary

# --------------------------- Helper functions --------------------
def assignChainsResseqs(parmedSystem, topPath, resseqs):
    """
    Assign chain IDs to residues in a parmed.Structure object
    according to the number of residues in each .itp inclusion.
    """
    residues = parmedSystem.residues

    if hasattr(parmedSystem, 'itps') and parmedSystem.itps:
        proteinITPs = [itp for itp in parmedSystem.itps if "Protein" in itp]
        chainIDs = [fn.split("_")[-1].split(".")[0] for fn in proteinITPs]
        counts = [len(parmed.load_file(os.path.join(topPath, itp), parametrize=False).residues)
                for itp in proteinITPs]

        idx = 0
        for cid, numRes in zip(chainIDs, counts):
            chainResidues = residues[idx: idx + numRes]
            for res in chainResidues:
                res.chain = cid
            idx += numRes

    for i, res in enumerate(residues):
        res.number = resseqs[i]

    missing_dihedrals = [dih for dih in parmedSystem.dihedrals if dih.funct is None]
    for dih in missing_dihedrals:
        dih.funct = 1

    missing_angles = [(i, a) for (i, a) in enumerate(parmedSystem.angles) if a.type is None]
    for i, angle in missing_angles:
        if angle.type is None:
            angle_type = parmed.AngleType(1.0, 109.5)
            # Replace the angle in the tracked list to ensure it's linked properly
            parmedSystem.angles[i].type = angle_type
        # sanity check
        if angle.type is None or not hasattr(angle.type, 'theteq'):
            print(f"Warning: angle {angle} still has no type")
