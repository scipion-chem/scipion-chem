# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
# *              Daniel Del Hoyo Gomez (ddelhoyo@cnb.csic.es)
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

"""
This package contains the protocols for
manipulation of atomic struct objects
"""

import os, subprocess
from subprocess import run
import pyworkflow.utils as pwutils
from pyworkflow.tests import DataSet
import pwem
from .bibtex import _bibtexStr
from .constants import *
from pwchem.install_helper import InstallHelper

_logo = 'pwchem_logo.png'
__version__ = DEFAULT_VERSION

class Plugin(pwem.Plugin):
	@classmethod
	def defineBinaries(cls, env):
		cls.addRDKitPackage(env)
		cls.addShapeitPackage(env)
		cls.addMGLToolsPackage(env)
		cls.addJChemPaintPackage(env)
		cls.addPyMolPackage(env)
		cls.addAliViewPackage(env)
		cls.addVMDPackage(env)

	@classmethod
	def _defineVariables(cls):
		cls._defineEmVar(MGL_DIC['home'], '{}-{}'.format(MGL_DIC['name'], MGL_DIC['version']))
		cls._defineEmVar(PYMOL_DIC['home'], '{}-{}'.format(PYMOL_DIC['name'], PYMOL_DIC['version']))
		cls._defineEmVar(JCHEM_DIC['home'], '{}-{}'.format(JCHEM_DIC['name'], JCHEM_DIC['version']))
		cls._defineEmVar(ALIVIEW_DIC['home'], '{}-{}'.format(ALIVIEW_DIC['name'], ALIVIEW_DIC['version']))
		cls._defineEmVar(VMD_DIC['home'], '{}-{}'.format(VMD_DIC['name'], VMD_DIC['version']))
		cls._defineEmVar(SHAPEIT_DIC['home'], '{}-{}'.format(SHAPEIT_DIC['name'], SHAPEIT_DIC['version']))

	@classmethod
	def getEnvActivationCommand(cls, packageDictionary, condaHook=True):
		""" This function returns the conda enviroment activation command for a given package. """
		return '{}conda activate {}-{}'.format(cls.getCondaActivationCmd() if condaHook else '', packageDictionary['name'], packageDictionary['version'])

######################## PACKAGES #########################
	@classmethod
	def addRDKitPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(RDKIT_DIC['name'], packageHome=cls.getVar(RDKIT_DIC['home']), packageVersion=RDKIT_DIC['version'])

		# Defining env path
		envPath = os.environ.get('PATH', "")  # keep path since conda likely in there

		# Installing package
		installer.addCommand(f'conda create --name {RDKIT_DIC["name"]}-{RDKIT_DIC["version"]} --file {cls.getEnvSpecsPath("rdkit")} -y', 'RDKIT_ENV_CREATED')\
			.addCommand('mkdir oddtModels', 'ODTMODELS_CREATED')\
			.addPackage(env, dependencies=['conda'], default=default, vars={'PATH': envPath} if envPath else None)
			
	@classmethod
	def addPyMolPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(PYMOL_DIC['name'], packageHome=cls.getVar(PYMOL_DIC['home']), packageVersion=PYMOL_DIC['version'])

		# Installing package
		installer.getExtraFile('https://pymol.org/installers/PyMOL-' + PYMOL_DIC['version'] + '_496-Linux-x86_64-py37.tar.bz2', 'PYMOL_DOWNLOADED')\
			.addCommand('tar -jxf PyMOL-' + PYMOL_DIC['version'] + '_496-Linux-x86_64-py37.tar.bz2', 'PYMOL_EXTRACTED')\
			.addCommand('rm PyMOL-2.5.5_496-Linux-x86_64-py37.tar.bz2', 'TAR_REMOVED')\
			.addPackage(env, dependencies=['tar', 'wget'], default=default)

	@classmethod
	def addMGLToolsPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(MGL_DIC['name'], packageHome=cls.getVar(MGL_DIC['home']), packageVersion=MGL_DIC['version'])

		# Defining file names
		tarFile = cls.getDefTar(MGL_DIC)

		# Installing package
		installer.getExtraFile('https://ccsb.scripps.edu/download/532/', 'MGLTOOLS_DOWNLOADED', fileName=tarFile)\
			.addCommand(f'tar -xf {tarFile} --strip-components 1 && rm {tarFile}', 'MGLTOOLS_EXTRACTED')\
			.addCommand(cls.getDefPath(MGL_DIC, 'install.sh'), 'MGLTOOLS_INSTALLED')\
			.addPackage(env, dependencies=['wget', 'tar'], default=default)

	@classmethod
	def addJChemPaintPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(JCHEM_DIC['name'], packageHome=cls.getVar(JCHEM_DIC['home']), packageVersion=JCHEM_DIC['version'])

		# Defining filename to download
		jchemFile = f"jchempaint-{JCHEM_DIC['version']}.jar"

		installer.getExtraFile(f"https://sourceforge.net/projects/cdk/files/JChemPaint/{JCHEM_DIC['version']}/{jchemFile}/download", 'JCHEM_DOWNLOADED', fileName=jchemFile)\
			.addCommand(f'chmod +x {jchemFile}', 'JCHEM_INSTALLED')\
			.addPackage(env, dependencies=['wget'], default=default)

	@classmethod
	def addShapeitPackage(cls, env, default=True):
		# Instantiating plip install helper
		plipInstaller = InstallHelper(PLIP_DIC['name'], packageHome=cls.getVar(SHAPEIT_DIC['home']), packageVersion=PLIP_DIC['version'])

		# Generating installation commands
		plipInstaller.getCondaEnvCommand(requirementsFile=False)\
			.addCondaPackages(['openbabel', 'swig', 'plip'], channel='conda-forge')\
			.addCondaPackages(['clustalo'], channel='bioconda', targetName='CLUSTALO_INSTALLED')
		
		# Instantiating shape it install helper
		shapeItInstaller = InstallHelper(SHAPEIT_DIC['name'], packageHome=cls.getVar(SHAPEIT_DIC['home']), packageVersion=SHAPEIT_DIC['version'])

		# Importing commands from plip and rdkit installers
		shapeItInstaller.importCommandList(plipInstaller.getCommandList())

		# Defining binaries folder name
		binariesDirectory = SHAPEIT_DIC['name']

		# Installing package
		shapeItInstaller.getCloneCommand('https://github.com/rdkit/shape-it.git', binaryFolderName=binariesDirectory)\
			.addCommand(f'{cls.getEnvActivationCommand(RDKIT_DIC, condaHook=True)} && cmake -DCMAKE_INSTALL_PREFIX=. -DOPENBABEL3_INCLUDE_DIR=$CONDA_PREFIX/include/openbabel3 -DOPENBABEL3_LIBRARIES=$CONDA_PREFIX/lib/libopenbabel.so -Bbuild .', 'MAKEFILES_BUILT', workDir=binariesDirectory)\
			.addCommand(f'cd {binariesDirectory}/build && make', 'SHAPEIT_COMPILED')\
			.addPackage(env, dependencies=['git', 'conda', 'cmake', 'make'], default=default)

	@classmethod
	def addAliViewPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(ALIVIEW_DIC['name'], packageHome=cls.getVar(ALIVIEW_DIC['home']), packageVersion=ALIVIEW_DIC['version'])

		# Defining filename
		fileName = cls.getDefTar(ALIVIEW_DIC)

		# Installing package
		installer.getExtraFile('https://ormbunkar.se/aliview/downloads/linux/linux-version-1.28/aliview.tgz', 'ALIVIEW_DOWNLOADED', fileName=fileName)\
			.addCommand(f'tar -xf {fileName} && rm {fileName}', 'ALIVIEW_EXTRACTED')\
			.addCommand(f"conda create --name {BIOCONDA_DIC['name']}-{BIOCONDA_DIC['version']} --file {cls.getEnvSpecsPath('bioconda')} -y", 'BIOCONDA_ENV_CREATED')\
			.addPackage(env, dependencies=['wget', 'conda'], default=default)

	@classmethod
	def addVMDPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(VMD_DIC['name'], packageHome=cls.getVar(VMD_DIC['home']), packageVersion=VMD_DIC['version'])

		installer.getCondaEnvCommand(requirementsFile=False).addCondaPackages(['vdm'], channel='conda-forge')\
			.addPackage(env, dependencies=['conda'], default=default)

	##################### RUN CALLS ######################
	@classmethod
	def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False):
		""" Run rdkit command from a given protocol. """
		scriptName = cls.getScriptsDir(scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(env), 'python', scriptName)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			subprocess.check_call(fullProgram + args, cwd=cwd, shell=True)

	@classmethod
	def runShapeIt(cls, protocol, program, args, cwd=None):
		""" Run shapeit command from a given protocol (it must be run from the shape-it dir, where the lib file is) """
		progDir = cls.getProgramHome(SHAPEIT_DIC)
		progFile = os.path.join(os.path.abspath(progDir), args.split('-s ')[1].split()[0])
		outFile = os.path.join(os.path.abspath(cwd), os.path.basename(progFile))
		protocol.runJob(program, args, env=cls.getEnviron(), cwd=progDir)
		os.rename(progFile, outFile)

	@classmethod
	def runJChemPaint(cls, protocol, cwd=None):
		""" Run jchempaint command from a given protocol. """
		protocol.runJob('java -jar {}'.format(cls.getProgramHome(JCHEM_DIC, 'jchempaint-{}.jar'.format(JCHEM_DIC['version']))), arguments='', env=cls.getEnviron(), cwd=cwd)

	@classmethod
	def runOPENBABEL(cls, protocol, program="obabel ", args=None, cwd=None, popen=False):
		""" Run openbabel command from a given protocol. """
		fullProgram = '%s && %s' % (cls.getEnvActivationCommand(PLIP_DIC, condaHook=True), program)
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads=1)
		else:
			run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True)

	@classmethod
	def runPLIP(cls, args, cwd=None):
		""" Run PLIP command from a given protocol. """
		fullProgram = '%s && %s ' % (cls.getCondaActivationCmd(), cls.getEnvActivationCommand(PLIP_DIC, condaHook=True), 'plip')
		run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True)

  ##################### UTILS ###########################
	@classmethod
	def getPluginHome(cls, path=""):
		import pwchem
		fnDir = os.path.split(pwchem.__file__)[0]
		return os.path.join(fnDir, path)

	# Default version paths
	@classmethod
	def getDefPath(cls, programDic, path=''):
		return os.path.join(pwem.Config.EM_ROOT, '{}-{}'.format(programDic['name'], programDic['version']), path)

	# In use paths
	@classmethod
	def getProgramHome(cls, programDic, path=''):
		return os.path.join(cls.getVar(programDic['home']), path)

	@classmethod
	def getScriptsDir(cls, scriptName):
		return cls.getPluginHome('scripts/%s' % scriptName)

	@classmethod
	def getEnvSpecsPath(cls, env):
		envFile = '/{}_env_spec.txt'.format(env) if env else ''
		return cls.getPluginHome('envs%s' % envFile)

	@classmethod
	def getMGLEnviron(cls):
		""" Create the needed environment for MGL Tools programs. """
		environ = pwutils.Environ(os.environ)
		pos = pwutils.Environ.BEGIN
		environ.update({
			'PATH': cls.getProgramHome(programDic=MGL_DIC, path='bin')
		}, position=pos)
		return environ

	@classmethod
	def getODDTModelsPath(cls, path=''):
		return os.path.abspath(os.path.join(cls.getVar(RDKIT_DIC['home']), 'oddtModels', path))

	@classmethod
	def getDefTar(cls, programDic, ext='tgz'):
		return os.path.join(cls.getDefPath(programDic), '{}-{}.{}'.format(programDic['name'], programDic['version'], ext))

DataSet(name='smallMolecules', folder='smallMolecules',
					files={
						'mix': 'mix/',
						'mol2': 'mol2/',
						'pdb': 'pdb/',
						'sdf': 'sdf/',
						'smi': 'smi/'})