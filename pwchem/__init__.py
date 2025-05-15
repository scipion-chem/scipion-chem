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

# General imports
import os, subprocess
from subprocess import run

# Scipion core imports
from pyworkflow import VarTypes
import pyworkflow.utils as pwutils
from pyworkflow.tests import DataSet
from scipion.install.funcs import InstallHelper
import pwem

# Plugin imports
from .bibtex import _bibtexStr
from .constants import *
from . import version

_logo = 'pwchem_logo.png'
__version__ = version.__version__

class Plugin(pwem.Plugin):
	@classmethod
	def defineBinaries(cls, env):
		cls.addRDKitPackage(env)
		cls.addOpenbabelPackage(env)
		cls.addMGLToolsPackage(env)
		cls.addJChemPaintPackage(env)
		cls.addAliViewPackage(env)
		cls.addVMDPackage(env)
		cls.addMDTrajPackage(env)
		cls.addDEAPPackage(env)
		cls.addRanxPackage(env)

	@classmethod
	def _defineVariables(cls):
		# Package home directories
		cls._defineEmVar(RDKIT_DIC['home'], cls.getEnvName(RDKIT_DIC))
		cls._defineEmVar(MGL_DIC['home'], cls.getEnvName(MGL_DIC))
		cls._defineEmVar(JCHEM_DIC['home'], cls.getEnvName(JCHEM_DIC))
		cls._defineEmVar(ALIVIEW_DIC['home'], cls.getEnvName(ALIVIEW_DIC))
		cls._defineEmVar(VMD_DIC['home'], cls.getEnvName(VMD_DIC))
		cls._defineEmVar(OPENBABEL_DIC['home'], cls.getEnvName(OPENBABEL_DIC))

		# Common enviroments
		cls._defineVar('RDKIT_ENV_ACTIVATION', cls.getEnvActivationCommand(RDKIT_DIC))
		cls._defineVar('BIOCONDA_ENV_ACTIVATION', cls.getEnvActivationCommand(BIOCONDA_DIC))
		cls._defineVar('OPENABEL_ENV_ACTIVATION', cls.getEnvActivationCommand(OPENBABEL_DIC))
		cls._defineVar(MAX_MOLS_SET, 1000000, var_type=VarTypes.INTEGER,
									 description='Maximum size for a SetOfSmallMolecules with 1 file per molecule to avoid memory '
															 'and IO overuse')

########################### ENVIROMENT MANIPULATION COMMON FUNCTIONS ###########################
	@classmethod
	def getEnvName(cls, packageDictionary):
		""" This function returns the name of the conda enviroment for a given package. """
		return '{}-{}'.format(packageDictionary['name'], packageDictionary['version'])

	@classmethod
	def getEnvActivationCommand(cls, packageDictionary, condaHook=True):
		""" This function returns the conda enviroment activation command for a given package. """
		return '{}conda activate {}'.format(cls.getCondaActivationCmd() if condaHook else '', cls.getEnvName(packageDictionary))
	
	@classmethod
	def getEnvPath(cls, packageDictionary=None, innerPath=None):
		""" This function returns the root path for the env given it's dictionary. """
		# Command to extract the base path
		condaCmd = f'{cls.getCondaActivationCmd()} conda activate && echo $CONDA_PREFIX'

		# Getting the base path
		basePath = subprocess.run(condaCmd, shell=True, stdout=subprocess.PIPE).stdout.strip().decode('utf-8')

		# Getting env base path
		envBasePath = os.path.join(basePath, 'envs', cls.getEnvName(packageDictionary)) if packageDictionary else basePath

		# Returning env path
		return os.path.join(envBasePath, innerPath) if innerPath else envBasePath

######################## PACKAGES #########################
	@classmethod
	def addRDKitPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(RDKIT_DIC['name'], packageHome=cls.getVar(RDKIT_DIC['home']), packageVersion=RDKIT_DIC['version'])

		# Defining env path
		env_path = os.environ.get('PATH', "")  # keep path since conda likely in there

		# Installing package
		rdkitEnvName = cls.getEnvName(RDKIT_DIC)
		installer.addCommand(f'conda create -c conda-forge --name {rdkitEnvName} '
												 f'{RDKIT_DIC["name"]}={RDKIT_DIC["version"]} oddt=0.7 python=3.10 -y', 'RDKIT_ENV_CREATED')\
			.addCommand('mkdir oddtModels', 'ODTMODELS_CREATED')\
			.addPackage(env, dependencies=['conda'], default=default, vars={'PATH': env_path} if env_path else None)
			
	@classmethod
	def addMGLToolsPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(MGL_DIC['name'], packageHome=cls.getVar(MGL_DIC['home']), packageVersion=MGL_DIC['version'])

		# Defining file names
		tar_file = cls.getDefTar(MGL_DIC)

		# Installing package
		installer.getExtraFile('https://ccsb.scripps.edu/download/532/', 'MGLTOOLS_DOWNLOADED', fileName=tar_file)\
			.addCommand(f'tar -xf {tar_file} --strip-components 1 && rm {tar_file}', 'MGLTOOLS_EXTRACTED')\
			.addCommand('export DISPLAY= && {}'.format(cls.getDefPath(MGL_DIC, 'install.sh')), 'MGLTOOLS_INSTALLED')\
			.addPackage(env, dependencies=['wget', 'tar'], default=default)

	@classmethod
	def addJChemPaintPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(JCHEM_DIC['name'], packageHome=cls.getVar(JCHEM_DIC['home']), packageVersion=JCHEM_DIC['version'])

		# Defining filename to download
		jchem_file = f"jchempaint-{JCHEM_DIC['version']}.jar"

		installer.getExtraFile(f"https://sourceforge.net/projects/cdk/files/JChemPaint/{JCHEM_DIC['version']}/{jchem_file}/download", 'JCHEM_DOWNLOADED', fileName=jchem_file)\
			.addCommand(f'chmod +x {jchem_file}', 'JCHEM_INSTALLED')\
			.addPackage(env, dependencies=['wget'], default=default)

	@classmethod
	def addOpenbabelPackage(cls, env, default=True):
		# Instantiating openbabel install helper
		openbabel_installer = InstallHelper(OPENBABEL_DIC['name'], packageHome=cls.getVar(SHAPEIT_DIC['home']), packageVersion=OPENBABEL_DIC['version'])

		# Generating installation commands
		openbabel_installer.getCondaEnvCommand()\
			.addCondaPackages(['openbabel', 'swig', 'plip', 'pdbfixer', 'pymol-open-source'], channel='conda-forge') \
			.addCondaPackages(['clustalo'], channel='bioconda', targetName='CLUSTALO_INSTALLED')\
			.addPackage(env, dependencies=['git', 'conda', 'cmake', 'make'], default=default)
		
		# # Instantiating shape it install helper
		# shape_it_installer = InstallHelper(SHAPEIT_DIC['name'], packageHome=cls.getVar(SHAPEIT_DIC['home']), packageVersion=SHAPEIT_DIC['version'])
		#
		# # Importing commands from openbabel and rdkit installers
		# shape_it_installer.importCommandList(openbabel_installer.getCommandList())
		#
		# # Defining binaries folder name
		# binaries_directory = SHAPEIT_DIC['name']
		#
		# # Installing package
		# shape_it_installer.getCloneCommand('https://github.com/rdkit/shape-it.git', binaryFolderName=binaries_directory)\
		# 	.addCommand(f'{cls.getEnvActivationCommand(RDKIT_DIC)} && cmake -DCMAKE_INSTALL_PREFIX=. -DOPENBABEL3_INCLUDE_DIR=$CONDA_PREFIX/include/openbabel3 -DOPENBABEL3_LIBRARIES=$CONDA_PREFIX/lib/libopenbabel.so -Bbuild .', 'MAKEFILES_BUILT', workDir=binaries_directory)\
		# 	.addCommand(f'cd {binaries_directory}/build && make', 'SHAPEIT_COMPILED')\
		# 	.addPackage(env, dependencies=['git', 'conda', 'cmake', 'make'], default=default)

	@classmethod
	def addAliViewPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(ALIVIEW_DIC['name'], packageHome=cls.getVar(ALIVIEW_DIC['home']), packageVersion=ALIVIEW_DIC['version'])

		# Defining filename
		file_name = cls.getDefTar(ALIVIEW_DIC)

		# Installing package
		installer.getExtraFile(cls.getAliviewUrl(), 'ALIVIEW_DOWNLOADED', fileName=file_name)\
			.addCommand(f'tar -xf {file_name} && rm {file_name}', 'ALIVIEW_EXTRACTED')\
			.addCommand(f"conda create --name {BIOCONDA_DIC['name']}-{BIOCONDA_DIC['version']} --file {cls.getEnvSpecsPath('bioconda')} -y", 'BIOCONDA_ENV_CREATED')\
			.addPackage(env, dependencies=['wget', 'conda'], default=default)

	@classmethod
	def addVMDPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(VMD_DIC['name'], packageHome=cls.getVar(VMD_DIC['home']), packageVersion=VMD_DIC['version'])

		installer.getCondaEnvCommand().addCondaPackages(['vmd'], channel='conda-forge')\
			.addPackage(env, dependencies=['conda'], default=default)

	@classmethod
	def addMDTrajPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(MDTRAJ_DIC['name'], packageHome=cls.getVar(MDTRAJ_DIC['home']), packageVersion=MDTRAJ_DIC['version'])

		installer.getCondaEnvCommand().addCondaPackages(['mdtraj', 'matplotlib', 'acpype'], channel='conda-forge')\
			.addPackage(env, dependencies=['conda'], default=default)

	@classmethod
	def addDEAPPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(DEAP_DIC['name'], packageHome=cls.getVar(DEAP_DIC['home']),
															packageVersion=DEAP_DIC['version'])

		scipionEnvPath = cls.getEnvPath(innerPath='envs/scipion3/lib/python3.8/site-packages/grape')
		installer.addCommand(f'{cls.getCondaActivationCmd()}conda activate scipion3 && conda install conda-forge::deap -y')\
			.addCommand('git clone https://github.com/bdsul/grape.git') \
			.addCommand(f'mv grape {scipionEnvPath}') \
			.addPackage(env, dependencies=['conda', 'git'], default=default)

	@classmethod
	def addRanxPackage(cls, env, default=True):
		# Instantiating install helper
		installer = InstallHelper(RANX_DIC['name'], packageHome=cls.getVar(RANX_DIC['home']),
															packageVersion=RANX_DIC['version'])

		installer.getCondaEnvCommand(RANX_DIC['name'], binaryVersion=RANX_DIC['version'], pythonVersion='3.10').\
			addCommand(f'{cls.getEnvActivationCommand(RANX_DIC)} && pip install ranx', 'RANKX_INSTALLED') \
			.addPackage(env, dependencies=['conda', 'pip'], default=default)

	##################### RUN CALLS ######################
	@classmethod
	def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False, wait=True, scriptDir=None, pyStr='python'):
		""" Run a script from a given protocol using a specific environment """
		scriptName = cls.getScriptsDir(scriptName) if scriptDir == None else os.path.join(scriptDir, scriptName)
		fullProgram = '%s && %s %s' % (cls.getEnvActivationCommand(env), pyStr, scriptName)

		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)
		else:
			if wait:
				subprocess.check_call(f'{fullProgram} {args}', cwd=cwd, shell=True)
			else:
				subprocess.Popen(f'{fullProgram} {args}', cwd=cwd, shell=True)

	@classmethod
	def runCondaCommand(cls, protocol, args, condaDic, program, cwd=None, popen=False, silent=True):
		""" General function to run conda commands """
		fullProgram = f'{cls.getEnvActivationCommand(condaDic)} && {program} '
		if not popen:
			protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads=1)
		else:
			kwargs = {}
			if silent:
				kwargs = {"stdout": subprocess.DEVNULL, "stderr": subprocess.DEVNULL}
			run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True, **kwargs)

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
	def runOPENBABEL(cls, protocol, program="obabel ", args=None, cwd=None, popen=False, silent=True):
		""" Run openbabel command from a given protocol. """
		cls.runCondaCommand(protocol, args, OPENBABEL_DIC, program, cwd, popen, silent)

	@classmethod
	def runACPYPE(cls, protocol, program='acpype', args=None, cwd=None, popen=False, silent=True):
		""" Run ACPYPE command from a given protocol. """
		cls.runCondaCommand(protocol, args, MDTRAJ_DIC, program, cwd, popen, silent)

	@classmethod
	def runPLIP(cls, args, cwd=None):
		""" Run PLIP command from a given protocol. """
		fullProgram = '%s && %s ' % (cls.getEnvActivationCommand(OPENBABEL_DIC), 'plip')
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
		oddtModelsDir = os.path.abspath(os.path.join(cls.getVar(RDKIT_DIC['home']), 'oddtModels'))
		if not os.path.exists(oddtModelsDir):
			os.mkdir(oddtModelsDir)
		return os.path.join(oddtModelsDir, path)

	@classmethod
	def getDefTar(cls, programDic, ext='tgz'):
		return os.path.join(cls.getDefPath(programDic), '{}-{}.{}'.format(programDic['name'], programDic['version'], ext))

	@classmethod
	def getAliviewUrl(cls, version='1.28'):
		return f'http://www.ormbunkar.se/aliview/downloads/linux/linux-versions-all/linux-version-{version}/aliview.tgz'

DataSet(name='smallMolecules', folder='smallMolecules',
					files={
						'mix': 'mix/',
						'mol2': 'mol2/',
						'pdb': 'pdb/',
						'sdf': 'sdf/',
						'smi': 'smi/'})
