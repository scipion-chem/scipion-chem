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

_logo = 'pwchem_logo.png'
RDKIT = 'rdkit'

class Plugin(pwem.Plugin):
    _rdkitHome = os.path.join(pwem.Config.EM_ROOT, RDKIT)

    @classmethod
    def defineBinaries(cls, env):
        #PLIP environment (with pymol bundle)
        cls.addPLIPPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addRDKitPackage(env, default=bool(cls.getCondaActivationCmd()))
        # cls.addShapeitPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addMGLToolsPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addJChemPaintPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addPyMolPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addAliViewPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addVMDPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def _defineVariables(cls):
        cls._defineVar("RDKIT_ENV_ACTIVATION", 'conda activate rdkit-env')
        cls._defineVar("PLIP_ENV_ACTIVATION", 'conda activate plip-env')
        cls._defineVar("VMD_ENV_ACTIVATION", 'conda activate vmd-env')
        cls._defineEmVar(MGL_DIC['home'], '{}-{}'.format(MGL_DIC['name'], MGL_DIC['version']))
        cls._defineEmVar(PYMOL_DIC['home'], '{}-{}'.format(PYMOL_DIC['name'], PYMOL_DIC['version']))
        cls._defineEmVar(JCHEM_DIC['home'], '{}-{}'.format(JCHEM_DIC['name'], JCHEM_DIC['version']))
        cls._defineEmVar(ALIVIEW_DIC['home'], '{}-{}'.format(ALIVIEW_DIC['name'], ALIVIEW_DIC['version']))
        # cls._defineEmVar(SHAPEIT_DIC['home'], '{}-{}'.format(SHAPEIT_DIC['name'], SHAPEIT_DIC['version']))

    @classmethod
    def getEnvActivation(cls, env):
      activation = cls.getVar("{}_ENV_ACTIVATION".format(env.upper()))
      return activation

    @classmethod
    def getRDKitEnvActivation(cls):
      activation = cls.getVar("RDKIT_ENV_ACTIVATION")
      return activation

    @classmethod
    def getPLIPEnvActivation(cls):
      activation = cls.getVar("PLIP_ENV_ACTIVATION")
      return activation


######################## PACKAGES #########################
    @classmethod
    def addPLIPPackage(cls, env, default=False):
      PLIP_INSTALLED = 'plip_installed'

      installationCmd = cls.getCondaActivationCmd()
      installationCmd += 'conda create --name plip-env --file {} && '.format(cls.getEnvSpecsPath('plip'))
      installationCmd += 'touch {}'.format(PLIP_INSTALLED)

      installationCmd = [(installationCmd, PLIP_INSTALLED)]
      env.addPackage(PLIP_DIC['name'], version=PLIP_DIC['version'],
                     tar='void.tgz',
                     commands=installationCmd,
                     default=True)

    @classmethod
    def addPyMolPackage(cls, env, default=False):
      PYMOL_INSTALLED = 'pymol_installed'
      pymol_commands = 'wget https://pymol.org/installers/PyMOL-2.5.2_293-Linux-x86_64-py37.tar.bz2 -O {} && '.\
        format(cls.getDefTar(PYMOL_DIC, ext='tar.bz2'))
      pymol_commands += 'tar -jxf {} --strip-components 1 && rm {} &&'.\
        format(*[cls.getDefTar(PYMOL_DIC, ext='tar.bz2')]*2)
      pymol_commands += 'touch ' + PYMOL_INSTALLED
      pymol_commands = [(pymol_commands, PYMOL_INSTALLED)]
      env.addPackage(PYMOL_DIC['name'], version=PYMOL_DIC['version'],
                     tar='void.tgz',
                     commands=pymol_commands,
                     default=True)

    @classmethod
    def addMGLToolsPackage(cls, env, default=False):
        MGL_INSTALLED = "initMGLtools.sh"
        mgl_commands = 'wget https://ccsb.scripps.edu/download/548/ -O {} --no-check-certificate && '. \
            format(cls.getDefTar(MGL_DIC))
        mgl_commands += 'tar -xf {} --strip-components 1 && rm {} &&'.format(*[cls.getDefTar(MGL_DIC)] * 2)
        mgl_commands += 'cp install.sh install.bash && sed -i "s/bin\/sh/bin\/bash/g" install.bash && '
        mgl_commands += '{} && '.format(cls.getDefPath(MGL_DIC, 'install.bash'))
        mgl_commands += 'touch ' + MGL_INSTALLED
        mgl_commands = [(mgl_commands, MGL_INSTALLED)]

        env.addPackage(MGL_DIC['name'], version=MGL_DIC['version'],
                       tar='void.tgz',
                       commands=mgl_commands,
                       default=True)

    @classmethod
    def addJChemPaintPackage(cls, env, default=False):
        JCHEM_INSTALLED = 'jchem_installed'
        jchem_commands = 'wget https://github.com/downloads/JChemPaint/jchempaint/jchempaint-3.3-1210.jar -O {} && '.\
          format(cls.getDefPath(JCHEM_DIC, 'jchempaint-{}.jar'.format(JCHEM_DIC['version'])))
        jchem_commands += 'chmod +x {} && '.format(cls.getDefPath(JCHEM_DIC, 'jchempaint-{}.jar'.
                                                                  format(JCHEM_DIC['version'])))
        jchem_commands += ' touch %s' % JCHEM_INSTALLED

        jchem_commands = [(jchem_commands, JCHEM_INSTALLED)]

        env.addPackage(JCHEM_DIC['name'], version=JCHEM_DIC['version'],
                       tar='void.tgz',
                       commands=jchem_commands,
                       default=True)

    @classmethod
    def addRDKitPackage(cls, env, default=False):
        RDKIT_INSTALLED = 'rdkit_installed'

        installationCmd = cls.getCondaActivationCmd()
        installationCmd += 'conda create --name rdkit-env --file {} && '.format(cls.getEnvSpecsPath('rdkit'))
        installationCmd += 'mkdir oddtModels && touch %s' % RDKIT_INSTALLED

        rdkit_commands = [(installationCmd, RDKIT_INSTALLED)]

        envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage(RDKIT,
                       tar='void.tgz',
                       commands=rdkit_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def addShapeitPackage(cls, env, default=False):
      SHAPEIT_INSTALLED = 'shapeit_installed'

      installationCmd = ' git clone https://github.com/rdkit/shape-it.git && cd shape-it && mkdir build && cd build &&'
      installationCmd += ' cmake -DCMAKE_INSTALL_PREFIX=. .. && make && make install && cd ../.. &&'
      installationCmd += ' mv shape-it/* . && rm -rf shape-it && mv build/shape-it . && mv build/libshapeit_lib.so . &&'
      installationCmd += ' touch %s' % SHAPEIT_INSTALLED
      installationCmd = [(installationCmd, SHAPEIT_INSTALLED)]

      env.addPackage(SHAPEIT_DIC['name'], version=SHAPEIT_DIC['version'],
                     tar='void.tgz',
                     commands=installationCmd,
                     neededProgs=cls.getDependencies(),
                     default=default)

    @classmethod
    def addAliViewPackage(cls, env, default=False):
      SEQS_INSTALLED = 'aliview_installed'
      seqs_commands = 'wget https://ormbunkar.se/aliview/downloads/linux/linux-version-1.28/aliview.tgz -O {} ' \
                      '--no-check-certificate && '.format(cls.getDefTar(ALIVIEW_DIC))
      seqs_commands += 'tar -xf {} && rm {} &&'.format(*[cls.getDefTar(ALIVIEW_DIC)]*2)
      seqs_commands += ' conda install -y -c bioconda clustalo &&'
      seqs_commands += ' conda install -y -c bioconda muscle &&'
      seqs_commands += ' conda install -y -c bioconda mafft &&'
      seqs_commands += ' conda install -y -c bioconda emboss &&'
      seqs_commands += ' touch %s' % SEQS_INSTALLED

      seqs_commands = [(seqs_commands, SEQS_INSTALLED)]

      env.addPackage(ALIVIEW_DIC['name'], version=ALIVIEW_DIC['version'],
                     tar='void.tgz',
                     commands=seqs_commands,
                     default=True)

    @classmethod
    def addVMDPackage(cls, env, default=False):
      VMD_INSTALLED = 'vmd_installed'

      installationCmd = cls.getCondaActivationCmd()
      installationCmd += ' conda create -y -c conda-forge -n vmd-env vmd &&'
      installationCmd += ' touch %s' % VMD_INSTALLED
      vmd_commands = [(installationCmd, VMD_INSTALLED)]

      env.addPackage(VMD_DIC['name'], version=VMD_DIC['version'],
                     tar='void.tgz',
                     commands=vmd_commands,
                     neededProgs=cls.getDependencies(),
                     default=default)


    ##################### RUN CALLS #################33333
    @classmethod
    def runScript(cls, protocol, scriptName, args, env, cwd=None, popen=False):
      """ Run rdkit command from a given protocol. """
      scriptName = cls.getScriptsDir(scriptName)
      fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getEnvActivation(env), 'python', scriptName)
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
    def runRDKit(cls, protocol, program, args, cwd=None):
      """ Run rdkit command from a given protocol. """
      fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getRDKitEnvActivation(), program)
      protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runJChemPaint(cls, protocol, cwd=None):
      """ Run jchempaint command from a given protocol. """
      protocol.runJob('java -jar {}'.format(cls.getProgramHome(JCHEM_DIC, 'jchempaint-{}.jar'.
                                                               format(JCHEM_DIC['version']))),
                      arguments='', env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runOPENBABEL(cls, protocol, program="obabel ", args=None, cwd=None, popen=False):
      """ Run openbabel command from a given protocol. """
      fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getPLIPEnvActivation(), program)
      if not popen:
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads=1)
      else:
        run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True)

    @classmethod
    def runPLIP(cls, args, cwd=None):
        """ Run PLIP command from a given protocol. """
        fullProgram = '%s %s && %s ' % (cls.getCondaActivationCmd(), cls.getPLIPEnvActivation(), 'plip')
        run(fullProgram + args, env=cls.getEnviron(), cwd=cwd, shell=True)


  ##################### UTILS ###########################

    @classmethod
    def getDependencies(cls):
      # try to get CONDA activation command
      condaActivationCmd = cls.getCondaActivationCmd()
      neededProgs = []
      if not condaActivationCmd:
        neededProgs.append('conda')

      return neededProgs

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
    def getEnvSpecsPath(cls, env=None):
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
      return os.path.abspath(os.path.join(cls._rdkitHome, 'oddtModels', path))

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