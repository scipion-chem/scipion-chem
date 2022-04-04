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

import os
from subprocess import run
import pyworkflow.utils as pwutils
from pyworkflow.tests import DataSet
import pwem
from .bibtex import _bibtexStr

from .constants import *

_logo = 'tool.png'
JCHEM, JCHEM_DEFAULT_VERSION = 'jchempaint', '3.3'
PYMOL, PYMOL_DEFAULT_VERSION = 'pymol', '2.5'
PLIP, PLIP_DEFAULT_VERSION = 'plip', '2.2'
RDKIT = 'rdkit'
ALIVIEW, ALIVIEW_DEFAULT_VERSION = 'aliview',  '1.28'

class Plugin(pwem.Plugin):
    _rdkitHome = os.path.join(pwem.Config.EM_ROOT, RDKIT)

    @classmethod
    def defineBinaries(cls, env):
        #PLIP environment (with pymol bundle)
        cls.addPLIPPackage(env, default=bool(cls.getCondaActivationCmd()))

        cls.addRDKitPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addMGLToolsPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addJChemPaintPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addPyMolPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addAliViewPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def _defineVariables(cls):
        cls._defineVar("RDKIT_ENV_ACTIVATION", 'conda activate rdkit-env')
        cls._defineVar("PLIP_ENV_ACTIVATION", 'conda activate plip-env')
        cls._defineEmVar('MGL_HOME', 'mgltools-1.5.6')
        cls._defineEmVar(PYMOL_HOME, 'pymol-2.5')

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
      plip_commands = 'conda create -y -n plip-env python=3.6 && '
      plip_commands += '%s %s && ' % (cls.getCondaActivationCmd(), cls.getPLIPEnvActivation())
      #Installing openbabel
      plip_commands += 'conda install -y -c conda-forge/label/cf202003 openbabel && '
      #Installing swig
      plip_commands += 'conda install -y -c conda-forge swig && '
      #Installing Pymol
      plip_commands += 'conda install -y -c schrodinger -c conda-forge pymol &&'
      #Installing PLIP
      plip_commands += 'conda install -y -c conda-forge plip && touch {}'.format(PLIP_INSTALLED)

      plip_commands = [(plip_commands, PLIP_INSTALLED)]
      env.addPackage('plip', version='2.2',
                     tar='void.tgz',
                     commands=plip_commands,
                     default=True)

    @classmethod
    def addPyMolPackage(cls, env, default=False):
      PYMOL_INSTALLED = 'pymol_installed'
      pymol_commands = 'wget https://pymol.org/installers/PyMOL-2.5.2_293-Linux-x86_64-py37.tar.bz2 -O {} && '.\
        format(cls.getPyMolTar())
      pymol_commands += 'tar -jxf {} --strip-components 1 && '.format(cls.getPyMolTar())
      pymol_commands += 'touch ' + PYMOL_INSTALLED
      pymol_commands = [(pymol_commands, PYMOL_INSTALLED)]
      env.addPackage(PYMOL, version=PYMOL_DEFAULT_VERSION,
                     tar='void.tgz',
                     commands=pymol_commands,
                     default=True)

    @classmethod
    def addMGLToolsPackage(cls, env, default=False):
      MGL_INSTALLED = "initMGLtools.sh"
      mgl_commands = 'wget https://ccsb.scripps.edu/download/548/ -O {} --no-check-certificate && '. \
        format(cls.getMGLPath('MGLTools.tar.gz'))
      mgl_commands += 'tar -xf {} --strip-components 1 && '.format(cls.getMGLPath('MGLTools.tar.gz'))
      mgl_commands += '{} && '.format(cls.getMGLPath('install.sh'))
      mgl_commands += 'touch ' + MGL_INSTALLED
      mgl_commands = [(mgl_commands, MGL_INSTALLED)]

      env.addPackage('mgltools', version='1.5.6',
                     tar='void.tgz',
                     commands=mgl_commands,
                     default=True)

    @classmethod
    def addJChemPaintPackage(cls, env, default=False):
        JCHEM_INSTALLED = 'jchem_installed'
        jchem_commands = 'wget https://github.com/downloads/JChemPaint/jchempaint/jchempaint-3.3-1210.jar -O {} && '.\
          format(cls.getJChemPath())
        jchem_commands += 'chmod +x {} && '.format(cls.getJChemPath())
        jchem_commands += ' touch %s' % JCHEM_INSTALLED

        jchem_commands = [(jchem_commands, JCHEM_INSTALLED)]

        env.addPackage('jchempaint', version='3.3',
                       tar='void.tgz',
                       commands=jchem_commands,
                       default=True)

    @classmethod
    def addRDKitPackage(cls, env, default=False):
        RDKIT_INSTALLED = 'rdkit_installed'

        installationCmd = cls.getCondaActivationCmd()
        installationCmd += ' conda create -y -c conda-forge -n rdkit-env rdkit oddt &&'
        installationCmd += ' mkdir oddtModels && touch %s' % RDKIT_INSTALLED

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
    def addAliViewPackage(cls, env, default=False):
      SEQS_INSTALLED = 'aliview_installed'
      seqs_commands = 'wget https://ormbunkar.se/aliview/downloads/linux/linux-version-1.28/aliview.tgz -O {} ' \
                      '--no-check-certificate && '.format(cls.getAliViewTar())
      seqs_commands += 'tar -xf {} && rm {} &&'.format(cls.getAliViewTar(), cls.getAliViewTar())
      seqs_commands += ' touch %s' % SEQS_INSTALLED

      seqs_commands = [(seqs_commands, SEQS_INSTALLED)]

      env.addPackage(ALIVIEW, version=ALIVIEW_DEFAULT_VERSION,
                     tar='void.tgz',
                     commands=seqs_commands,
                     default=True)


    ##################### RUN CALLS #################33333

    @classmethod
    def runRDKit(cls, protocol, program, args, cwd=None):
        """ Run rdkit command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getRDKitEnvActivation(), program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runRDKitScript(cls, protocol, scriptName, args, cwd=None):
      """ Run rdkit command from a given protocol. """
      scriptName = cls.getScriptsDir(scriptName)
      fullProgram = '%s %s && %s %s' % (cls.getCondaActivationCmd(), cls.getRDKitEnvActivation(), 'python', scriptName)
      protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def runJChemPaint(cls, protocol, cwd=None):
      """ Run jchempaint command from a given protocol. """
      protocol.runJob('java -jar {}'.format(cls.getJChemPath()), arguments='', env=cls.getEnviron(), cwd=cwd)

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

    @classmethod
    def getMGLPath(cls, path=''):
      return os.path.join(cls.getVar('MGL_HOME'), path)

    @classmethod
    def getJChemPath(cls):
        softwareHome = os.path.join(pwem.Config.EM_ROOT, JCHEM + '-' + JCHEM_DEFAULT_VERSION)
        return softwareHome + '/' + JCHEM + '-' + JCHEM_DEFAULT_VERSION + '.jar'

    @classmethod
    def getAliViewPath(cls):
      softwareHome = os.path.join(pwem.Config.EM_ROOT, ALIVIEW + '-' + ALIVIEW_DEFAULT_VERSION)
      return softwareHome

    @classmethod
    def getAliViewTar(cls):
      return cls.getAliViewPath() + '/' + ALIVIEW + '-' + ALIVIEW_DEFAULT_VERSION + '.tgz'

    @classmethod
    def getPyMolDir(cls):
      softwareHome = os.path.join(pwem.Config.EM_ROOT, PYMOL + '-' + PYMOL_DEFAULT_VERSION)
      return softwareHome

    @classmethod
    def getScriptsDir(cls, scriptName):
        return cls.getPluginHome('scripts/%s' % scriptName)

    @classmethod
    def getPyMolTar(cls):
        return cls.getPyMolDir() + '/' + PYMOL + '-' + PYMOL_DEFAULT_VERSION + '.tar.bz2'

    @classmethod
    def getPyMolPath(cls):
      return cls.getPyMolDir() + '/' + 'pymol'

    @classmethod
    def getMGLEnviron(cls):
      """ Create the needed environment for MGL Tools programs. """
      environ = pwutils.Environ(os.environ)
      pos = pwutils.Environ.BEGIN
      environ.update({
        'PATH': cls.getMGLPath('bin')
      }, position=pos)
      return environ

    @classmethod
    def getODDTModelsPath(cls, path=''):
      return os.path.abspath(os.path.join(cls._rdkitHome, 'oddtModels', path))


DataSet(name='smallMolecules', folder='smallMolecules',
            files={
              'mix': 'mix/',
              'mol2': 'mol2/',
              'pdb': 'pdb/',
              'sdf': 'sdf/',
              'smi': 'smi/'})