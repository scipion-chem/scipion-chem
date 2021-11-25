# **************************************************************************
# *
# * Authors:     Carlos Oscar Sorzano (coss@cnb.csic.es)
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
import pyworkflow.utils as pwutils
from pyworkflow.tests import DataSet
import pwem
from .bibtex import _bibtexStr

_logo = 'tool.png'

class Plugin(pwem.Plugin):
    @classmethod
    def defineBinaries(cls, env):
        # The activation command should be something like
        # CONDA_ACTIVATION_CMD=eval "$(path-to-conda shell.bash hook)"
        # in the .config/scipion/scipion.conf file
        cls.addRDKitPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addopenbabelPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addMGLToolsPackage(env, default=bool(cls.getCondaActivationCmd()))
        cls.addJChemPaintPackage(env, default=bool(cls.getCondaActivationCmd()))

    @classmethod
    def _defineVariables(cls):
        cls._defineVar("RDKIT_ENV_ACTIVATION", 'conda activate my-rdkit-env')
        cls._defineVar("OPENBABEL_ENV_ACTIVATION", 'conda activate openbabel-env')
        cls._defineEmVar('MGL_HOME', 'mgltools-1.5.6')

    @classmethod
    def getRDKitEnvActivation(cls):
        activation = cls.getVar("RDKIT_ENV_ACTIVATION")
        return activation

    @classmethod
    def getOpenbabelEnvActivation(cls):
      activation = cls.getVar("OPENBABEL_ENV_ACTIVATION")
      return activation

    @classmethod
    def addMGLToolsPackage(cls, env, default=False):
      env.addPackage('mgltools', version='1.5.6',
                     url='http://mgltools.scripps.edu/downloads/downloads/tars/releases/REL1.5.6/mgltools_x86_64Linux2_1.5.6.tar.gz',
                     buildDir='mgltools_x86_64Linux2_1.5.6',
                     commands=[("./install.sh", "initMGLtools.sh")],
                     default=True)

    @classmethod
    def addJChemPaintPackage(cls, env, default=False):
        jchem_commands = 'wget https://github.com/downloads/JChemPaint/jchempaint/jchempaint-3.3-1210.jar -O {} && '.\
          format(cls.getJChemPath())
        jchem_commands += 'chmod +x {}'.format(cls.getJChemPath())

        env.addPackage('jchempaint', version='3.3',
                       commands=jchem_commands,
                       default=True)

    @classmethod
    def addopenbabelPackage(cls, env, default=False):
      OPENBABEL_INSTALLED = 'openbabel_installed'

      # try to get CONDA activation command
      installationCmd = cls.getCondaActivationCmd()

      # Create the environment
      installationCmd += ' conda create -y -c conda-forge -n openbabel-env openbabel &&'

      # Flag installation finished
      installationCmd += ' touch %s' % OPENBABEL_INSTALLED

      openbabel_commands = [(installationCmd, OPENBABEL_INSTALLED)]

      envPath = os.environ.get('PATH', "")
      installEnvVars = {'PATH': envPath} if envPath else None
      env.addPackage('openbabel',
                     tar='void.tgz',
                     commands=openbabel_commands,
                     neededProgs=cls.getDependencies(),
                     default=default,
                     vars=installEnvVars)

    @classmethod
    def runOPENBABEL(cls, protocol, program="obabel", args=None, cwd=None):
      """ Run openbabel command from a given protocol. """
      fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getOpenbabelEnvActivation(), program)
      protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd, numberOfThreads=1)

    @classmethod
    def getDependencies(cls):
        # try to get CONDA activation command
        condaActivationCmd = cls.getCondaActivationCmd()
        neededProgs = []
        if not condaActivationCmd:
            neededProgs.append('conda')

        return neededProgs

    @classmethod
    def addRDKitPackage(cls, env, default=False):
        RDKIT_INSTALLED = 'rdkit_installed'

        # try to get CONDA activation command
        installationCmd = cls.getCondaActivationCmd()

        # Create the environment
        installationCmd += ' conda create -y -c conda-forge -n my-rdkit-env rdkit &&'

        # Flag installation finished
        installationCmd += ' touch %s' % RDKIT_INSTALLED

        rdkit_commands = [(installationCmd, RDKIT_INSTALLED)]

        envPath = os.environ.get('PATH', "")  # keep path since conda likely in there
        installEnvVars = {'PATH': envPath} if envPath else None
        env.addPackage('rdkit',
                       tar='void.tgz',
                       commands=rdkit_commands,
                       neededProgs=cls.getDependencies(),
                       default=default,
                       vars=installEnvVars)

    @classmethod
    def getPluginHome(cls, path=""):
        import pwchem
        fnDir = os.path.split(pwchem.__file__)[0]
        return os.path.join(fnDir,path)

    @classmethod
    def runRDKit(cls, protocol, program, args, cwd=None):
        """ Run rdkit command from a given protocol. """
        fullProgram = '%s %s && %s' % (cls.getCondaActivationCmd(), cls.getRDKitEnvActivation(), program)
        protocol.runJob(fullProgram, args, env=cls.getEnviron(), cwd=cwd)

    @classmethod
    def getMGLPath(cls, path=''):
      return os.path.join(cls.getVar('MGL_HOME'), path)

    @classmethod
    def getJChemPath(cls):
      return cls.getPluginHome('jChemPaint_3.3.jar')

    @classmethod
    def getMGLEnviron(cls):
      """ Create the needed environment for MGL Tools programs. """
      environ = pwutils.Environ(os.environ)
      pos = pwutils.Environ.BEGIN
      environ.update({
        'PATH': cls.getMGLPath('bin')
      }, position=pos)
      return environ


DataSet(name='smallMolecules', folder='smallMolecules',
            files={
              'mix': 'mix/',
              'mol2': 'mol2/',
              'pdb': 'pdb/',
              'sdf': 'sdf/',
              'smi': 'smi/'})