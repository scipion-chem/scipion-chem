from typing import List, Tuple
import pwem, os

class InstallHelper():
    """
    This class is intended to be used to ease the plugin installation process.
    """
    # Global variables
    DEFAULT_VERSION = '1.0'

    def __init__(self):
        """
        Constructor for the InstallHelper class.
        """
        # Private list of tuples containing commands with targets
        self.__commandList = []
    
    #--------------------------------------- PRIVATE FUNCTIONS ---------------------------------------#
    def __getTargetCommand(self, targetName: str) -> str:
        """
        This private function returns the neccessary command to create a target file given its name.
        Targets are always in uppercase and underscore format.

        Parameters:
        targetName (str): Name of the target file.

        Returns:
        (str): The command needed to create the target file.
        """
        return 'touch {}'.format(targetName)
    
    def __getBinaryEnvName(self, protocolName: str, version: str=DEFAULT_VERSION, binaryName: str=None) -> str:
        """
        This function returns the env name for a given protocol and repo.

        Parameters:
        protocolName (str): Name of the protocol.
        version (str): Binary's version.
        repoName (str): Optional. Name of the binary inside the protocol. Intended for protocols whose binaries' name differs from protocol's.

        Returns:
        (str): The enviroment name for this binary.
        """
        return (binaryName if binaryName else protocolName) + "-" + version
    
    def __getEnvActivationCommand(self, protocolName: str, binaryName: str=None, binaryVersion: str=DEFAULT_VERSION) -> str:
        """
        Returns the conda activation command for the given enviroment.

        Parameters:
        protocolName (str): Name of the protocol.
        binaryName (str): Optional. Name of the binary inside the protocol. Intended for protocols whose binaries' name differs from protocol's.
        binaryVersion (str): Optional. Version of the binary inside the protocol.
        """
        return "conda activate " + self.__getBinaryEnvName(protocolName, binaryVersion, binaryName)

    #--------------------------------------- PUBLIC FUNCTIONS ---------------------------------------#
    def getCommandList(self) -> List[Tuple[str, str]]:
        """
        This function returns the list of commands with targets for debugging purposes.

        Returns:
        (list[tuple[str, str]]): Command list with target files.
        """
        return self.__commandList

    def addCommand(self, command: str, targetName: str, workDir: str='', protocolPath: str=''):
        """
        This function adds the given command with target to the command list.

        Parameters:
        command (str): Command to be added.
        targetName (str): Name of the target file to be produced after commands are completed successfully.
        workDir (str): Optional. Directory where the command will be executed from.
        protocolPath (str): Optional. Protocol's root directory where target files are stored.
        """
        # Getting work directory
        workDirCmd = 'cd {} && '.format(workDir) if workDir else ''
        goBackCmd = ' && cd {}'.format(protocolPath if protocolPath else '-')

        command = (workDirCmd + command + goBackCmd) if workDir else command
        self.__commandList.append((command + " && {}".format(self.__getTargetCommand(targetName)), targetName))
        return self
    
    def addCommands(self, protocolName: str, commandList: List[str], binaryName: str=None, workDir:str='', protocolPath: str='', targetNames: List[str]=[]):
        """
        This function adds the given commands with targets to the command list.

        Parameters:
        protocolName (str): Name of the protocol.
        commandList (list[str]): List containing the commands to add.
        binaryName (str): Optional. Name of the binary.
        workDir (str): Optional. Directory where the commands will be executed from.
        protocolPath (str): Optional. Protocol's root directory where target files are stored.
        targetNames (list[str]): Optional. List containing the name of the target files for this commands.
        """
        # Defining binary name
        binaryName = binaryName if binaryName else protocolName

        # Defining default target name preffix
        defaultTargetPreffix = '{}_EXTRA_COMMAND_'.format(binaryName.upper())

        # Executing commands
        for idx in range(len(commandList)):
            targetName = targetNames[idx] if targetNames else (defaultTargetPreffix + str(idx))
            self.addCommand(commandList[idx], targetName, workDir, protocolPath)

        return self
    
    def getCloneCommand(self, protocolName: str, protocolHome: str, url: str, binaryFolderName: str=None, targeName: str=None):
        """
        This function creates the neccessary command to clone a repository from Github.

        Parameters:
        protocolName (str): Name of the protocol.
        protocolHome (str): Path to the protocol. It can be absolute or relative to current directory.
        url (str): URL to the git repository.
        binaryFolderName (str): Optional. Name of the binary directory.
        targetName (str): Optional. Name of the target file for this command.
        """
        # Defining binary name
        binaryFolderName = binaryFolderName if binaryFolderName else protocolName

        # Defining target name
        targeName = targeName if targeName else '{}_CLONED'.format(binaryFolderName.upper())

        # Adding command
        self.addCommand('git clone {} {}'.format( url, binaryFolderName), targeName, workDir=protocolHome)

        return self
    
    def getCondaEnvCommand(self, protocolName: str, binaryPath: str=None, binaryName: str=None, binaryVersion: str=DEFAULT_VERSION, pythonVersion: str=None, requirementsFile: bool=True,
                           requirementFileName: str='requirements.txt', requirementList: List[str]=[], extraCommands: List[str]=[], targetName: str=None):
        """
        This function creates the command string for creating a Conda enviroment and installing required dependencies for a given binary inside a protocol.

        Parameters:
        protocolName (str): Name of the protocol.
        binaryPath (str): Path to the binary. It can be absolute or relative to current directory.
        binaryName (str): Optional. Name of the binary.
        binaryVersion (str): Optional. Binary's version.
        pythonVersion (str): Optional. Python version needed for the protocol.
        requirementsFile (bool): Optional. Defines if a requirements file exists.
        requirementFileName (bool): Optional. Name of the requirements file.
        requirementList (list[str]): Optional. List of python packages to be installed. Can be used together with requirements file, but packages cannot be repeated.
        extraCommands (list[str]): Optional. List of extra conda-related commands to execute within the conda enviroment.
        targetName (str): Optional. Name of the target file for this command.
        """
        # Binary name definition
        binaryName = binaryName if binaryName else protocolName

        # Conda env creation
        createEnvCmd = 'conda create -y -n {}{}'.format(self.__getBinaryEnvName(protocolName, binaryVersion, binaryName), (' python={}'.format(pythonVersion)) if pythonVersion else '')

        # Requirements installation
        pipInstallCmd = 'conda install pip -y'
        requirementPrefixCmd = '$CONDA_PREFIX/bin/pip install'
        installWithFile = requirementPrefixCmd + ' -r ' + requirementFileName if requirementsFile else ''
        installManual = ' '.join(requirementList)
        installManual = (requirementPrefixCmd + " " + installManual) if installManual else ''
        finalInstallCmd = (' && ' + pipInstallCmd) if (installWithFile or installManual) else ''
        if finalInstallCmd:
            finalInstallCmd += ' && {}'.format(installWithFile) if installWithFile else ''
            finalInstallCmd += ' && {}'.format(installManual) if installManual else ''
        
        # Defining target name
        targetName = targetName if targetName else '{}_CONDA_ENV_CREATED'.format(binaryName.upper())
        
        # Adding conda commands
        self.addCommand('{} {} && {}{}{}{}{}'\
            .format(pwem.Plugin.getCondaActivationCmd(),
            createEnvCmd,
            self.__getEnvActivationCommand(protocolName, binaryName, binaryVersion),
            ' && cd {}'.format(binaryPath) if binaryPath else '',
            finalInstallCmd,
            " && ".join(extraCommands),
            ' && cd ..' if binaryPath else ''),
            targetName)

        return self
    
    def addCondaPackages(self, protocolName: str, packets: List[str], binaryName: str=None, binaryVersion: str=DEFAULT_VERSION, channel: str=None, targetName: str=None):
        """
        This function returns the command used for installing extra packages in a conda enviroment.

        Parameters:
        protocolName (str): Name of the protocol.
        packets (list[str]): List of conda packages to install.
        binaryName (str): Optional. Name of the binary.
        binaryVersion (str): Optional. Binary's version.
        channel (str): Optional. Channel to download the package from.
        targetName (str): Optional. Name of the target file for this command.
        """
        # Defining binary name
        binaryName = binaryName if binaryName else protocolName

        # Defininig target name
        targetName = targetName if targetName else '{}_CONDA_PACKAGES_INSTALLED'.format(binaryName.upper())

        # Adding installation command
        command = "{} {} && conda install -y {}".format(pwem.Plugin.getCondaActivationCmd(), self.__getEnvActivationCommand(protocolName, binaryName, binaryVersion), ' '.join(packets))
        if channel:
            command += " -c {}".format(channel)
        self.addCommand(command, targetName)

        return self
    
    def getExtraFile(self, url: str, targetName: str, location: str=".", workDir: str=''):
        """
        This function creates the command to download with wget the file in the given link into the given path.
        The downloaded file will overwrite a local one if they have the same name.
        This is done to overwrite potential corrupt files whose download was not fully completed.

        Parameters:
        url (str): URL of the resource to download.
        targetName (str): Name of the target file for this command.
        location (str): Optional. Location where the file will be downloaded.
        workDir (str): Optional. Directory where the file will be downloaded from from.
        """
        # Getting filename for wget
        fileName = os.path.basename(url)
        mkdirCmd = "mkdir -p {} && ".format(location) if location else ''
        location = location if location else '.'

        downloadCmd = "{}wget -O {}/{} {}".format(mkdirCmd, location, fileName, url)
        self.addCommand(downloadCmd, targetName, workDir)
    
        return self
    
    def getExtraFiles(self, protocolName: str, fileList: List[Tuple[str, str]], binaryName: str=None, workDir: str='', targetNames: List[str]=None):
        """
        This function creates the command to download with wget the file in the given link into the given path.
        The downloaded file will overwrite a local one if they have the same name.
        This is done to overwrite potential corrupt files whose download was not fully completed.

        Parameters:
        protocolName (str): Name of the protocol.
        fileList (list[tuple[str, str]]): List containing files to be downloaded. Example: [(url1, path1), (url2, path2)]
        binaryName (str): Optional. Name of the binary.
        Each file is a list contaning url and location to download it. Paths can be an empty string for default location.
        workDir (str): Optional. Directory where the files will be downloaded from.
        targetNames (list[str]): Optional. List containing the name of the target files for this commands.
        """
        # Defining binary name
        binaryName = binaryName if binaryName else protocolName

        # Default preffix for target names
        defaultTargetPreffix = "{}_FILE_".format(binaryName.upper())

        # For each file in the list, download file
        for idx in range(len(fileList)):
            targetName = targetNames[idx] if targetNames else (defaultTargetPreffix + str(idx))
            self.getExtraFile(fileList[idx][0], targetName, location=fileList[idx][1], workDir=workDir)
    
        return self
    
    def addProtocolPackage(self, env, protocolName: str, protocolVersion: str=DEFAULT_VERSION, dependencies: List[str]=[], default: bool=True):
        """
        This function adds the given protocol to scipion installation with some provided parameters.
        
        Parameters:
        env: Scipion enviroment.
        protocolName (str): Name of the protocol.
        protocolVersion (str): Protocol version.
        dependencies (list[str]): Optional. List of dependencies the protocol has.
        default (bool): Optional. Defines if this protocol version is automatically installed with the plugin.
        Intended for cases where multiple versions of the same protocol coexist in the same plugin.
        """
        env.addPackage(protocolName, version=protocolVersion, tar='void.tgz', commands=self.__commandList, neededProgs=dependencies, default=default)