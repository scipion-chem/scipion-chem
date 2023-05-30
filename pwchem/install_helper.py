from typing import List, Tuple, Dict
import pwem, os

class InstallHelper():
    """
    ### This class is intended to be used to ease the plugin installation process.

    #### Usage:
    InstallHelper class needs to be instanciated before it can be used.
    After that, commands can be chained together to run them in the defined order.
    The last command always needs to be addPackage().

    #### Example:
    installer = InstallHelper() # Instanciating class\n
    installer.getCloneCommand('test-package', '/home/user/myCustomPath', 'github.com/myRepo') # Cloning GitHub repository\n
    installer.getCondaenvCommand('test-package') # Creating conda enviroment\n
    installer.addPackage(env, 'test-package') # Install package\n

    #### It can also be done in a single line:
    installer.getCloneCommand('test-package', '/home/user/myCustomPath', 'github.com/myRepo').getCondaenvCommand('test-package').addPackage(env, 'test-package')\n

    #### If you want to check the command strings you are producing, use the function getCommandList() instead of addPackage() and assign it to a variable so you can print it.
    """
    # Global variables
    DEFAULT_VERSION = '1.0'

    def __init__(self, packageName: str, packageHome: str=None, packageVersion: str=DEFAULT_VERSION):
        """
        ### Constructor for the InstallHelper class.

        #### Parameters:
        packageName (str): Name of the package.
        packageHome (str): Optional. Path to the package. It can be absolute or relative to current directory.
        packageVersion (str): Optional. Package version.
        """
        # Private list of tuples containing commands with targets
        self.__commandList = []
        
        # Package name, version, and home
        self.__packageName = packageName
        self.__packageVersion = packageVersion
        self.__packageHome = packageHome if packageHome else os.path.join(pwem.Config.EM_ROOT, packageName + '-' + packageVersion)
    
    #--------------------------------------- PRIVATE FUNCTIONS ---------------------------------------#
    def __getTargetCommand(self, targetName: str) -> str:
        """
        ### This private function returns the neccessary command to create a target file given its name.
        ### Targets are always in uppercase and underscore format.

        #### Parameters:
        targetName (str): Name of the target file.

        #### Returns:
        (str): The command needed to create the target file.
        """
        return 'touch {}'.format(targetName)
    
    def __getBinaryEnvName(self, binaryName: str, binaryVersion: str=DEFAULT_VERSION) -> str:
        """
        ### This function returns the env name for a given package and repo.

        #### Parameters:
        binaryName (str): Name of the binary inside the package.
        binaryVersion (str): Optional. Binary's version.

        #### Returns:
        (str): The enviroment name for this binary.
        """
        return binaryName + "-" + binaryVersion
    
    def __getEnvActivationCommand(self, binaryName: str, binaryVersion: str=DEFAULT_VERSION) -> str:
        """
        ### Returns the conda activation command for the given enviroment.

        #### Parameters:
        binaryName (str): Name of the binary inside the package.
        binaryVersion (str): Optional. Version of the binary inside the package.

        #### Returns:
        (str): The enviroment activation command.
        """
        return "conda activate " + self.__getBinaryEnvName(binaryName, binaryVersion=binaryVersion)
    
    def __getBinaryNameAndVersion(self, binaryName: str=None, binaryVersion: str=None)  -> Tuple[str, str]:
        """
        ### Returns the binary name and version from an optionally introduced binary name and version.

        #### Parameters:
        binaryName (str): Name of the binary inside the package.
        binaryVersion (str): Optional. Version of the binary inside the package.

        #### Returns:
        tuple(str, str): The binary name and binary version.
        """
        binaryName = binaryName if binaryName else self.__packageName
        binaryVersion = binaryVersion if binaryVersion else self.__packageVersion
        return binaryName, binaryVersion
    
    #--------------------------------------- PUBLIC FUNCTIONS ---------------------------------------#
    def getCommandList(self) -> List[Tuple[str, str]]:
        """
        ### This function returns the list of commands with targets for debugging purposes or to export into another install helper.

        #### Returns:
        (list[tuple[str, str]]): Command list with target files.

        #### Usage:
        commandList = installer.getCommandList()
        """
        return self.__commandList
    
    def importCommandList(self, commandList: List[Tuple[str, str]]):
        """
        ### This function inserts the given formatted commands from another install helper into the current one.

        #### Parameters:
        commandList (list[tuple[str, str]]): List of commands generated by an install helper.

        #### Usage:
        installer1 = InstallHelper('package1', packageHome='/home/user/package2', packageVersion='1.0')
        installer1.addCommand('cd /home', 'CHANGED_DIRECTORY')
        installer2 = InstallHelper('package2', packageHome='/home/user/package2', packageVersion='1.0')
        installer2.importCommandList(installer1.getCommandList())

        #### Note:
        Argument 'packageHome' of the first installer must be the same as second installer.
        """
        # Adding given commands to current list
        self.__commandList.extend(commandList)
        return self
    
    def addCommand(self, command: str, targetName: str, workDir: str=''):
        """
        ### This function adds the given command with target to the command list.
        ### The target file needs to be located inside packageHome's directory so Scipion can detect it.

        #### Parameters:
        command (str): Command to be added.
        targetName (str): Name of the target file to be produced after commands are completed successfully.
        workDir (str): Optional. Directory where the command will be executed from.

        #### Usage:
        installer.addCommand('python3 myScript.py', 'MYSCRIPT_COMPLETED', workDir='/home/user/Documents/otherDirectory')

        #### This function call will generate the following commands:
        cd /home/user/Documents/otherDirectory && python3 myScript.py && touch /home/user/scipion/software/em/test-package-1.0/MYSCRIPT_COMPLETED
        """
        # Getting work directory
        workDirCmd = 'cd {} && '.format(workDir) if workDir else ''

        # Getting target name
        fullTargetName = os.path.join(self.__packageHome, targetName)

        command = (workDirCmd + command) if workDir else command
        self.__commandList.append((command + " && {}".format(self.__getTargetCommand(fullTargetName)), targetName))
        return self
    
    def addCommands(self, commandList: List[str], binaryName: str=None, workDir:str='', targetNames: List[str]=[]):
        """
        ### This function adds the given commands with targets to the command list.

        #### Parameters:
        commandList (list[str]): List containing the commands to add.
        binaryName (str): Optional. Name of the binary. Default is package name.
        workDir (str): Optional. Directory where the commands will be executed from.
        targetNames (list[str]): Optional. List containing the name of the target files for this commands.

        #### Usage:
        installer.addCommands(['python3 myScript.py', 'ls'], binaryName='myBinary', workDir='/home/user/Documents/otherDirectory',
            targetNames=['MYSCRIPT_COMPLETED', 'DIRECTORY_LISTED'])

        #### This function call will generate the following commands:
        cd /home/user/Documents/otherDirectory && python3 myScript.py && touch /home/user/scipion/software/em/test-package-1.0/MYSCRIPT_COMPLETED\n
        cd /home/user/Documents/otherDirectory && ls && touch /home/user/scipion/software/em/test-package-1.0/DIRECTORY_LISTED
        """
        # Defining binary name
        binaryName = self.__getBinaryNameAndVersion(binaryName=binaryName)[0]

        # Defining default target name preffix
        defaultTargetPreffix = '{}_EXTRA_COMMAND_'.format(binaryName.upper())

        # Executing commands
        for idx in range(len(commandList)):
            targetName = targetNames[idx] if targetNames else (defaultTargetPreffix + str(idx))
            self.addCommand(commandList[idx], targetName, workDir=workDir)

        return self
    
    def getCloneCommand(self, url: str, binaryFolderName: str='', targeName: str=None):
        """
        ### This function creates the neccessary command to clone a repository from Github.

        #### Parameters:
        url (str): URL to the git repository.
        binaryFolderName (str): Optional. Name of the binary directory.
        targetName (str): Optional. Name of the target file for this command.

        #### Usage:
        installer.getCloneCommand('https://github.com/myRepo.git', binaryFolderName='myCustomBinary', targeName='BINARY_CLONED')

        #### This function call will generate the following command:
        cd /home/user/scipion/software/em/test-package-1.0 && git clone https://github.com/myRepo.git myCustomBinary && touch BINARY_CLONED
        """
        # Defining target name
        targeName = targeName if targeName else '{}_CLONED'.format(binaryFolderName.upper())

        # Modifying binary name with a space for the command
        binaryFolderName = (' ' + binaryFolderName) if binaryFolderName else ''

        # Adding command
        self.addCommand('git clone {}{}'.format(url, binaryFolderName), targeName, workDir=self.__packageHome)

        return self
    
    def getCondaEnvCommand(self, binaryName: str=None, binaryPath: str=None, binaryVersion: str=None, pythonVersion: str=None, requirementsFile: bool=True,
                           requirementFileName: str='requirements.txt', requirementList: List[str]=[], extraCommands: List[str]=[], targetName: str=None):
        """
        ### This function creates the command string for creating a Conda enviroment and installing required dependencies for a given binary inside a package.

        #### Parameters:
        binaryName (str): Optional. Name of the binary. Default is package name.
        binaryPath (str): Optional. Path to the binary. It can be absolute or relative to current directory.
        binaryVersion (str): Optional. Binary's version. Default is package version.
        pythonVersion (str): Optional. Python version needed for the package.
        requirementsFile (bool): Optional. Defines if a Python requirements file exists.
        requirementFileName (bool): Optional. Name of the Python requirements file.
        requirementList (list[str]): Optional. List of Python packages to be installed. Can be used together with requirements file, but packages cannot be repeated.
        extraCommands (list[str]): Optional. List of extra conda-related commands to execute within the conda enviroment.
        targetName (str): Optional. Name of the target file for this command.

        #### Usage:
        installer.getCondaEnvCommand(binaryName='myBinary', binaryPath='/home/user/scipion/software/em/test-package-1.0/myBinary', binaryVersion='1.5', pythonVersion='3.11',
            requirementsFile=True, requirementFileName='requirements.txt', requirementList=['torch==1.2.0', 'numpy'],
            extraCommands=['conda info --envs'], targetName='CONDA_ENV_CREATED')

        #### This function call will generate the following command:
        eval "$(/home/user/miniconda/bin/conda shell.bash hook)"&& conda create -y -n myBinary-1.5 python=3.11 && conda activate myBinary-1.5 &&
        cd /home/user/scipion/software/em/test-package-1.0/myBinary && conda install pip -y && $CONDA_PREFIX/bin/pip install -r requirements.txt &&
        $CONDA_PREFIX/bin/pip install torch==1.2.0 numpyconda info --envs && cd /home/user/scipion/software/em/test-package-1.0 && touch CONDA_ENV_CREATED
        #### The path in the first command (eval ...) might vary, depending on the value of CONDA_ACTIVATION_CMD in your scipion.conf file.
        """
        # Binary name and version definition
        binaryName, binaryVersion = self.__getBinaryNameAndVersion(binaryName=binaryName, binaryVersion=binaryVersion)

        # Conda env creation
        createEnvCmd = 'conda create -y -n {}{}'.format(self.__getBinaryEnvName(binaryName, binaryVersion=binaryVersion), (' python={}'.format(pythonVersion)) if pythonVersion else '')

        # Command to install pip
        pipInstallCmd = 'conda install pip -y'

        # Command prefix for Python packages installation
        requirementPrefixCmd = '$CONDA_PREFIX/bin/pip install'

        # Requirements file name
        requirementFileName = os.path.join(binaryPath, requirementFileName) if requirementFileName and binaryPath else requirementFileName

        # Command for installing Python packages with requirements file
        installWithFile = (requirementPrefixCmd + ' -r ' + requirementFileName) if requirementsFile else ''

        # Command for installing Python packages manually
        installManual = ' '.join(requirementList)
        installManual = (requirementPrefixCmd + " " + installManual) if installManual else ''

        # Only install pip and Python packages if requiremenst file or manual list has been provided
        pythonCommands = (' && ' + pipInstallCmd) if (installWithFile or installManual) else ''
        if pythonCommands:
            pythonCommands += ' && {}'.format(installWithFile) if installWithFile else ''
            pythonCommands += ' && {}'.format(installManual) if installManual else ''
        
        # Defining target name
        targetName = targetName if targetName else '{}_CONDA_ENV_CREATED'.format(binaryName.upper())
        
        # Crafting final command string
        command = pwem.Plugin.getCondaActivationCmd() + ' ' + createEnvCmd                          # Basic commands: hook and env creation
        command += ' && ' + self.__getEnvActivationCommand(binaryName, binaryVersion=binaryVersion) # Env activation
        command += ' && cd {}'.format(binaryPath) if binaryPath else ''                             # cd to binary path if proceeds
        command += pythonCommands                                                                   # Python related commands
        command += " && ".join(extraCommands)                                                       # Extra conda commands
        command += ' && cd {}'.format(self.__packageHome) if binaryPath else ''                     # Return to package's root directory
        
        # Adding command
        self.addCommand(command, targetName)
        return self
    
    def addCondaPackages(self, packages: List[str], binaryName: str=None, binaryVersion: str=None, channel: str=None, targetName: str=None):
        """
        ### This function returns the command used for installing extra packages in a conda enviroment.

        #### Parameters:
        binaryName (str): Name of the binary. Default is package name.
        packages (list[str]): List of conda packages to install.
        binaryVersion (str): Optional. Binary's version. Default is package version.
        channel (str): Optional. Channel to download the package from.
        targetName (str): Optional. Name of the target file for this command.

        #### Usage:
        installer.addCondaPackages(packages=['pytorch==1.1.0', 'cudatoolkit=10.0'], binaryName='myBinary',
            binaryVersion='1.5', channel='conda-forge', targetName='CONDA_PACKAGES_INSTALLED')

        #### This function call will generate the following command:
        eval "$(/home/user/miniconda/bin/conda shell.bash hook)"&& conda activate myBinary-1.5 &&
        conda install -y pytorch==1.1.0 cudatoolkit=10.0 -c conda-forge && touch CONDA_PACKAGES_INSTALLED
        #### The path in the first command (eval ...) might vary, depending on the value of CONDA_ACTIVATION_CMD in your scipion.conf file.
        """
        # Binary name and version definition
        binaryName, binaryVersion = self.__getBinaryNameAndVersion(binaryName=binaryName, binaryVersion=binaryVersion)

        # Defininig target name
        targetName = targetName if targetName else '{}_CONDA_PACKAGES_INSTALLED'.format(binaryName.upper())

        # Adding installation command
        command = "{} {} && conda install -y {}".format(pwem.Plugin.getCondaActivationCmd(), self.__getEnvActivationCommand(binaryName, binaryVersion=binaryVersion), ' '.join(packages))
        if channel:
            command += " -c {}".format(channel)
        self.addCommand(command, targetName)

        return self
    
    def getExtraFile(self, url: str, targetName: str, location: str=".", workDir: str='', fileName: str=None):
        """
        ### This function creates the command to download with wget the file in the given link into the given path.
        ### The downloaded file will overwrite a local one if they have the same name.
        ### This is done to overwrite potential corrupt files whose download was not fully completed.

        #### Parameters:
        url (str): URL of the resource to download.
        targetName (str): Name of the target file for this command.
        location (str): Optional. Location where the file will be downloaded. It can be absolute or relative to current directory.
        workDir (str): Optional. Directory where the file will be downloaded from.
        fileName (str): Optional. Name of the file after the download. Use intended for cases when expected name differs from url name.

        #### Usage:
        installer.getExtraFile('https://site.com/myfile.tar', 'FILE_DOWNLOADED', location='/home/user/scipion/software/em/test-package-1.0/subdirectory', workDir='/home/user', fileName='test.tar')

        #### This function call will generate the following command:
        cd /home/user && mkdir -p /home/user/scipion/software/em/test-package-1.0/subdirectory &&
        wget -O /home/user/scipion/software/em/test-package-1.0/subdirectory/test.tar https://site.com/myfile.tar && touch /home/user/scipion/software/em/test-package-1.0/FILE_DOWNLOADED
        """
        # Getting filename for wget
        fileName = fileName if fileName else os.path.basename(url)
        mkdirCmd = "mkdir -p {} && ".format(location) if location else ''

        downloadCmd = "{}wget -O {} {}".format(mkdirCmd, os.path.join(location, fileName), url)
        self.addCommand(downloadCmd, targetName, workDir=workDir)
    
        return self

    def getExtraFiles(self, fileList: List[Dict[str, str]], binaryName: str=None, workDir: str='', targetNames: List[str]=None):
        """
        ### This function creates the command to download with wget the file in the given link into the given path.
        ### The downloaded file will overwrite a local one if they have the same name.
        ### This is done to overwrite potential corrupt files whose download was not fully completed.

        #### Parameters:
        fileList (list[dict[str, str, str]]): List containing files to be downloaded. Example: [{'url': url1, 'path': path1, 'name': 'test.tar'}, {'url': url2, 'path': path2, 'name': 'test2.tar'}]
        binaryName (str): Optional. Name of the binary.
        Each file is a list contaning url and location to download it. Paths can be an empty string for default location.
        workDir (str): Optional. Directory where the files will be downloaded from.
        targetNames (list[str]): Optional. List containing the name of the target files for this commands.

        #### Usage:
        installer.getExtraFiles(
            [
                {'url': 'https://site.com/myfile.tar', 'path': '/home/user/scipion/software/em/test-package-1.0/subdirectory1', 'name': 'test.tar'},
                {'url': 'https://site.com/myfile.tar2', 'path': '/home/user/scipion/software/em/test-package-1.0/subdirectory2', 'name': 'test2.tar2'}
            ],
            binaryName='myBinary', workDir='/home/user', targetNames=['DOWNLOADED_FILE_1', 'DOWNLOADED_FILE_2'])

        #### This function call will generate the following commands:
        cd /home/user && mkdir -p /home/user/scipion/software/em/test-package-1.0/subdirectory1 &&
        wget -O /home/user/scipion/software/em/test-package-1.0/subdirectory1/test.tar https://site.com/myfile.tar && touch /home/user/scipion/software/em/test-package-1.0/DOWNLOADED_FILE_1
        
        cd /home/user && mkdir -p /home/user/scipion/software/em/test-package-1.0/subdirectory2 &&
        wget -O /home/user/scipion/software/em/test-package-1.0/subdirectory2/test2.tar2 https://site.com/myfile.tar2 && touch /home/user/scipion/software/em/test-package-1.0/DOWNLOADED_FILE_2
        """
        # Defining binary name
        binaryName = self.__getBinaryNameAndVersion(binaryName=binaryName)[0]

        # Default preffix for target names
        defaultTargetPreffix = "{}_FILE_".format(binaryName.upper())

        # For each file in the list, download file
        for idx in range(len(fileList)):
            # Checking if file dictionary contains url
            if 'url' not in fileList[idx]:
                raise KeyError("ERROR: Download url has not been set for at least one file. You can create the appropiate dictionary calling function getFileDict.")
            
            # Getting proper file dictionary
            kwargs = {}
            if 'name' in fileList[idx]:
                kwargs['name'] = fileList[idx]['name']
            if 'path' in fileList[idx]:
                kwargs['path'] = fileList[idx]['path']
            downloadable = fileList[idx] if ('path' in fileList[idx] and 'name' in fileList[idx]) else self.getFileDict(fileList[idx]['url'], **kwargs)

            targetName = targetNames[idx] if targetNames else (defaultTargetPreffix + str(idx))
            self.getExtraFile(downloadable['url'], targetName, location=downloadable['path'], workDir=workDir, fileName=downloadable['name'])
    
        return self
    
    def addPackage(self, env, dependencies: List[str]=[], default: bool=True, **kwargs):
        """
        ### This function adds the given package to scipion installation with some provided parameters.
        
        #### Parameters:
        env: Scipion enviroment.
        dependencies (list[str]): Optional. List of dependencies the package has.
        default (bool): Optional. Defines if this package version is automatically installed with the plugin.
        **kwargs: Optional. Other possible keyword parameters that will be directly passed to env.addPackage.
        Intended for cases where multiple versions of the same package coexist in the same plugin.

        #### Usage:
        installer.addPackage(env, dependencies=['wget', 'conda'], default=True)
        """
        env.addPackage(self.__packageName, version=self.__packageVersion, tar='void.tgz', commands=self.__commandList, neededProgs=dependencies, default=default, **kwargs)
    
    #--------------------------------------- PUBLIC UTILS FUNCTIONS ---------------------------------------#
    def getFileDict(self, url: str, path: str='.', fileName: str=None) -> Dict[str, str]:
        """ This function generates the dictionary for a downloadable file. """
        # Getting file name
        fileName = fileName if fileName else os.path.basename(url)

        # Returning dictionary
        return {'url': url, 'path': path, 'name': fileName}