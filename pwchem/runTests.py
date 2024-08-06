import subprocess, argparse, multiprocessing, sys, json, os, copy
from collections.abc import Callable
from typing import List, Union, Tuple, Dict

################################## INTERNAL UTILS FUNCTIONS ##################################
def removeLooseDependencyTests(testList: List[str], dependantTests: Dict) -> Tuple[List[str], Dict]:
	"""
	### This function removes all dependant tests where at least one of its dependencies cannot be met.

	#### Params:
	- testList (list(str)): List of tests.
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (list(str)): Updated list of tests.
	- (dict): Updated list of tests that depend on other tests and can be run.
	"""
	# Get test name and dependencies for each test
	for test, deps in list(dependantTests.items()):
		# Check if each dependency is in the main list
		for dep in deps:
			if dep not in testList:
				# If dependency cannot be met, delete test from list
				if test in testList:
					testList.remove(test)
					printSkippingTest(test, skipType='other', reason=f"Missing dependency with test {dep}")
				del dependantTests[test]

	return testList, dependantTests

def removeLooseDependencyTestsCascade(testList: List[str], dependantTests: Dict) -> Tuple[List[str], Dict]:
	"""
	### This function removes all dependant tests where at least one of its dependencies cannot be met.
	### It also handles cascading removals.

	#### Params:
	- testList (list(str)): List of tests.
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (list(str)): Updated list of tests.
	- (dict): Updated list of tests that depend on other tests and can be run.
	"""
	while True:
		# Call removal function until it stops removing tests
		initialLen = len(testList)
		testList, dependantTests = removeLooseDependencyTests(testList, dependantTests)
		if len(testList) == initialLen:
			break

	return testList, dependantTests

def removeDependentTests(testList: List[str], dependantTests: Dict) -> Tuple[List[str], Dict]:
	"""
	### This function removes the tests that depend on other tests that cannot be ran.

	#### Params:
	- testList (list(str)): List of tests.
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (list(str)): Updated list of tests.
	- (dict): List of tests that deppend on other tests and can be run.
	"""
	# Remove tests from the dictionary where at least one of its dependencies is not in the main test list
	testList, dependantTests = removeLooseDependencyTestsCascade(testList, dependantTests)

	# Find all tests with a circular dependency
	circularTests = {}
	for test in dependantTests:
		path = findCircularDependency(test, dependantTests)
		if path and test in path:
			circularTests[test] = path
	
	# Remove all tests with a circular dependency
	for test in circularTests.keys():
		if test in testList:
			testList.remove(test)
			depStr = ' --> '.join(circularTests[test])
			printSkippingTest(test, skipType='other', reason=f"It has a circular dependency: {depStr}")
		del dependantTests[test]

	# Remove tests that depend on removed tests (cascading removals)
	testList, dependantTests = removeLooseDependencyTestsCascade(testList, dependantTests)

	return testList, dependantTests

def getFirstExecutionLevel(dependantTests: Dict) -> List[str]:
	"""
	### This function returns the tests of the first execution level given a dictionary with test dependencies.

	#### Params:
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (list(str)): List of the tests of the first level.
	"""
	# Initialize list
	levelList = []

	# Create a copy of the dictionary to iterate over it
	dependantTestsCopy = copy.deepcopy(dependantTests)

	# For every test, check if at least one of its dependencies is another test of the list
	for test, deps in list(dependantTestsCopy.items()):
		if not any(dep in dependantTests.keys() for dep in deps):
			levelList.append(test)
			del dependantTestsCopy[test]
	
	return levelList

################################## GENERAL UTILS FUNCTIONS ##################################
def colorStr(string: str, color: str) -> str:
	"""
	### This function returns the input string wrapped in the specified color.

	#### Params:
	- string (str): Text to encase in a different color.
	- color (str): Color to encase the string into. Values can be "green", "yellow", "red", or "blue.

	#### Return:
	- (str): Text encased in given color.
	"""
	if color == 'green':
		return f"\033[92m{string}\033[0m"
	elif color == 'yellow':
		return f"\033[93m{string}\033[0m"
	elif color == 'red':
		return f"\033[91m{string}\033[0m"
	elif color == 'blue':
		return f"\033[94m{string}\033[0m"
	else:
		return string

def runJob(cmd: List[str]) -> Tuple[int, str]:
	"""
	### This function runs the given command.

	#### Params:
	- cmd (list(str)): Command to run formatted in ["executable", "params"] format.

	#### Return:
	- (int): Return code.
	- (str): Output text of the command.
	"""
	try:
		result = subprocess.run(cmd, check=True, capture_output=True, text=True)
		# Capture return code and message
		retCode = result.returncode
		output = result.stdout if retCode == 0 else result.stderr
	except subprocess.CalledProcessError as e:
		# Unknown error while running command
		output = e.stderr
		retCode = 1
	return retCode, output

def runInParallel(func: Callable, *args, paramList: List[str], jobs: int) -> List[str]:
	"""
	### This function creates a pool of workers to run the given function in parallel.
	### Also returns a list with the failed commands.

	#### Params:
	- func (Callable): Function to apply in parallel.
	- paramList (list(str)): List of param strings. One for each parallel call.
	- jobs (int): Number of parallel threads to use.

	#### Return:
	- (list(str)): List of failed commands.
	"""
	# Create a pool of worker processes
	nJobs = len(paramList) if len(paramList) < jobs else jobs
	pool = multiprocessing.Pool(processes=nJobs)

	# Apply the given function to the given param list using the pool
	results = [pool.apply_async(func, args=(param,*args,)) for param in paramList]

	# Initializing list of failed commands
	failedCommands = []

	# Check if any process encountered an error
	for result in results:
		if result.get():
			failedCommands.append(result.get())
	
	# Close the pool to release resources
	pool.close()

	# Join the pool, waiting for all processes to complete
	pool.join()

	# Only return list of failed commands
	return failedCommands

def runTest(test: str, scipionExecutable: str, testPrefix: str) -> Union[str, None]:
	"""
	### This function receives a test and runs it.

	#### Params:
	- test (str): Test name.
	- scipionExecutable (str): Path to Scipion's executable.
	- testPrefix (str): Prefix of the test.

	#### Return:
	- (None | str): None if everything worked properly, or the test name otherwise.
	"""
	printAndFlush(f"Running test {test}...")
	# Running test
	retCode, output = runJob([scipionExecutable, testPrefix + test])
	if retCode == 0:
		printAndFlush(colorStr(f"Test {test} OK", color='green'))
	else:
		# Detect failed test
		printAndFlush(output)
		printAndFlush(colorStr(f"Test {test} failed with above message.", color='red'))
		return test

def runTestBatch(testList: List[str], jobs: int, scipion: str, testPrefix: str, plugin: str) -> List[str]:
	"""
	### This function runs the given list of tests and returns the failed ones.

	#### Params:
	- testList (list(str)): List of tests to run.
	- testList (int): Max number of jobs to use.
	- scipion (str): Path to Scipion's executable.
	- testPrefix (str): Prefix of the test.
	- plugin (str): Module name of the plugin.

	#### Return:
	- (list(str)): List of all the failed tests.
	"""
	# Showing message with number of tests
	testSize = len(testList)
	nJobs = testSize if testSize < jobs else jobs
	testWord = 'test' + ('s' if testSize > 1 else '')
	batchMessage = f' in batches of {nJobs} processes' if testSize > 1 else ''
	printAndFlush(colorStr(f"Running a total of {testSize} {testWord} for {plugin}{batchMessage}...", color='blue'))

	# Run all the tests in parallel
	return runInParallel(runTest, scipion, testPrefix, paramList=testList, jobs=nJobs)

def printAndFlush(message: str):
	"""
	### This function prints a given message flushing stdout.

	#### Params:
	- message (str): Text to print.
	"""
	print(message, flush=True)

def printSkippingTest(test: str, skipType: str='', dependency: str='', reason: str=''):
	"""
	### This function prints the message when a test is skipped.
	#### Note: param dependency is only used when skipType is 'dependency', and param reason is only used when skipType is 'other'.

	#### Params:
	- test (str): Test to skip.
	- skipType (str): Optional. Type of predefined reason to skip the test.
	- dependency (str): Optional. Possible dependency of the test with another package.
	- reason (str): Optional. Custom reason to skip a test.
	"""
	# Defining base message
	baseMessage = f"Skipping test {test}. "
	if skipType != 'other' or reason:
		baseMessage += "Reason: "

	if skipType == 'gpu':
		# Printing custom message for gpu tests
		printAndFlush(colorStr(f"{baseMessage}Needs GPU.", color='yellow'))
	elif skipType == 'dependency':
		# Printing custom message for dependency tests
		dependencySpec = f" with plugin {dependency}" if dependency else ''
		printAndFlush(colorStr(f"{baseMessage}Unmet dependency{dependencySpec}.", color='yellow'))
	elif skipType == 'other':
		reason = reason if reason else "No reason provided"
		printAndFlush(colorStr(f"{baseMessage}{reason}.", color='yellow'))
	else:
		# Error message when skip type is not recogniced
		printFatalError(f"ERROR: Test skip type \"{skipType}\" not recognized.")

def findCircularDependency(test: str, dependantTests: List[Dict], path: List[str]=None) -> Union[List[str], bool]:
	"""
	### This function returns a list of all the tests involved in a circular path, or False if there is no such path.

	#### Params:
	- test (str): Name of the test being checked for a circular dependency.
	- dependantTests (dict): Dictionary containing tests that depend on other tests.
	- path (list(str)): Optional: List of tests that might form a circular dependency.

	#### Return:
	- (list(str)): List with the test forming a circular path, False if there is no such path.
	"""
	# Initialize circular path
	path = [] if path is None else path

	# If test is in path, we have a circular path
	if test in path:
		# Only return the part of the path that starts with current test
		# and ends with such test
		return path[path.index(test):] + [test]
	# Add test to path (visit node)
	path.append(test)

	# Check if every dependency is in a circular path
	for dep in dependantTests.get(test, []):
		circularPath = findCircularDependency(dep, dependantTests, path=path[:])
		# If the return is not False, then it's a non-empty list with a circular path
		if circularPath:
			return circularPath
	return False

def printFatalError(message: str):
	"""
	### This function prints the given message in red, flushes stdout, and exits with code 1.

	#### Params:
	- message (str): Error text to show.
	"""
	printAndFlush(colorStr(message, color='red'))
	sys.exit(1)

def getCondaActivationCmd(scipionExecutable: str) -> Union[str, None]:
	"""
	### This function reads the scipion config file and returns the conda activation command.
	
	#### Params:
	- scipionExecutable (str): Path to Scipion's executable.

	#### Return:
	- (str | None): Conda activation command if exists, else None.
	"""
	# Path to scipion config file
	configFile = os.path.join(os.path.dirname(scipionExecutable), 'config', 'scipion.conf')

	# Read file and return conda activation command
	with open(configFile, 'r') as f:
		for line in f:
			if line.startswith('CONDA_ACTIVATION_CMD'):
				return line.split('=')[1].strip()

def testPythonCommand(scipion: str, pythonCommand: str) -> bool:
	"""
	### This function executes the given Python command within scipion3 env and returns True if it succeeded or False if it failed.

	#### Params:
	- scipion (str): Path to Scipion's executable.
	- pythonCommand (str): Python command to test.

	#### Return:
	- (bool): True if command succeeded, False otherwise.
	"""
	command = f"{getCondaActivationCmd(scipion)} && conda activate scipion3 && python -c '{pythonCommand}'"
	result = subprocess.run(command, shell=True)
	return bool(result.returncode)

################################## MAIN EXECUTION FUNCTIONS ##################################
def getAllTests(scipion: str, pluginModule: str, testPrefix: str) -> List[str]:
	"""
	### This function finds the full list of tests from a given module.

	#### Params:
	- scipion (str): Path to Scipion's executable.
	- pluginModule (str): Module name of the plugin.
	- testPrefix (str): Prefix for the test name.

	#### Return:
	- (list(str)): List of tests names.
	"""
	# Construct the test discovery command
	command = f"{scipion} test --grep {pluginModule}."

	# Run the shell command and capture the output
	try:
		output = subprocess.check_output(command, shell=True, text=True)
	except subprocess.CalledProcessError:
		printFatalError("ERROR: Test search command failed. Check line above for more detailed info.")

	# Define test command string variables
	scipionTestsStartingSpaces = '   '

	# Separate lines into a list
	lines = output.split('\n')

	# For each line, keep only the ones about individual tests in a reduced form
	filteredLines = []
	for line in lines:
		if line.startswith(scipionTestsStartingSpaces):
			filteredLines.append(line.replace(f'{scipionTestsStartingSpaces}scipion3 {testPrefix}', ''))
	
	# If no tests were found, check if module was not found or if plugin has no tests
	if not filteredLines:
		# If import caused an error, module was not found
		if not testPythonCommand(scipion, f"import {pluginModule}"):
			printFatalError(f"ERROR: No tests were found for module {args.plugin}. Are you sure this module is properly installed?")
	
	# Return full list of tests
	return filteredLines

def readTestDataFile(testDataFilePath: str) -> Tuple[Union[List[str], None], Union[Dict, None], Union[List[Dict], None]]:
	"""
	### This function returns a list with the necessary datasets for the tests, as well as an object with the different tests and the situations where to skip them.

	#### Params:
	- testDataFilePath (str): Path to testData json file.

	#### Return:
	- (list(str) | None): List of dataset names if there are any, None otherwise.
	- (dict | None): Dictionary containing all skippable tests if there are any, None otherwise.
	- (list(dict)): List of tests dependant on other tests.
	"""
	# Read the JSON data from the file
	try:
		with open(testDataFilePath, 'r') as file:
			dataFile = json.load(file)
			return dataFile.get("datasets", []), dataFile.get("skippable", {}), dataFile.get("test-dependencies", {})
	except FileNotFoundError:
		printAndFlush(colorStr("No skippable tests file found, running all.", color='yellow'))
		return None, None, None
	except PermissionError:
		printFatalError(f"ERROR: Permission denied to open file '{testDataFilePath}'.")
	except json.JSONDecodeError as e:
		printFatalError(f"ERROR: Invalid JSON format in file '{testDataFilePath}':\n{e}")
	except Exception as e:
		printFatalError(f"An unexpected error occurred:\n{e}")
	
def sortDependentTests(dependantTests: Dict) -> Dict:
	"""
	### This function sorts the tests that depend on other tests into the order they must be ran.

	#### Params:
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (dict): Dictionary of tests that deppend on other tests, sorted in execution order.
	"""
	# Check if received dictionary contains data
	if not dependantTests:
		return {}

	# Initialize empty dictionary
	executionLevels = {}

	# Obtain execution level in a loop until there are no more tests
	level = 0
	while dependantTests:
		levelList = getFirstExecutionLevel(dependantTests)
		for test in levelList:
			executionLevels[level] = levelList
			del dependantTests[test]
		level += 1

	return executionLevels

def cleanAndSortDependentTests(testList: List[str], dependantTests: Dict) -> Tuple[List[str], Dict]:
	"""
	### This function sorts the tests that depend on other tests into the order they must be ran, skipping the ones that cannot be ran.

	#### Params:
	- testList (list(str)): List of tests.
	- dependantTests (dict): Dictionary containing tests that depend on other tests.

	#### Return:
	- (list(str)): Updated list of tests.
	- (dict): Dictionary of tests that deppend on other tests and can be run, sorted in execution order.
	"""
	# Remove tests that do not meet dependency criteria
	testList, dependantTests = removeDependentTests(testList, dependantTests)

	# Build execution levels, to run after the main list
	executionLevels = sortDependentTests(dependantTests)

	# Remove tests in the execution levels dictionary from the main list
	testWithLevels = []
	for level in executionLevels.values():
		for test in level:
			testWithLevels.append(test)
	testList = [test for test in testList if test not in testWithLevels]
	
	return testList, executionLevels

def downloadDatset(dataset: str, scipionExecutable: str) -> Union[str, None]:
	"""
	### This function downloads a given dataset for scipion tests.

	#### Params:
	- dataset (str): Dataset name.
	- scipionExecutable (str): Path to Scipion's executable.

	#### Return:
	- (str | None): Dataset name if command failed, None otherwise.
	"""
	printAndFlush(f"Downloading dataset {dataset}...")
	retCode, output = runJob([scipionExecutable, f" testdata --download {dataset}"])
	if retCode == 0:
		printAndFlush(colorStr(f"Dataset {dataset} download OK", color='green'))
	else:
		# Detect failed test
		printAndFlush(output)
		printAndFlush(colorStr(f"Dataset {dataset} download failed with the above message.", color='red'))
		return dataset
	
def removeGPUTests(testList: List[str], noGPU: bool, gpuSkippableTests: List[str]) -> List[str]:
	"""
	### This function removes the GPU tests if flag noGPU has been set.

	#### Params:
	- testList (list(str)): List of tests.
	- noGPU (bool): Whether to skip GPU tests or not.
	- gpuSkippableTests (list(str)): List of tests skippable in noGPU mode.

	#### Return:
	- (list(str)): List of tests without GPU based ones if they need to be removed.
	"""
	# Check if noGPU flag was set, only remove if so
	if noGPU:
		for gpuTest in gpuSkippableTests:
			# For each test, remove from list if it was there (dev might have mispelled the test on the .json file)
			if gpuTest in testList:
				printSkippingTest(gpuTest, skipType='gpu')
				testList.remove(gpuTest)
	
	# Return modified list of tests
	return testList

def removeDependecyTests(scipion: str, testList: List[str], dependenciesSkippableTests: List[Dict]) -> List[str]:
	"""
	### This function removes the depencendy related tests if the dependencies are not met.

	#### Params:
	- scipion (str): Path to Scipion's executable.
	- testList (list(str)): List of tests.
	- dependenciesSkippableTests (list(dict)): List of tests skippable if dependencies are not met.

	#### Return:
	- (list(str)): List of tests without dependency limited ones if dependencies are not met.
	"""
	# Removing dependency tests if dependencies are not present
	for dependency in dependenciesSkippableTests:
		# Getting plugin and module name
		dependencyTestKey = dependency.get("name", None)
		dependencyTestModule = dependency.get("module", None)

		# Try to import module if provided
		if dependencyTestModule:
			# Creating import command to run within scipion3 conda env
			if testPythonCommand(scipion, f"import {dependencyTestModule}"):
				continue

		# If no module was provided or import raised a ModuleNotFoundError exception, skip tests
		for dependencyTest in dependency.get("tests", []):
			if dependencyTest in testList:
				printSkippingTest(dependencyTest, skipType='dependency', dependency=dependencyTestKey)
				testList.remove(dependencyTest)
	
	# Return modified list of tests
	return testList

def removeOtherTests(testList: List[str], otherSkippableTests: List[Dict]) -> List[str]:
	"""
	### This function removes the other skippable tests.

	#### Params:
	- testList (list(str)): List of tests.
	- otherSkippableTests (list(dict)): List of other skippable tests.

	#### Return:
	- (list(str)): List of tests without other skippable ones.
	"""
	for otherTest in otherSkippableTests:
		# Getting test name and reason
		testName = otherTest.get("test", None)
		testReason = otherTest.get("reason", None)

		# If test exists, skip
		if testName and testName in testList:
			printSkippingTest(testName, skipType='other', reason=testReason)
			testList.remove(testName)
	
	# Return modified test list
	return testList

def removeSkippableTests(scipion: str, testList: List[str], noGPU: bool, gpuSkippableTests: List[str], dependenciesSkippableTests: List[Dict], otherSkippableTests: List[Dict]) -> List[str]:
	"""
	### This function removes from the list of all tests, the ones that have to be skipped.

	#### Params:
	- scipion (str): Path to Scipion's executable.
	- testList (list(str)): List of tests.
	- noGPU (bool): Whether to skip GPU tests or not.
	- gpuSkippableTests (list(str)): List of tests skippable in noGPU mode.
	- dependenciesSkippableTests (list(dict)): List of tests skippable if dependencies are not met.
	- otherSkippableTests (list(dict)): List of other skippable tests.

	#### Return:
	- (list(str)): List of tests without skipped ones.
	"""
	# Remove GPU skippable tests
	testList = removeGPUTests(testList, noGPU, gpuSkippableTests)

	# Removing dependency tests if dependencies are not present
	testList = removeDependecyTests(scipion, testList, dependenciesSkippableTests)

	# Removing other tests for reasons stated
	testList = removeOtherTests(testList, otherSkippableTests)
	
	return testList

def getResultDictionary(testList: List[str], failedTests: List[str]) -> Dict:
	"""
	### This function returns a dictionary with all the passed/failed tests grouped by their origin file.

	#### Params:
	- testList (list(str)): List of tests.
	- failedTests (list(str)): List of failed tests.

	#### Return:
	- (dict): Dictionary of results.
	"""
	groupedTests = {}

	# Group the tests based on the file they are declared on
	for test in testList:
		key, value = test.split('.', 1)
		groupedTests.setdefault(key, []).append(value)

	# Create a new dictionary to separate from the grouped tests,
	# the ones failing from the ones passing
	results = {}
	for key, valueList in groupedTests.items():
		# Initialize the 'failed' and 'passed' lists for each key
		results[key] = {'failed': [], 'passed': []}
		
		# Check if each test is in the failed tests list
		for test in valueList:
			resultListKey = 'failed' if f"{key}.{test}" in failedTests else 'passed'
			results[key][resultListKey].append(test)
	
	return results

def printSummary(results: Dict):
	"""
	### This function prints the sumarry of the test results.

	#### Params:
	- results (dict): Dictionary of results.
	"""
	printAndFlush("SUMMARY:")
	for testGroup in results:
		# Calculating passed and failed tests for each group
		passed = results.get(testGroup, {}).get("passed", [])
		failed = results.get(testGroup, {}).get("failed", [])
		total = len(passed) + len(failed)
		printAndFlush(f"{testGroup}: [{len(passed)} / {total}]")
		if failed:
			# If there are any failed tests, show which
			printAndFlush(colorStr(f"\tFailed tests: {' '.join(failed)}", color='red'))

################################## MAIN FUNCTION ##################################
def main(args: Dict):
	"""
	### Main execution function.

	#### Params:
	- args (dict): Dictionary containing all command-line arguments.
	"""
	# Defining test prefix
	testPrefix = f'tests {args.plugin}.tests.'

	# Getting full list of tests
	filteredLines = getAllTests(args.scipion, args.plugin, testPrefix)

	# If test list is empty, plugin has no tests
	if not filteredLines:
		# This case is considered a sucess since nothing actually failed
		printAndFlush(colorStr(f"Module {args.plugin} has not tests. Nothing to run.", color='yellow'))
		sys.exit(0)

	# Obtaining datasets, skippable tests, and dependant tests according to situation
	datasets, allSkippableTests, dependantTests = readTestDataFile(args.testData)

	if allSkippableTests:
		gpuSkippableTests = allSkippableTests.get('gpu', [])
		dependenciesSkippableTests = allSkippableTests.get('dependencies', [])
		otherSkippableTests = allSkippableTests.get('others', [])

		# Removing skippable tests
		filteredLines = removeSkippableTests(args.scipion, filteredLines, args.noGPU, gpuSkippableTests, dependenciesSkippableTests, otherSkippableTests)

	# Obtain tests with dependencies and add them to a separete list
	# (or skip them if their dependencies cannot be met or they have a circular dependency)
	if dependantTests:
		filteredLines, dependantTests = cleanAndSortDependentTests(filteredLines, dependantTests)

	# Downloading in parallel required datasets if there are any
	if datasets:
		printAndFlush(colorStr(f"Downloading {len(datasets)} datasets...", color="blue"))
		failedDownloads = runInParallel(downloadDatset, args.scipion, paramList=datasets, jobs=len(datasets))

		# Check if there were any errors
		if failedDownloads:
			printFatalError("The download of at least one dataset ended with errors. Exiting.")

	# Showing initial message abount main tests and run them
	if dependantTests:
		printAndFlush(colorStr("Initial run of non-dependent tests.", color="blue"))
	failedTests = runTestBatch(filteredLines, args.jobs, args.scipion, testPrefix, args.plugin)

	# Run all dependent test in batches for each level and collect failed ones
	for level in dependantTests.keys():
		# Showing batch execution message and run it
		printAndFlush(colorStr(f"Batch of dependent tests {level+1}/{len(dependantTests.keys())}.", color="blue"))
		levelFailedTests = runTestBatch(dependantTests[level], args.jobs, args.scipion, testPrefix, args.plugin)
		failedTests.extend(levelFailedTests)
		filteredLines.extend(dependantTests[level])

	# Get results grouped by origin file and separated into passed and failed
	results = getResultDictionary(filteredLines, failedTests)

	# Print summary of passed/failed tests
	printSummary(results)

	# Check if an error occurred
	if failedTests:
		printFatalError("Some tests ended with errors. Exiting.")

	# Message if all tests succeeded
	printAndFlush(colorStr("\nAll test passed!", color='green'))

if __name__ == "__main__":
	""" Calls main function when executed. """
	# Parse the command-line arguments
	epilog = "Example 1: python script.py /path/to/scipion pwchem -j 2"
	epilog += "\nExample 2: python script.py /path/to/scipion pwchem -noGPU"
	parser = argparse.ArgumentParser(
		epilog=epilog,
		formatter_class=argparse.RawDescriptionHelpFormatter
	)
	parser.add_argument("scipion", help="Path to Scipion executable, relative or absolute")
	parser.add_argument("plugin", help="Name of the plugin's Python module")
	parser.add_argument("-j", "--jobs", type=int, default=multiprocessing.cpu_count(), help="Number of jobs. Defaults to max available")
	parser.add_argument("-noGPU", action='store_true', help="If set, no tests that need a GPU will run. Use it in enviroments where a GPU cannot be accessed.")
	parser.add_argument("-testData", default='', help="Location of the test data JSON file.")
	args = parser.parse_args()

	# Set the test data file path
	if args.testData:
		args.testData = os.path.expanduser(args.testData)

	# Call main function
	main(args)
