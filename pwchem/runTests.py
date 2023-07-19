# General imports
import subprocess, argparse, multiprocessing, sys, json, os

################################## UTILS FUNCTIONS ##################################
def colorStr(string, color):
	""" This function returns the input string wrapped in the specified color. """
	if color == 'green':
		return f"\033[92m{string}\033[0m"
	elif color == 'yellow':
		return f"\033[93m{string}\033[0m"
	elif color == 'red':
		return f"\033[91m{string}\033[0m"
	else:
		return string

def runInParallel(func, *args, paramList, jobs):
	"""
	This function creates a pool of workers to run the given function in parallel.
	Also returns a list with the failed commands.
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

	# Return list of failed commands
	return failedCommands

def runTest(test, scipionExecutable, testPrefix):
	""" This function receives a test and runs it. """
	printAndFlush(f"Running test {test}...")
	try:
		# Running test
		result = subprocess.run([scipionExecutable, testPrefix + test], check=True, capture_output=True, text=True)
		if result.returncode == 0:
			printAndFlush(colorStr(f"Test {test} OK", color='green'))
	except subprocess.CalledProcessError as e:
		# Detect failed test
		printAndFlush(e.stderr)
		printAndFlush(colorStr(f"Test {test} failed with above message.", color='red'))
		return test

def printAndFlush(message):
	""" This function prints a given message and inmediately flushes stdout. """
	print(message)
	sys.stdout.flush()

def printSkippingTest(test, skipType='', dependency='', reason=''):
	"""
	This function prints the message when a test is skipped.
	Note: param dependency is only used when skipType is 'dependency',
	and param reason is only used when skipType is 'other'.
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

def printFatalError(message):
	""" This function prints the given message in red, flushes stdout, and exits with code 1. """
	printAndFlush(colorStr(message, color='red'))
	sys.exit(1)

def getCondaActivationCmd(scipionExecutable):
	""" This function reads the scipion config file and returns the conda activation command. """
	# Path to scipion config file
	configFile = os.path.join(os.path.dirname(scipionExecutable), 'config', 'scipion.conf')

	# Read file and return conda activation command
	with open(configFile, 'r') as f:
		for line in f:
			if line.startswith('CONDA_ACTIVATION_CMD'):
				return line.split('=')[1].strip()

################################## MAIN EXECUTION FUNCTIONS ##################################
def getAllTests(scipion, pluginModule, testPrefix):
	""" This function finds the full list of tests from a given module. """
	# Construct the test discovery command
	command = f"{scipion} test --grep {pluginModule}"

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
	
	# If no tests were found, module was not found
	if not filteredLines:
		printFatalError(f"ERROR: No tests were found for module {args.plugin}. Are you sure this module is properly installed?")
	
	# Return full list of tests
	return filteredLines

def readTestDataFile(testDataFilePath):
	"""
	This function returns a list with the necessary datasets for the tests, as well as
	an object with the different tests and the situations where to skip them.
	"""
	# Read the JSON data from the file
	try:
		with open(testDataFilePath, 'r') as file:
			dataFile = json.load(file)
			return dataFile.get("datasets", []), dataFile.get("skippable", {})
	except FileNotFoundError:
		printAndFlush(colorStr("No skippable tests file found, running all.", color='yellow'))
		return None, None
	except PermissionError:
		printFatalError(f"ERROR: Permission denied to open file '{testDataFilePath}'.")
	except json.JSONDecodeError as e:
		printFatalError(f"ERROR: Invalid JSON format in file '{testDataFilePath}':\n{e}")
	except Exception as e:
		printFatalError(f"An unexpected error occurred:\n{e}")
	
def downloadDatset(dataset, scipionExecutable):
	""" This function downloads a given dataset for scipion tests. """
	printAndFlush(colorStr(f"Downloading dataset {dataset}...", color='yellow'))
	try:
		result = subprocess.run([scipionExecutable, f" testdata --download {dataset}"], check=True, capture_output=True, text=True)
		if result.returncode == 0:
				printAndFlush(colorStr(f"Dataset {dataset} download OK", color='green'))
	except subprocess.CalledProcessError as e:
		# Detect failed test
		printAndFlush(e.stderr)
		printAndFlush(colorStr(f"Dataset {dataset} download failed with the above message.", color='red'))
		return dataset
	
def removeGPUTests(testList, noGPU, gpuSkippableTests):
	""" This function removes the GPU tests if flag noGPU has been set. """
	# Check if noGPU flag was set, only remove if so
	if noGPU:
		for gpuTest in gpuSkippableTests:
			# For each test, remove from list if it was there (dev might have mispelled the test on the .json file)
			if gpuTest in testList:
				printSkippingTest(gpuTest, skipType='gpu')
				testList.remove(gpuTest)
	
	# Return modified list of tests
	return testList

def removeDependecyTests(testList, dependenciesSkippableTests):
	""" This function removes the depencendy related tests if the dependencies are not met. """
	# Removing dependency tests if dependencies are not present
	for dependency in dependenciesSkippableTests:
		# Getting plugin and module name
		dependencyTestKey = dependency.get("name", None)
		dependencyTestModule = dependency.get("module", None)

		# Try to import module if provided
		if dependencyTestModule:
			# Creating import command to run within scipion3 conda env
			command = f"{getCondaActivationCmd(args.scipion)} && conda activate scipion3 && python -c 'import {dependencyTestModule}' 2>/dev/null && echo 1 || echo 0"
			sucess = int(subprocess.check_output(command, shell=True).decode().replace('\n', ''))
			if sucess:
				continue

		# If no module was provided or import raised a ModuleNotFoundError exception, skip tests
		for dependencyTest in dependency.get("tests", []):
			if dependencyTest in testList:
				printSkippingTest(dependencyTest, skipType='dependency', dependency=dependencyTestKey)
				testList.remove(dependencyTest)
	
	# Return modified list of tests
	return testList

def removeOtherTests(testList, otherSkippableTests):
	""" This function removes the other skippable tests. """
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

def removeSkippableTests(testList, noGPU, gpuSkippableTests, dependenciesSkippableTests, otherSkippableTests):
	""" This function removes from the list of all tests, the ones that have to be skipped. """
	# Remove GPU skippable tests
	testList = removeGPUTests(testList, noGPU, gpuSkippableTests)

	# Removing dependency tests if dependencies are not present
	testList = removeDependecyTests(testList, dependenciesSkippableTests)

	# Removing other tests for reasons stated
	testList = removeOtherTests(testList, otherSkippableTests)
	
	# Return new list of tests
	return testList

def getResultDictionary(testList, failedTests):
	""" This function returns a dictionary with all the passed/failed tests grouped by their origin file. """
	# Initialize an empty dictionary to store the grouped tests
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
	
	# Return result dictionary
	return results

def printSummary(results):
	""" This function prints the sumarry of the test results. """
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
def main(args):
	""" Main execution function. """
	# Defining test prefix
	testPrefix = f'tests {args.plugin}.tests.'

	# Getting full list of tests
	filteredLines = getAllTests(args.scipion, args.plugin, testPrefix)

	# Obtaining datasets and skippable tests according to situation
	datasets, allSkippableTests = readTestDataFile(args.testData)

	if allSkippableTests:
		gpuSkippableTests = allSkippableTests.get('gpu', [])
		dependenciesSkippableTests = allSkippableTests.get('dependencies', [])
		otherSkippableTests = allSkippableTests.get('others', [])

		# Removing skippable tests
		filteredLines = removeSkippableTests(filteredLines, args.noGPU, gpuSkippableTests, dependenciesSkippableTests, otherSkippableTests)

	# Downloading in parallel required datasets if there are any
	if datasets:
		failedDownloads = runInParallel(downloadDatset, args.scipion, paramList=datasets, jobs=args.jobs)

		# Check if there were any errors
		if failedDownloads:
			printFatalError("The download of at least one dataset ended with errors. Exiting.")

	# Showing initial message with number of tests
	nJobs = len(filteredLines) if len(filteredLines) < args.jobs else args.jobs
	printAndFlush(colorStr(f"Running a total of {len(filteredLines)} tests for {args.plugin} in batches of {nJobs} processes...", color='yellow'))

	# Run all the tests in parallel
	failedTests = runInParallel(runTest, args.scipion, testPrefix, paramList=filteredLines, jobs=args.jobs)

	# Get results grouped by orogin file and separated into passed and failed
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
	parser.add_argument("-j", "--jobs", type=int, help="Number of jobs. Defaults to max available")
	parser.add_argument("-noGPU", action='store_true', help="If set, no tests that need a GPU will run. Use it in enviroments where a GPU cannot be accessed.")
	parser.add_argument("-testData", help="Location of the test data JSON file.")
	args = parser.parse_args()

	# Set the number of jobs
	args.jobs = args.jobs if args.jobs else multiprocessing.cpu_count()

	# Set the test data file path
	args.testData = os.path.expanduser(args.testData) if args.testData else ''

	# Call main function
	main(args)