# General imports
import subprocess, argparse, multiprocessing, sys, json, os, importlib

# Global variables
skippableTestFileName = 'skippableTests.json'

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

def runTest(test):
	""" This function receives a test and runs it. """
	printAndFlush(f"Running test {test}...")
	try:
		# Running test
		result = subprocess.run([args.scipion, testPrefix + test], check=True, capture_output=True, text=True)
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

def readSkippableTestsFile():
	""" This function returns an object with the different tests and the situations where to skip them. """
	# Get file's location
	currentLocation = os.path.dirname(os.path.abspath(__file__))
	skippableTestFile = os.path.join(currentLocation, skippableTestFileName)

	# Read the JSON data from the file
	try:
		with open(skippableTestFile, 'r') as file:
			return json.load(file)
	except FileNotFoundError:
		printAndFlush(colorStr("No skippable tests file found, running all.", color='yellow'))
	except PermissionError:
		printAndFlush(colorStr(f"Error: Permission denied to open file '{skippableTestFile}'.", color='red'))
		sys.exit(1)
	except json.JSONDecodeError as e:
		printAndFlush(colorStr(f"Error: Invalid JSON format in file '{skippableTestFile}':\n{e}", color='red'))
		sys.exit(1)
	except Exception as e:
		printAndFlush(colorStr(f"An unexpected error occurred:\n{e}", color='red'))
		sys.exit(1)

# Parse the command-line arguments
parser = argparse.ArgumentParser(
	epilog="Example: python script.py /path/to/scipion pwchem -j 2",
	formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("scipion", help="Path to Scipion executable, relative or absolute")
parser.add_argument("plugin", help="Name of the plugin's Python module")
parser.add_argument("-j", "--jobs", type=int, help="Number of jobs. Defaults to max available")
parser.add_argument("-noGPU", action='store_true', help="If set, no tests that need a GPU will run. Use it in enviroments where a GPU cannot be accessed.")
args = parser.parse_args()

# Set the number of jobs
args.jobs = args.jobs if args.jobs else multiprocessing.cpu_count()

# Construct the test discovery command
command = f"{args.scipion} test --grep {args.plugin}"

# Run the shell command and capture the output
try:
	output = subprocess.check_output(command, shell=True, text=True)
except subprocess.CalledProcessError:
	printAndFlush(colorStr("ERROR: Test search command failed. Check line above for more detailed info.", color='red'))
	sys.exit(1)

# Define test command string variables
scipionTestsStartingSpaces = '   '
testPrefix = f'tests {args.plugin}.tests.'

# Separate lines into a list
lines = output.split('\n')

# For each line, keep only the ones about individual tests in a reduced form
filteredLines = []
for line in lines:
	if line.startswith(scipionTestsStartingSpaces):
		filteredLines.append(line.replace(f'{scipionTestsStartingSpaces}scipion3 {testPrefix}', ''))

# If no tests were found, module was not found
if len(filteredLines) == 0:
	printAndFlush(colorStr(f"ERROR: No tests were found for module {args.plugin}. Are you sure this module is properly installed?", color='red'))
	sys.exit(1)

# Obtaining skippable tests according to situation
allSkippableTests = readSkippableTestsFile()
gpuSkippableTests = allSkippableTests.get('gpu', [])
dependenciesSkippableTests = allSkippableTests.get('dependencies', [])
otherSkippableTests = allSkippableTests.get('others', [])
#print(filteredLines)

# Removing GPU skippable tests if no GPU flag provided by args
if args.noGPU:
	for gpuTest in gpuSkippableTests:
		if gpuTest in filteredLines:
			printAndFlush(colorStr(f"Skipping test {gpuTest}. Reason: GPU.", color='yellow'))
			filteredLines.remove(gpuTest)

# Removing dependency tests if dependencies are not present
for dependency in dependenciesSkippableTests:
	# Getting plugin and module name
	dependencyTestKey = dependency.get("name", None)
	dependencyTestModule = dependency.get("module", None)

	# Try to import module if provided
	try:
		if dependencyTestModule:
			importlib.import_module(dependencyTestModule)
			continue
	except ModuleNotFoundError:
		pass

	# If no module was provided or import raised a ModuleNotFoundError exception, skip tests
	for dependencyTest in dependency.get("tests", []):
		if dependencyTest in filteredLines:
			printAndFlush(colorStr(f"Skipping test {dependencyTest}. Reason: Unmet dependency with plugin {dependencyTestKey}.", color='yellow'))
			filteredLines.remove(dependencyTest)

# Removing other tests for reasons stated
for otherTest in otherSkippableTests:
	# Getting test name and reason
	testName = otherTest.get("test", None)
	testReason = otherTest.get("reason", None)

	# If test exists, skip
	if testName and testName in filteredLines:
		reason = f"Reason: {testReason}" if testReason else "No reason provided"
		printAndFlush(colorStr(f"Skipping test {testName}. {reason}.", color='yellow'))
		filteredLines.remove(testName)

# Showing initial message with number of tests
printAndFlush(colorStr(f"Running a total of {len(filteredLines)} tests for {args.plugin} in batches of {args.jobs} processes...", color='yellow'))

# Create a shared flag to indicate error occurrence
manager = multiprocessing.Manager()
errorFlag = manager.Value('i', 0)

# Create a pool of worker processes
pool = multiprocessing.Pool(processes=args.jobs)

# Apply the runTest function to the filtered lines using the pool
results = [pool.apply_async(runTest, args=(test,)) for test in filteredLines]

# Check if any process encountered an error
errorFlag = False
for result in results:
	if result.get():
		errorFlag = True

# Close the pool to release resources
pool.close()

# Join the pool, waiting for all processes to complete
pool.join()

# Check if an error occurred
if errorFlag:
	printAndFlush(colorStr("Some tests ended with errors. Exiting.", color='red'))
	sys.exit(1)

# Message if all tests succeeded
printAndFlush(colorStr("\nAll test passed!", color='green'))