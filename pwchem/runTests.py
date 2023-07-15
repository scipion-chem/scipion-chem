# General imports
import subprocess, argparse, multiprocessing, sys

# Function to colorize a string
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

# Function to run a test
def runTest(test):
    """ This function receives a test and runs it. """
    print(f"Running test {test}...")
    try:
        # Running test
        result = subprocess.run([args.scipion, testPrefix + test], check=True, capture_output=True, text=True)
        if result.returncode == 0:
            print(colorStr(f"Test {test} OK", color='green'))
    except subprocess.CalledProcessError as e:
        # Detect failed test
        print(e.stderr)
        print(colorStr(f"Test {test} failed with above message.", color='red'))
        return test

# Parse the command-line arguments
parser = argparse.ArgumentParser(
    epilog="Example: python script.py /path/to/scipion pwchem -j 2",
    formatter_class=argparse.RawDescriptionHelpFormatter
)
parser.add_argument("scipion", help="Path to Scipion executable, relative or absolute")
parser.add_argument("plugin", help="Name of the plugin's Python module")
parser.add_argument("-j", "--jobs", type=int, help="Number of jobs. Defaults to max available")
args = parser.parse_args()

# Set the number of jobs
args.jobs = args.jobs if args.jobs else multiprocessing.cpu_count()

# Construct the test discovery command
command = f"{args.scipion} test --grep {args.plugin}"

# Run the shell command and capture the output
try:
    output = subprocess.check_output(command, shell=True, text=True)
except subprocess.CalledProcessError:
    print(colorStr("ERROR: Test search command failed. Check line above for more detailed info.", color='red'))
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

# TEST:
#filteredLines = ['tests_databases.TestIdentifyLigands']

# If no tests were found, module was not found
if len(filteredLines) == 0:
    print(colorStr(f"ERROR: No tests were found for module {args.plugin}. Are you sure this module is properly installed?", color='red'))
    sys.exit(1)

# Showing initial message with number of tests
print(colorStr(f"Running a total of {len(filteredLines)} tests for {args.plugin} in batches of {args.jobs} processes...", color='yellow'))

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
    print(colorStr("Some tests ended with errors. Exiting.", color='red'))
    sys.exit(1)

# Message if all tests succeeded
print(colorStr("\nAll test passed!", color='green'))