import os, time, re
from subprocess import Popen, PIPE
from pwem import Plugin
from scipion.utils import getScipionHome

import pwchem

def clearPrintLines(n):
    for j in range(n):
        print(LINE_UP, end=LINE_CLEAR)


def parseTests(testFile):
    testNames, testName = [], None
    with open(testFile) as f:
        for line in f:
            if line.startswith('class '):
                testName = line.split()[1].split('(')[0]
            elif line.strip().startswith('def test') and testName:
                testNames.append(testName)
                testName = None
    return testNames

def parseOutputTests(outDir, packageName):
    failed, passed, totalTime = [], {}, 0
    with open(os.path.join(outDir, packageName + '.txt')) as f:
        for line in f:
            sline = line.strip()
            if re.findall('\[ *FAILED *\]', line) and re.findall('Test.*\.', line):
                failed.append(sline.split()[-1])
            elif re.findall('\[ *RUN *OK *\]', line):
                passed[sline.split()[-3]] = sline.split()[-2].replace('(', '')
            elif re.findall('\[=*\]', line):
                totalTime += float(sline.split()[-2].replace('(', ''))
    return {'PASSED': passed, 'FAILED': failed, 'TIME': totalTime}

if __name__ == "__main__":
    fnDir = os.path.split(pwchem.__file__)[0]
    outDir = os.path.join(fnDir, 'allTests')
    if not os.path.exists(outDir):
        os.mkdir(outDir)

    scipion3 = os.path.join(getScipionHome(), 'scipion3')
    devnull = open(os.devnull, 'wb') #python >= 2.4
    ps = {}
    for test_file in os.listdir('tests'):
        if not test_file.startswith('__'):
            packageName = test_file.replace('.py', '')
            testNames = parseTests('tests/' + test_file)

            commands = ['{} conda activate scipion3'.format(Plugin.getCondaActivationCmd())]
            for tName in testNames:
                commands += ['{} tests pwchem.tests.{}.{}'.format(scipion3, packageName, tName)]
            command = ' ; '.join(commands)

            with open("{}/{}.txt".format(outDir, packageName), "wb") as err:
                p = Popen(command, shell=True, stderr=err, stdout=devnull)

            ps[packageName] = p

    finished, first = [], True
    keepCheck = True
    LINE_UP = '\033[1A'
    LINE_CLEAR = '\x1b[2K'

    print('STATUS:')
    while keepCheck:
        if not first:
            clearPrintLines(len(ps))
        else:
            first = False

        for packageName in ps:
            status = '\033[93m Running \033[0m'
            if ps[packageName].poll() != None:
                status = '\033[92m Finished \033[0m'
                if not packageName in finished:
                    finished.append(packageName)
            print('Tests {} : {}'.format(packageName, status))

        if len(finished) == len(ps):
            keepCheck = False
        else:
            time.sleep(10)

    clearPrintLines(len(ps))
    for packageName in ps:
        print('Tests {} : \033[92m Finished \033[0m'.format(packageName))

    print('All tests finished')

    print('\nSUMMARY')
    sumDic = {}
    for pName in ps:
        sumDic[pName] = parseOutputTests(outDir, pName)

        nPassed = len(sumDic[pName]['PASSED'])
        nTotal = nPassed + len(sumDic[pName]['FAILED'])
        ttime = round(sumDic[pName]['TIME'], 2)

        sumStr = '\nTests {}: \n\tPassed: {} / {}\tTime: {} secs'.\
            format(pName, nPassed, nTotal, ttime)
        if nPassed != nTotal:
            sumStr += '\n\tFailed: {}'.format(', '.join(sumDic[pName]['FAILED']))
        print(sumStr)