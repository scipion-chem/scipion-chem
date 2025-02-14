import sys, csv, json
from ranx import Run, fuse

from utils import parseParams

def writeResults(combinedRun, outputPath):

    with open(outputPath, mode='w', newline='', encoding='utf-8') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        header = ['Rank', 'Mut', 'Score']
        writer.writerow(header)
        combinedScores = combinedRun.get_doc_ids_and_scores()

        for i, mutation in enumerate(combinedScores):
            for mutName, scoreComb in mutation.items():
                scoreComb = round(float(scoreComb), 4)
                writer.writerow([i+1] + [mutName] + [scoreComb])


#################################################################################################################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <normMethod> <fusionMethod> <kwargs> <runDics>'''
    paramsDic = parseParams(sys.argv[1], listParams=['ligandFiles'], sep='::')
    runDics = json.loads(paramsDic['runDics'])
    
    lruns = []
    for nameRun, runDic in runDics.items():
        lruns.append(Run({'q_1': runDic}, name=nameRun))

    combinedRun = fuse(runs=lruns,
                       norm=paramsDic['normMethod'],
                       method=paramsDic['fusionMethod'],
                       params=eval(paramsDic['kwargs']))
    
    writeResults(combinedRun, paramsDic['outputPath'])
