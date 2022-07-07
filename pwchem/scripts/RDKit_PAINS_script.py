import sys
from rdkit import Chem
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

###Funcion para procesar los archivos de input

def preprocessLigands(ligandsFiles):
    mols_dict = {}

    for molFile in ligandsFiles:
        if molFile.endswith('.mol2'):
            m = Chem.MolFromMol2File(molFile)
            mols_dict[m] = molFile

        elif molFile.endswith('.mol'):
            m = Chem.MolFromMolFile(molFile)
            mols_dict[m] = molFile

        elif molFile.endswith('.pdb'):
            m = Chem.MolFromPDBFile(molFile)
            mols_dict[m] = molFile

        elif molFile.endswith('.smi'):
            f = open(molFile, "r")
            firstline = next(f)
            m = Chem.MolFromSmiles(str(firstline))
            mols_dict[m] = molFile

        elif molFile.endswith('.sdf'):
            suppl = Chem.SDMolSupplier(molFile)
            for mol in suppl:
                mols_dict[mol] = molFile

    mols = list(mols_dict.keys())

    return mols_dict, mols


def parseParams(paramsFile):
    paramsDic = {}
    with open(paramsFile) as f:
        for line in f:
            key, value = line.strip().split(':')
            if key == 'ligandFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic


#################################################################################################################
# initialize filter


def rdkit_filter(mols_dict):
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    matches = []
    clean = []
    for key, item in mols_dict.items():
        entry = catalog.GetFirstMatch(key)  # Get the first matching PAINS
        if entry is not None:
            # store PAINS information
            name = item.split("/")[-1]
            matches.append([item, entry.GetDescription().capitalize()])

        else:
            # collect indices of molecules without PAINS
            clean.append(item)

    return matches, clean

######################################################
def pains_filt(files, dic):
    mols_dict, mols = preprocessLigands(files)

    list_pains = []
    list_no_pains = []
    for mol, name in mols_dict.items():
        counter = 0
        for k, v in dic.items():
            subs = Chem.MolFromSmarts(k)
            if subs != None:
                if mol.HasSubstructMatch(subs):
                    counter +=1
                else:
                    pass
        if counter == 0:
            list_no_pains.append(name)
        else:
            list_pains.append(name)

    return set(list_pains), set(list_no_pains)



#################################################################################################################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    files = paramsDic['ligandFiles']

    if "painsFile" in paramsDic.keys():
        pains_file = str(paramsDic['painsFile'])

        inf = open(pains_file, "r")
        sub_strct = [line.rstrip().split(" ") for line in inf]
        smarts = [line[0] for line in sub_strct]
        desc = [line[1] for line in sub_strct]
        dic = dict(zip(smarts, desc))

        pains, no_pains = pains_filt(files, dic)

        with open(paramsDic['outputPath'], 'w') as f:

            f.write("#The following molecules have passed the PAINS filter:" + "\n")
            no_pains = list(no_pains)
            for molecule in no_pains:
                f.write(str(molecule) + "\n")

        with open("with_pains.txt", "w") as a:
            a.write("#Molecules that have NOT passed the filter " + "\n")
            pains = list(pains)
            for molecule_p in pains:
                a.write(str(molecule_p) + "\n")

        a.close()

    else:
        mols_dict, mols = preprocessLigands(files)
        matches, clean = rdkit_filter(mols_dict)

        #print(matches)

        with open(paramsDic['outputPath'], 'w') as f:
            f.write("#The following molecules have been passed the PAINS filter:" + "\n")
            for element in clean:
                f.write(element + "\n")

        f.close()

        with open("with_pains.txt", "w") as a:
            a.write("#Molecules that have NOT passed the filter " + "\n")
            for element1 in matches:
                a.write(element1[0] + "," + element1[1] + "\n")

        a.close()

#################################################################################################################

