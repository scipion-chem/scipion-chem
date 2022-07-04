import math
import sys
import pandas as pd
from chembl_webresource_client.new_client import new_client

#############################################################################################################

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


#################Ad bioactivity values #############
def convert_ic50_to_pic50(IC50_value):
    pIC50_value = 9 - math.log10(IC50_value)
    return pIC50_value

def convert_ec50_to_pec50(EC50_value):
    pEC50_value = 9 - math.log10(EC50_value)
    return pEC50_value


#############################################################################################################
if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    uniprot_id = str(paramsDic['uniprotID'])
    target_type = str(paramsDic['targetType'])
    assay_type = str(paramsDic['assayType'])
    relation = str(paramsDic['relation'])
    bioactivity = str(paramsDic['bioactivity'])
    structures_number = int(paramsDic['numStructures'])

    targets_api = new_client.target
    compounds_api = new_client.molecule
    bioactivities_api = new_client.activity

    targets = targets_api.get(target_components__accession=uniprot_id).only("target_chembl_id", "organism", "pref_name", "target_type")

    targets = pd.DataFrame.from_records(targets)

    #############COGER EL PRIMER RESULTADO DE LA BÚSQUEDA################

    list_indice = []

    for indice_fila, fila in targets.iterrows():
        if fila[3] == target_type:
            list_indice.append(indice_fila)

    if len(list_indice) > 1:
        target = targets.iloc[int(list_indice[0])]

    elif len(list_indice) == 0:
        print("Cannot continue as the target type is not found")
        sys.exit()

    else:
        target = targets.iloc[int(list_indice[0])]

        #####################################################################

    chembl_id = target.target_chembl_id

    ##############################Get bioactivity data########################################

    try:
        bioactivities = bioactivities_api.filter(target_chembl_id=chembl_id, type=bioactivity, relation=relation, assay_type=assay_type[1:2]).only(
        "activity_id",
        "assay_chembl_id",
        "assay_description",
        "assay_type",
        "molecule_chembl_id",
        "type",
        "standard_units",
        "relation",
        "standard_value",
        "target_chembl_id",
        "target_organism",)

    except:
        print("The ChEMBL query does not exist")
        sys.exit()

    bioactivities_df = pd.DataFrame.from_records(bioactivities)

    ###################################################### ARRGLAR COSITAS DEL DATAFRAME #################

    bioactivities_df["units"].unique()
    bioactivities_df.drop(["units", "value"], axis=1, inplace=True)
    bioactivities_df = bioactivities_df.astype({"standard_value": "float64"})

    bioactivities_df.dropna(axis=0, how="any", inplace=True)
    bioactivities_df['standard_units'].unique()
    bioactivities_df = bioactivities_df[bioactivities_df["standard_units"] == "nM"]

    bioactivities_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)
    bioactivities_df.reset_index(drop=True, inplace=True)
    bioactivities_df.rename(columns={"standard_value": bioactivity, "standard_units": "units"}, inplace=True)

    ######################################### GET compound data ##############################

    compounds_provider = compounds_api.filter(molecule_chembl_id__in=list(bioactivities_df["molecule_chembl_id"])).only(("molecule_chembl_id", "molecule_structures"))
    compounds = list(compounds_provider)
    compounds_df = pd.DataFrame.from_records(compounds,)

    ####################### PROCESS COMPOUNDS DF ##################

    compounds_df.dropna(axis=0, how="any", inplace=True)
    compounds_df.drop_duplicates("molecule_chembl_id", keep="first", inplace=True)

    canonical_smiles = []

    for i, compounds in compounds_df.iterrows():
        try:
            canonical_smiles.append(compounds["molecule_structures"]["canonical_smiles"])
        except KeyError:
            canonical_smiles.append(None)

    compounds_df["smiles"] = canonical_smiles
    compounds_df.drop("molecule_structures", axis=1, inplace=True)
    compounds_df.dropna(axis=0, how="any", inplace=True)

    ##############OBTAIN THE FINAL DATASET ########################
    # Merge DataFrames
    output_df = pd.merge(
        bioactivities_df[["molecule_chembl_id", bioactivity, "units"]],
        compounds_df,
        on="molecule_chembl_id",)

    # Reset row indices
    output_df.reset_index(drop=True, inplace=True)

    #################################################################################################
    if bioactivity == "IC50":
       output_df["bioactivity2"] = output_df.apply(lambda x: convert_ic50_to_pic50(x.IC50), axis=1)
       output_df = output_df.sort_values(by=['bioactivity2'], ascending=False)
       output_df.reset_index(drop=True, inplace=True)

    elif bioactivity == "EC50":
        output_df["bioactivity2"] = output_df.apply(lambda x: convert_ec50_to_pec50(x.EC50), axis=1)
        output_df = output_df.sort_values(by=['bioactivity2'], ascending=False)
        output_df.reset_index(drop=True, inplace=True)

    #elif bioactivity == "Ki":
        #output_df = output_df.sort_values(by=[bioactivity], ascending=True)
        #output_df.reset_index(drop=True, inplace=True)


    else:
        output_df.rename(columns = {"bioactivity", "bioactivity2"})
        output_df = output_df.sort_values(by=["bioactivity2"], ascending=True)
        output_df.reset_index(drop=True, inplace=True)

    #####Obtain final files

    output_df = output_df[0:structures_number]

    #######CSV FILE#########
    #output_df.to_csv("EGFR_compounds.csv")

    dict_smiles = {}

    for i in output_df.index:
        chembl_id = output_df["molecule_chembl_id"][i]
        smiles = output_df["smiles"][i]
        bioactivity_f = output_df["bioactivity2"][i]
        dict_smiles[chembl_id] = [smiles, bioactivity_f]

    f = open("output.smi", "w")
    for key, item in dict_smiles.items():
        f.write(item[0] + "," + key + "," + str(item[1]) + "\n")

    f.close()