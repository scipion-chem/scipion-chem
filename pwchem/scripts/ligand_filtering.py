import os
import sys
from os import remove
import biotite.database.rcsb as rcsb
import pandas as pd
import pypdb
import requests
from datetime import datetime
from bs4 import BeautifulSoup
from opencadd.structure.core import Structure
from opencadd.structure.superposition.api import align, METHODS
from urllib.request import urlretrieve
import gzip
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import DBSCAN

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


#############################################################################################################


def describe_one_pdb_id(pdb_id):
    """Fetch meta information from PDB."""
    described = pypdb.describe_pdb(pdb_id)
    if described is None:
        print(f"! Error while fetching {pdb_id}, retrying ...")
        raise ValueError(f"Could not fetch PDB id {pdb_id}")
    return described


def get_ligands(pdb_id):
    info = pypdb.get_info(pdb_id)
    nonpolymers = info.get("rcsb_entry_info", {}).get("nonpolymer_bound_components", [])
    #print(nonpolymers)
    ligands = {}
    #Por ejemplo si tenemos como un nonpolymer 8AM, la direccion que buscamos es http://ligand-expo.rcsb.org/reports/8/8AM/
    for ligand_expo_id in nonpolymers:
    #buscamos la dirección
        r = requests.get(
            f"http://ligand-expo.rcsb.org/reports/{ligand_expo_id[0]}/{ligand_expo_id}/"
            )

        # Postprocess some known values
        r.raise_for_status()
        html = BeautifulSoup(r.text, features="html.parser")
        info = {}
        for table in html.find_all("table"):
            for row in table.find_all("tr"):
                cells = row.find_all("td")
                if len(cells) != 2:
                    continue
                key, value = cells
                if key.string and key.string.strip():
                    info[key.string.strip()] = "".join(value.find_all(text=True))

        # Postprocess some known values
        info["Molecular weight"] = float(info["Molecular weight"].split()[0])
        info["Formal charge"] = int(info["Formal charge"])
        info["Atom count"] = int(info["Atom count"])
        info["Chiral atom count"] = int(info["Chiral atom count"])
        info["Bond count"] = int(info["Bond count"])
        info["Aromatic bond count"] = int(info["Aromatic bond count"])
        ligands[ligand_expo_id] = info

    return ligands, nonpolymers



def obtain_pdb_gz(pdb_code):
    download_folder = '.'
    compressed = True
    filename = '%s.pdb' % pdb_code[:4]
    if compressed:
        filename = '%s.gz' % filename
        url = "http://files.rcsb.org/download/%s" % filename
        destination_file = os.path.join(download_folder, filename)
        urlretrieve(url, destination_file)


def obtain_mass_center(files, nonpolymer_dict):
    coord_list = []
    coord_x_list = []
    coord_y_list = []
    coord_z_list = []
    id_list = []

    for file1 in files:
        ligand = (nonpolymer_dict[file1[0:-7]][0])
        interest_lines = []
        f = gzip.open(file1, 'rb')

        for line in f:
            line = line.decode().strip()

            if ligand in line and "HETATM" in line:
                interest_lines.append(line.rsplit())

        sum_x = 0
        sum_y = 0
        sum_z = 0
        suma_masas = 0

        for element in interest_lines:
            x = float(element[5]) * float(masas[element[-1]])
            y = float(element[6]) * float(masas[element[-1]])
            z = float(element[7]) * float(masas[element[-1]])

            suma_masas += float(masas[element[-1]])
            sum_x += x
            sum_y += y
            sum_z += z

            coor_x = x/suma_masas
            coor_y = y/suma_masas
            coor_z = z/suma_masas


        coord_list.append([coor_x, coor_y, coor_z])
        coord_x_list.append(coor_x)
        coord_y_list.append(coor_y)
        coord_z_list.append(coor_z)
        id_list.append(file1[0:-7])

    return coord_list, coord_x_list, coord_y_list, coord_z_list, id_list



if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    uniprot_id = str(paramsDic['uniprotID'])
    experimental_method = str(paramsDic['experimentalMethod'])
    molecular_weight = float(paramsDic['molecularWeigth'])
    max_resolution = float(paramsDic['maxResolution'])
    polymer_count = int(paramsDic['polymerCount'])
    top_structures_number = int(paramsDic['topStructures'])
    before_deposition_date = str(paramsDic['date'])
    date_time_str = "01/01/" + str(before_deposition_date)
    date = datetime.strptime(date_time_str, '%d/%m/%Y')
    eps = float(paramsDic['eps'])
    min_samples = int(paramsDic['min_samples'])

    # Filter queries
    query_by_uniprot_id = rcsb.FieldQuery("rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession", exact_match=uniprot_id, )
    query_by_deposition_date = rcsb.FieldQuery("rcsb_accession_info.deposit_date", less=date)
    query_by_experimental_method = rcsb.FieldQuery("exptl.method", exact_match=experimental_method)
    query_by_resolution = rcsb.FieldQuery("rcsb_entry_info.resolution_combined", less_or_equal=max_resolution)
    query_by_polymer_count = rcsb.FieldQuery("rcsb_entry_info.deposited_polymer_entity_instance_count", equals=polymer_count)
    query_by_ligand_mw = rcsb.FieldQuery("chem_comp.formula_weight", molecular_definition=True, greater=molecular_weight)

    ###Total query

    if experimental_method == "None":
        query = rcsb.CompositeQuery(
        [query_by_uniprot_id, query_by_deposition_date, query_by_resolution,
        query_by_polymer_count, query_by_ligand_mw, ], "and", )


    else:
        query = rcsb.CompositeQuery(
        [query_by_uniprot_id, query_by_deposition_date, query_by_experimental_method, query_by_resolution,
        query_by_polymer_count, query_by_ligand_mw, ], "and", )

    # Do the query
    pdb_ids = rcsb.search(query)

    pdbs_data = [describe_one_pdb_id(pdb_id) for pdb_id in pdb_ids]

    # Obtain resolution dataframe
    resolution = pd.DataFrame(
        [[pdb_data["entry"]["id"], pdb_data["rcsb_entry_info"]["resolution_combined"][0]] for pdb_data in pdbs_data],
        columns=["pdb_id", "resolution"]).sort_values(by="resolution", ignore_index=True)
    selected_pdb_ids = resolution[:top_structures_number]["pdb_id"].to_list()
    ############################################################################################################################

    nonpolymer_dict = {}
    ligands_dict = {}
    for pdb_id in selected_pdb_ids:
        ligands, nonpolymers = get_ligands(pdb_id)
        ligands_dict[pdb_id] = ligands
        nonpolymer_dict[pdb_id] = nonpolymers

    ###########Tienen las estructuras ligando??##########################

    for element in selected_pdb_ids:
        if not nonpolymer_dict[element]:
            selected_pdb_ids.remove(element)

        else:
            pass

    ###################################################################
    masas = {'H': 1.007825, 'C': 12.01, 'O': 15.9994, 'N': 14.0067, 'S': 31.972071,
             'P': 30.973762, "F": 18.998403}

    for pdb_code in selected_pdb_ids:
        obtain_pdb_gz(pdb_code)

    files = []

    for file in os.listdir("."):
        if file.endswith(".pdb.gz"):
            files.append(file)

    coord_list, coord_x_list, coord_y_list, coord_z_list, id_list = obtain_mass_center(files, nonpolymer_dict)

    ####################################CLUSTERING###############################
    data_df = {'IDs': id_list, 'X': coord_x_list, "Y": coord_y_list, "Z": coord_z_list}
    df_off = pd.DataFrame(data_df)
    print(df_off)

    data = np.array(coord_list)
    model = DBSCAN(eps=eps, min_samples=min_samples)
    model.fit_predict(data)
    pred = model.fit_predict(data)

    print("number of cluster found: {}".format(len(set(model.labels_))))
    print('cluster for each point: ', model.labels_)

    # clustering = pd.DataFrame({'IDs':df_off['IDs'],'X':df_off['X'], "Y":df_off['Y'], "Z":df_off['Z'], 'Grupo': pred})
    clustering = pd.DataFrame({'IDs': df_off['IDs'], 'Grupo': pred})
    print(clustering)
    dict_clustering = dict(clustering.values)
    print(dict_clustering)

    new = {}
    for key, value in dict_clustering.items():
        if value in new.keys():
            new[value].append(key)

        else:
            new[value] = [key]
    #########################################################################

    ############################################################################
    columns = [
        "@structureId",
        "@chemicalID",
        "@type",
        "@molecularWeight",
        "chemicalName",
        "formula",
        "InChI",
        "InChIKey",
        "smiles",
    ]
    rows = []
    res = []

    for pdb_id, ligands in ligands_dict.items():
        if len(ligands.keys()) != 0:
            ligand_id, properties = max(ligands.items(), key=lambda kv: kv[1]["Molecular weight"])
            rows.append(
                [
                    pdb_id,
                    ligand_id,
                    properties["Component type"],
                    properties["Molecular weight"],
                    properties["Name"],
                    properties["Formula"],
                    properties["InChI descriptor"],
                    properties["InChIKey descriptor"],
                    properties["Stereo SMILES (OpenEye)"],
                ]
            )

        else:
            pass

    ligands = pd.DataFrame(rows, columns=columns)

    ############# HASTA AQUI PUEDO HACER UN CSV CON ESAS CARACTERÍSTICAS PARA LOS LIGANDOS DE INTERÉS, WIIIIIIII #######################
    ligands.to_csv('PDB_top_ligands.csv', header=True, index=False)
    ###ESTO FUNCIONA PERO MIRAR BIEN LO QUE ES POR SI PODEMOS QUITARLO
    pairs = dict(zip(ligands["@structureId"], ligands["@chemicalID"]))
    print(pairs)

    #########################################################################################
    structures = [Structure.from_pdbid(pdb_id) for pdb_id in pairs]

    #######EXTRACT THE PROTEIN AND THE LIGAND #################
    complexes = [Structure.from_atomgroup(structure.select_atoms(f"protein or resname {ligand}")) for structure, ligand
                 in zip(structures, pairs.values())]

    ####### EXTRACT THE LIGANDS ##############
    ligands = [Structure.from_atomgroup(complex_.select_atoms(f"resname {ligand}")) for complex_, ligand in
               zip(complexes, pairs.values())]
    print(ligands)

    ####OUTPUT --> Ligandos en archivos pdb

    pdb_ligands = []
    for ligand, pdb_id in zip(ligands, pairs.keys()):
        name = f"{pdb_id}_lig.pdb"
        ligand.write(name)
        pdb_ligands.append(name)

    ##################OBTENER SOLO LA PROTEINA DE CADA############
    #pdb_proteins = []
    #for id_ in pairs:
        #out = str(id_) + "_protein.pdb"
        #out_protein = clean_PDB(id_, out, waters=True, HETATM=True)
        #remove(str(id_).lower() + ".pdb.gz")
        #remove(str(id_) + ".pdb.gz")
        #pdb_proteins.append(out_protein)

    ############################################################
    diccionario_output = {}

    for key1, values2 in new.items():
        lista = []
        for value1 in values2:
            lista.append(os.path.abspath(str(value1) + "_lig.pdb"))

        diccionario_output[key1] = lista
    ############################################################

    f = open("output.txt", "w")
    for key, values in diccionario_output.items():
        f.write("Cluster " + str(key) + '\n')
        for value in values:
            f.write(str(value) + '\n')
        f.write("#" + '\n')


    f.close()

    f = gzip.open(files[0], 'r')
    file_content = f.read()
    file_content = file_content.decode('utf-8')
    f_out = open('protein.pdb', 'w+')
    f_out.write(file_content)
    f.close()
    f_out.close()

    cwd = os.getcwd()
    f2 = open("protein.txt", "w")
    f2.write(str(cwd) + "/" + 'protein.pdb')
    f2.close()

    for id_ in pairs:
        #remove(str(id_).lower() + ".pdb.gz")
        remove(str(id_) + ".pdb.gz")





