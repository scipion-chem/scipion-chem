# # -*- coding: utf-8 -*-
# # # **************************************************************************
# # # *
# # # * Authors: Alba Lomas Redondo (albalomasredon@gmail.com)
# # # *
# # # *
# * This program is free software; you can redistribute it and/or modify
# * it under the terms of the GNU General Public License as published by
# * the Free Software Foundation; either version 2 of the License, or
# * (at your option) any later version.
# *
# * This program is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# * GNU General Public License for more details.
# *
# * You should have received a copy of the GNU General Public License
# * along with this program; if not, write to the Free Software
# * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# * 02111-1307  USA
# *
# *  All comments concerning this program package may be sent to the
# *  e-mail address 'you@yourinstitution.email'
# *
# **************************************************************************
import sys, copy
import ast
import pickle
from rdkit import RDConfig, Chem, Geometry, DistanceGeometry
from rdkit.Chem import (
    ChemicalFeatures,
    rdDistGeom,
    Draw,
    rdMolTransforms,
    AllChem,
)
from rdkit.Chem.Pharm3D import Pharmacophore, EmbedLib
import os
import pandas as pd
import math

#estos paquetes creo que los tengo
import operator
import collections
from collections import OrderedDict, Counter
from sklearn import datasets, cluster

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
            if key == 'ligandFiles' or key == 'moleculesFiles':
                paramsDic[key] = value.strip().split()
            else:
                paramsDic[key] = value.strip()
    return paramsDic

##########################OBTENER LOS IDS DE LOS LIGANDOS######################

if __name__ == "__main__":
    '''Use: python <scriptName> <paramsFile> <outputDir>
    ParamsFile must include:
        <outputPath> <descritor> <receptorFile> <molFile1> <molFile2> ...'''
    paramsDic = parseParams(sys.argv[1])
    outDir = sys.argv[2]
    number = float(paramsDic['featuresNumber'])
    paths_ligands = paramsDic['ligandFiles']
    paths_molecules = paramsDic['moleculesFiles']
    #path_file = paramsDic['filePath']
    list_smiles = ast.literal_eval(paramsDic['ligandSmiles'])


#####################################################################
    ids = []
    molecules = []
    dict_molecules = {}

    for file in paths_ligands:
        # todo: more robust way to read id
        id = os.path.basename(file).split('.')[0]
        ids.append(id)
        #todo: read from other types of files
        molecule = Chem.MolFromPDBFile(file, removeHs=False)
        molecules.append(molecule)
        dict_molecules[id] = copy.copy(molecule)

    print('lenmol: ', len(molecules))
    with open(os.path.join(outDir, 'molecules.pkl'), 'wb') as outp:
        pickle.dump(molecules, outp, pickle.HIGHEST_PROTOCOL)

    print('Dict mol: ', dict_molecules)
    ################## ACONDICIONAR LA MOLÉCULA ########################
    # Load SMILES for PDB ligand structures
    #ligands = pd.read_csv(path_file)
    #dict_smiles = dict(zip(ligands["@structureId"], ligands["smiles"]))
    #print(dict_smiles)

    dict_smiles = {}
    for element in list_smiles:
        dict_smiles[element[0]] = element[1]
    print('Dic smiles. ', dict_smiles)

    dict_mol_reference = {}

    for key, value in dict_molecules.items():
        if key in dict_smiles.keys():
            smiles_mol = Chem.MolFromSmiles(dict_smiles[key])
            dict_mol_reference[key] = [dict_molecules[key], smiles_mol]

    print('Dict mol refe: ', dict_mol_reference)
    dict_mol_reference = OrderedDict(sorted(dict_mol_reference.items()))

    molecules_2 = []
    for key in dict_mol_reference:
        molecule, reference_molecule = dict_mol_reference[key]
        try:
            molecules_2 += [AllChem.AssignBondOrdersFromTemplate(reference_molecule, molecule)]
        except:
            print('AssignBondOrdersFromTemplate failed for ', key)

    ###################### RDKIT PHARMACOPHORE #######################
    fdefName = os.path.join(RDConfig.RDDataDir, 'BaseFeatures.fdef')
    factory = ChemicalFeatures.BuildFeatureFactory(fdefName)

    molecule_feature_frequencies = []
    for molecule2 in molecules_2:
        features = [feature.GetFamily() for feature in factory.GetFeaturesForMol(molecule2)]
        feature_frequency = collections.Counter(features)
        molecule_feature_frequencies.append(feature_frequency)

    feature_frequencies_df = (
        pd.DataFrame(molecule_feature_frequencies, index=[f"Mol{i}" for i, _ in enumerate(molecules_2, 1)],).fillna(0).astype(int))

    print(feature_frequencies_df.transpose())

    acceptors = []
    donors = []
    hydrophobics = []
    aromatics = []
    lumpedHydrophobe = []
    posIonizable = []

    for molecule_2 in molecules_2:
        acceptors.append(factory.GetFeaturesForMol(molecule_2, includeOnly="Acceptor"))
        donors.append(factory.GetFeaturesForMol(molecule_2, includeOnly="Donor"))
        hydrophobics.append(factory.GetFeaturesForMol(molecule_2, includeOnly="Hydrophobe"))
        aromatics.append(factory.GetFeaturesForMol(molecule_2, includeOnly="Aromatic"))
        lumpedHydrophobe.append(factory.GetFeaturesForMol(molecule_2, includeOnly="LumpedHydrophobe"))
        posIonizable.append(factory.GetFeaturesForMol(molecule_2, includeOnly="PosIonizable"))

    features = {
        "donor": donors,
        "acceptor": acceptors,
        "hydrophobic": hydrophobics,
        "posIonizable": posIonizable,
        "aromatic": aromatics,
        "lumpedHydrophobe": lumpedHydrophobe
    }

    #################################################################################
    features_coord = {}

    for key in features.keys():
        features_coord[key] = [list(item.GetPos()) for sublist in features[key] for item in sublist]

    # print(features_coord)

    #######HASTA AQUÍ TENEMOS LISTAS CON TODOS LAS FEATURES DE LOS TIPOS DESEADOS########
    # esto se tiene que definir como parametro

    kq = 7

    min_cluster_size = int(len(molecules) * 0.2)
    top_cluster_number = 4


    ############################################
    ##Se realiza el clustering con k-means

    def clustering(feature_coord, kq):
        k_means = None
        if feature_coord:
            k = math.ceil(len(feature_coord) / kq)
            if len(feature_coord) >= max([k, 2]):
                # Tailor-made adaption of k for hydrophobics in for the example in this talktorial
                k = 2 if k == 1 else k
                print(f"Clustering: \nVariable k in k-means: {k} of {len(feature_coord)} points\n")

                # Initialize k-means and compute clustering
                k_means = cluster.KMeans(n_clusters=k)
                k_means.fit(feature_coord)

        return k_means


    k_means = {}

    for coord in features_coord.keys():
        k_means[coord] = clustering(features_coord[coord], kq)


    ##Selecciona los clusters de una agrupación de k-means y devuelve los índices de los clusters seleccionados.

    def get_clusters(k_means, min_cluster_size, top_cluster_number):
        if k_means:

            feature_labels = k_means.labels_

            feature_labels_count = Counter(feature_labels)

            feature_labels_count = sorted(
                feature_labels_count.items(), key=operator.itemgetter(1), reverse=True
            )
            print(f"Clusters sorted by size: \n{feature_labels_count}\n")

            # Get number of the largest clusters, which are larger then the threshold (selected clusters)
            cluster_indices_sel = []

            for cluster_index, cluster_size in feature_labels_count:
                if cluster_size >= min_cluster_size and top_cluster_number > 0:
                    cluster_indices_sel.append(cluster_index)
                    top_cluster_number -= 1

            print(f"Cluster indices of selected clusters: \n{cluster_indices_sel}\n")

        else:
            cluster_indices_sel = None

        return cluster_indices_sel


    cluster_indices_sel_donors = get_clusters(k_means["donor"], min_cluster_size, top_cluster_number)
    cluster_indices_sel_acceptors = get_clusters(k_means["acceptor"], min_cluster_size, top_cluster_number)
    cluster_indices_sel_hydrophobic = get_clusters(k_means["hydrophobic"], min_cluster_size, top_cluster_number)
    cluster_indices_sel_posIonizable = get_clusters(k_means["posIonizable"], min_cluster_size, top_cluster_number)
    cluster_indices_sel_aromatics = get_clusters(k_means["aromatic"], min_cluster_size, top_cluster_number)
    cluster_indices_sel_lumpedHydrophobe = get_clusters(k_means["lumpedHydrophobe"], min_cluster_size, top_cluster_number)

    cluster_indices_sel = {
        "donor": cluster_indices_sel_donors,
        "acceptor": cluster_indices_sel_acceptors,
        "hydrophobic": cluster_indices_sel_hydrophobic,
        "posIonizable": cluster_indices_sel_posIonizable,
        "aromatic": cluster_indices_sel_aromatics,
        "lumpedHydrophobe": cluster_indices_sel_lumpedHydrophobe
        }

    with open(os.path.join(outDir, 'cluster_index.pkl'), 'wb') as outp:
        pickle.dump(cluster_indices_sel, outp, pickle.HIGHEST_PROTOCOL)

    # Recupera las coordenadas del centro del cluster para los clusters seleccionados.
    def get_selected_cluster_center_coords(k_means, cluster_indices_sel):
        cluster_centers_sel = " "
        if cluster_indices_sel:
            cluster_centers = k_means.cluster_centers_
            # Cast to list and then to pandas Series (for element selection by indices)
            cluster_centers = pd.Series(cluster_centers.tolist())
            # Select cluster centers by indices of selected clusters
            cluster_centers_sel = cluster_centers[cluster_indices_sel]

        else:
            pass

        return list(cluster_centers_sel)


    cluster_centers_sel = {}

    for center in cluster_indices_sel.keys():
        center2 = str(center).capitalize()
        cluster_centers_sel[center2] = get_selected_cluster_center_coords(k_means[center], cluster_indices_sel[center])

    # print(cluster_centers_sel)
    with open(os.path.join(outDir, 'cluster_centers.pkl'), 'wb') as outp:
        pickle.dump(cluster_centers_sel, outp, pickle.HIGHEST_PROTOCOL)

    final_dict = {}
    for key1, values in cluster_centers_sel.items():
        i = 0
        for value in values:
            if value == " ":
                pass
            else:
                i = i + 1
                key_new = key1 + str(i)
                final_dict[key_new] = tuple(value)

    print('Final dict: ', final_dict)
    ############################YA TENEMOS EL DICCIONARIO FINAL##################
    # ¿cuantas features se pueden comprobar?

    keys = final_dict.keys()
    features_pharmacophore = []
    for key in keys:
        features_pharmacophore.append(key[:-1])

    features_pharmacophore = list(dict.fromkeys(features_pharmacophore))
    print(features_pharmacophore)

    #########################################################################

    def create_Ph4Feats(final_dict, feature):
        samplePh4Feats = []
        for key, coord in final_dict.items():
            if key[:-1] == feature:
                print(feature)
                a = coord[0]
                b = coord[1]
                c = coord[2]

                samplePh4Feats.append(ChemicalFeatures.FreeChemicalFeature(str(key[:-1]), Geometry.Point3D(a, b, c)))
        pcophore = Pharmacophore.Pharmacophore(samplePh4Feats)

        return(pcophore)

    mols_molecules_dict, mols_molecules = preprocessLigands(paths_molecules)

    f = open(paramsDic['outputPath'], 'w')
    f.write("#The following molecules have been passed the pharmacophore filter:" + "\n")

    for key, item in mols_molecules_dict.items():
        contador = 0
        for e in features_pharmacophore:
            pcophore = create_Ph4Feats(final_dict, e)
            canMatch, allMatches = EmbedLib.MatchPharmacophoreToMol(key, factory, pcophore)
            print(canMatch)
            if canMatch == True:
                contador += 1

        if contador >= number:
            f.write(item + "\n")

