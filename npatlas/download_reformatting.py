#!/usr/bin/env python3

import pandas as pd
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Draw, AllChem
import csv


def reformat_npclassifier_json(npclassifier_path):
    """Receive path to NP Atlas json. Transpose and export as csv, excluding fingerprints"""
    with open(npclassifier_path) as f:
        npclassifier_data = pd.read_json(f)

    transposed_npclassifier = npclassifier_data.transpose()

    full_df_output_path = Path(npclassifier_path.parent, str(npclassifier_path.stem) + "_transposed.csv")

    with open(full_df_output_path, "w") as g:
        transposed_npclassifier[["class_results",
                                 "superclass_results",
                                 "pathway_results",
                                 "isglycoside"]].to_csv(g)

    return transposed_npclassifier


def filter_npclassifier_json(transposed_json, classification_type, classification_term):
    """Filter npclassifier data by field and term"""

    superclass_expanded = transposed_json[classification_type].apply(pd.Series)
    macrolide_subset = superclass_expanded[superclass_expanded.isin([classification_term]).any(axis=1)]
    macrolide_npaids = macrolide_subset.index.to_list()
    macrolide_data = transposed_json.filter(items=macrolide_npaids, axis=0)

    filtered_df_output_path = Path(npclassifier_path.parent,
                                   str(npclassifier_path.stem) + "_" + classification_term + ".csv")
    with open(filtered_df_output_path, "w") as h:
        macrolide_data[["class_results",
                        "superclass_results",
                        "pathway_results",
                        "isglycoside"]].to_csv(h)

    return macrolide_npaids


def print_selected_structures(npaid_list, atlas_path, set_name):
    """Create pdf of structures from any list of NPAIDs"""
    structure_dict = {}

    with open(atlas_path) as f:
        csv_f = csv.reader(f, delimiter='\t')
        next(f)
        for row in csv_f:
            structure_dict[row[0]] = row[10]

    mol_list = []
    for npaid in npaid_list:
        mol = Chem.MolFromSmiles(structure_dict[npaid])
        mol.SetProp("_Name", npaid)
        AllChem.Compute2DCoords(mol)
        Draw.MolToFile(mol,
                       '/Users/roger/Git/rgl_utilities/data/NP_Atlas/Structure_images/' + npaid + '.png',
                       legend=npaid)

    img = Draw.MolsToGridImage(mol_list,
                             molsPerRow=4,
                             subImgSize=(200,200),
                             legends=[x.GetProp("_Name") for x in mol_list])


if __name__ == "__main__":

    npclassifier_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/np_atlas_2020_06_npclassifier.json")
    npatlas_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/np_atlas_2020_06.tsv")

    reformatted_npclassifier_df = reformat_npclassifier_json(npclassifier_path)
    filtered_napids = filter_npclassifier_json(reformatted_npclassifier_df, "superclass_results", "Macrolides")
    print_selected_structures(filtered_napids, npatlas_path, "atlas_macrolides")
