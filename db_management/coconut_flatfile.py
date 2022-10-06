#!/usr/bin/env python3

"""Convert COCONUT output from JvS external SD drive into flat file"""

import json
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors
import csv
import sys


def reformat_jvs_final_json():
    """Reformat Jeff's input JSON (from JVS SSD) to original SNAP-MS TSV format"""

    with open(Path("/Users/roger/Git/rgl_utilities/data/coconut/inserted_FINAL.json")) as f:
        coconut_input = json.load(f)

    print(len(coconut_input))
    print(coconut_input[0])

    column_headers = [['npaid',
                       'compound_accurate_mass',
                       'compound_m_plus_h',
                       'compound_m_plus_na',
                       'compound_smiles',
                       'compound_names',
                       'npatlas_url']]

    coconut_list = []

    for entry in coconut_input:
        mol = Chem.MolFromSmiles(entry["smiles"])
        exact_mass = round(Chem.Descriptors.ExactMolWt(mol), 4)
        m_plus_h = exact_mass + 1.0073
        m_plus_na = exact_mass + 22.9898
        if type(entry["name"]) == 'str':
            name = entry["name"]
        else:
            name = "Unknown"
        coconut_list.append(["NPU" + str(entry["npuaid"]).zfill(6),
                             exact_mass,
                             m_plus_h,
                             m_plus_na,
                             entry["smiles"],
                             name,
                             None])

    with open(Path("/Users/roger/Git/rgl_utilities/data/coconut/coconut_snapms.tsv"),
              "w",
              newline='',
              encoding='utf-8') as g:
        csv_g = csv.writer(g, delimiter='\t')
        csv_g.writerows(column_headers + coconut_list)


def reformat_coconut_sdf():
    """Reformat COCONUT combined sdf output into original SNAP-MS TSV format"""

    coconut_sdf = Chem.SDMolSupplier("/Users/roger/Git/rgl_utilities/data/coconut/COCONUT.sdf")

    column_headers = [['npaid',
                       'compound_accurate_mass',
                       'compound_m_plus_h',
                       'compound_m_plus_na',
                       'compound_smiles',
                       'compound_names',
                       'npatlas_url']]

    coconut_list = []
    error_count = 0

    for mol in coconut_sdf:
        if mol:
            coconut_id = mol.GetProp('coconut_id')
            exact_mass = round(Chem.Descriptors.ExactMolWt(mol), 4)
            m_plus_h = exact_mass + 1.0073
            m_plus_na = exact_mass + 22.9898
            name = mol.GetProp('name')
            insert_data = [coconut_id,
                           exact_mass,
                           m_plus_h,
                           m_plus_na,
                           mol.GetProp('SMILES'),
                           name,
                           "https://coconut.naturalproducts.net/compound/coconut_id/" + coconut_id]
            coconut_list.append(insert_data)

        else:
            error_count += 1

    print("Number of compounds with RDKit errors: " + str(error_count))

    with open(Path("/Users/roger/Git/rgl_utilities/data/coconut/coconut_sdf_snapms.tsv"),
              "w",
              newline='',
              encoding='utf-8') as g:
        csv_g = csv.writer(g, delimiter='\t')
        csv_g.writerows(column_headers + coconut_list)


def test_sdf():
    """Create example RDKit SDF file"""
    smiles_list = ["CCCC", "CCCCO"]

    writer = Chem.SDWriter('/Users/roger/Git/rgl_utilities/data/coconut/out.sdf')

    for smiles in smiles_list:
        writer.write(Chem.MolFromSmiles(smiles))


if __name__ == "__main__":

    reformat_coconut_sdf()