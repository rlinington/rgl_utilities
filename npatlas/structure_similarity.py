#!/usr/bin/env python3

"""Tools to compare new molecules to the chemical space in the NP Atlas DB"""

import pandas as pd
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def generate_fingerprint(smiles):
    """Calculate molecular fingerprint using Morgan fingerprint radius 2"""

    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprint(mol, 2)


def calculate_similarity(atlas_fp, input_fp):
    """Calculate similarity score between two rdkit fingerprints using Dice similarity scoring"""

    return round(DataStructs.DiceSimilarity(atlas_fp, input_fp), 2)


def single_compound_structure_similarity(atlas_path: Path,
                                         compound_name: str,
                                         compound_smiles: str):
    """Calculate similarity scores between a test compound and all compounds in Atlas DB"""

    new_compound_fingerprint = generate_fingerprint(compound_smiles)

    input_df = pd.read_csv(atlas_tsv_path, delimiter='\t')
    input_df["molecule_fingerprint"] = input_df["compound_smiles"].map(generate_fingerprint)
    input_df["Dice_similarity_score"] = input_df["molecule_fingerprint"].apply(calculate_similarity,
                                                                               args=(new_compound_fingerprint,))
    input_df = input_df.drop(["molecule_fingerprint"], axis=1)

    output_path = Path(atlas_tsv_path.parent, compound_name + "_atlas_structure_similarity.tsv")
    input_df.to_csv(output_path, index=False, sep='\t', lineterminator='\n')


def multiple_compound_structure_similarity(ref_db: pd.DataFrame, test_compounds: list, analysis_name: str):
    """Create matrix containing similarity scores between reference DB (rows) and test compounds (columns)"""

    ref_db["molecule_fingerprint"] = ref_db["compound_smiles"].map(generate_fingerprint)

    for compound in test_compounds:
        new_compound_fingerprint = generate_fingerprint(compound[1])
        ref_db[compound[0]] = ref_db["molecule_fingerprint"].apply(calculate_similarity,
                                                                   args=(new_compound_fingerprint,))

    ref_db = ref_db.drop(["molecule_fingerprint"], axis=1)

    output_path = Path(atlas_tsv_path.parent, analysis_name + "_atlas_structure_similarities.tsv")
    ref_db.to_csv(output_path, index=False, sep='\t', lineterminator='\n')

    return ref_db


def year_range_similarity_matrix(atlas_db: pd.DataFrame, start_year: int, end_year: int):
    """Filter Atlas DB so that ref_db is everything before start year, and compound list is everything within start and
    end years (inclusive) then create matrix of similarity scores. Return dataframe and write to file"""
    filtered_db = atlas_db[atlas_db["original_reference_year"] < start_year]
    test_compounds = atlas_db[(atlas_db["original_reference_year"] >= start_year) &
                              (atlas_db["original_reference_year"] <= end_year)]
    test_compound_list = test_compounds[["npaid", "compound_smiles"]].values.tolist()

    similarity_matrix = multiple_compound_structure_similarity(filtered_db,
                                                               test_compound_list,
                                                               f"similarity_{start_year}_{end_year}")

    return similarity_matrix


def count_rare_compounds(similarity_matrix: pd.DataFrame, maximum_score: float):
    """Count the number of compounds with maximum similarity scores at or below set value. Print NPAIDs for all 'rare'
    compounds"""

    npa_matrix = similarity_matrix.filter(regex='NPA0')
    rare_compounds = npa_matrix[npa_matrix.columns[npa_matrix.max() <= maximum_score]]
    print(f"Rare compounds: {rare_compounds.columns}")
    print(f"Number of rare compounds with similarities <= {maximum_score}: {len(rare_compounds.columns)}")


def count_new_entries(old_db: pd.DataFrame, new_db: pd.DataFrame, entry_term: str):
    """Return count of the number of any field added to the DB between versions
    (e.g. how many genera are only in the new version of the DB?)"""

    new_genus_list = set(new_db[entry_term].values.tolist()).difference(set(old_db[entry_term].values.tolist()))
    print(new_genus_list)
    print(f"Number of new {entry_term}: {len(new_genus_list)}")


def fungal_bacterial_ratio(atlas_db: pd.DataFrame, start_year: int, end_year: int):
    """Determine the ratio of fungal to bacterium compounds for a given time period
    Args:
        atlas_db (pd.Dataframe): A Pandas df of the Atlas download
        start_year (int): The starting year of the analysis
        end_year (int): the end year of the analysis (inclusive)"""

    compounds = atlas_db[atlas_db["original_reference_year"].between(start_year, end_year, inclusive='both')]
    bacterial_count = len(compounds[compounds["origin_type"] == "Bacterium"])
    fungal_count = len(compounds[compounds["origin_type"] == "Fungus"])

    print(f"Fungal:bacterial ratio = {round(fungal_count/bacterial_count, 2)}:1")


def count_new_articles(atlas_db: pd.DataFrame, start_year: int, end_year: int):
    """Count number of new isolation papers added to DB in a given period"""

    compounds = atlas_db[atlas_db["original_reference_year"].between(start_year, end_year, inclusive='both')]
    papers_count = len(compounds["original_reference_doi"].unique())

    print(f"Number of new isolation papers added: {papers_count}")


if __name__ == "__main__":

    # Similarity scoring for one compound
    # atlas_tsv_path = Path("/Users/roger/Git_new/rgl_utilities/data/NP_Atlas/NPAtlas_download_2023_06.tsv")
    # new_compound_name = "Megapolypeptin B"
    # new_compound_smiles = "CC(C)C(=O)C(=O)NC(CC(=O)NC(C(=O)NC(C(=O)NC(CO)C(O)CC(N)=O)C(C)O)C(C)O)CC(O)CC/C=C/CC/C=C/CC/C=C/CCCC(C)OC(CC(=O)CCC(=O)O)C(=O)O"
    # single_compound_structure_similarity(atlas_tsv_path, new_compound_name, new_compound_smiles)

    # Rare compound count for year range
    # atlas_tsv_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_03.tsv")
    # input_df = pd.read_csv(atlas_tsv_path, delimiter='\t')
    # similarity_matrix = year_range_similarity_matrix(input_df, 2021, 2022)
    # count_rare_compounds(similarity_matrix, 0.5)

    # Count of new genera since previous DB release
    old_database = pd.read_csv("/Users/roger/Git/rgl_utilities/data/NP_Atlas/np_atlas_2021_08.tsv", delimiter='\t')
    new_database = pd.read_csv("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.tsv", delimiter='\t')
    count_new_entries(old_database, new_database, "genus")
    count_new_entries(old_database, new_database, "original_reference_doi")

    # Determine ratio of fungal:bacterial compounds for a given timespan
    fungal_bacterial_ratio(new_database, 2019, 2020)

