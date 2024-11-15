#! usr/bin/env python3

"""Tools to explore Atlas molecules based on NP Classifier designations"""

import pandas as pd
from pathlib import Path

ATLAS_JSON_PATH = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.json")


def extract_np_classifier_data(atlas_df):
    """Extract glycoside assignment from npclassifier results"""

    extracted_np_classifier_df = pd.json_normalize(atlas_df["npclassifier"])

    merged_df = atlas_df.join(extracted_np_classifier_df)
    merged_df["class_results"] = merged_df["class_results"].str.get(0)
    merged_df["pathway_results"] = merged_df["pathway_results"].str.get(0)
    merged_df["superclass_results"] = merged_df["superclass_results"].str.get(0)

    return merged_df


def convert_json(json_path):
    """Open json file as Pandas df

    Input: path to json
    Returns: pandas df"""

    with open(json_path, "rb") as f:
        atlas_df = pd.read_json(f)

    expanded_atlas_df = extract_np_classifier_data(atlas_df)
    expanded_atlas_df.rename(columns={'smiles': 'compound_smiles'}, inplace=True)

    return expanded_atlas_df


def extract_biosynthetic_class(atlas_df: pd.DataFrame, target_level: str, target_class: str) -> pd.DataFrame:
    """Filter Atlas df by biosynthetic class and return copy of filtered df"""

    filtered_df = atlas_df[atlas_df[target_level] == target_class].copy()

    return filtered_df


if __name__ == "__main__":
    atlas_df = convert_json(ATLAS_JSON_PATH)
    biosynthetic_class_df = extract_biosynthetic_class(atlas_df,
                                                       "pathway_results",
                                                       "Polyketides")
    print(biosynthetic_class_df)