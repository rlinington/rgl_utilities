#!/usr/bin/env python3

"""Convert NP Atlas json download to Pandas df"""

import pandas as pd


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

    return expanded_atlas_df


def convert_tsv(tsv_path):
    """Open tsv file as Pandas df

    Input: path to tsv
    Returns: pandas df"""

    with open(tsv_path) as f:
        atlas_df = pd.read_csv(f, sep="\t")

    return atlas_df
