#!/usr/bin/env python3

"""Convert NP Atlas json download to Pandas df"""

import json
import pandas as pd


def convert_json(json_path):
    """Open json file as Pandas df

    Input: path to json
    Returns: pandas df"""

    with open(json_path, "rb") as f:
        atlas_df = pd.read_json(f)

    return atlas_df


def convert_tsv(tsv_path):
    """Open tsv file as Pandas df

    Input: path to tsv
    Returns: pandas df"""

    with open(tsv_path) as f:
        atlas_df = pd.read_csv(f, sep="\t")

    return atlas_df
