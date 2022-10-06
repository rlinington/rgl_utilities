#!/usr/bin/env python3

"""Tools to explore UV data from legacy database files"""

import os
import pandas as pd


def extract_uv_data(compound_list):
    """Split text for UV data into int list"""


if __name__ == "__main__":

    import_file = os.path.join("data", "legacy_db", "AntiMarin_info.csv")
    import_df = pd.read_csv(import_file)

    names_uv = import_df[0]

