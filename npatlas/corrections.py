#!usr/bin/env python3

"""Tools to manage user depositions and corrections"""

import pandas as pd
import json
from pathlib import Path


def json_to_csv(json_path: Path):
    """Convert any json file into a csv copy and save"""
    with open(json_path) as f:
        input_data = pd.read_json(f)

    output_directory = json_path.parent
    filename = json_path.stem
    output_path = Path(output_directory, f"{filename}.csv")

    input_data.to_csv(output_path, lineterminator='\n', index=False, encoding='utf8')


if __name__ == "__main__":
    json_source_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/np_atlas_all_corrections_25-04-2024.json")
    json_to_csv(json_source_path)