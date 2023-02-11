#!/usr/bin/env python3

"""Tools for refactoring flat files for importation to NP Analyst v1"""

import pandas as pd
from pathlib import Path
import pickle


def ms2analyte_to_npanalyst(flat_grid_df, analyte_dict):
    """Tool to rename MS2Analyte output grid with correct headers from NP Analyst MZMine import"""

    # headers = flat_grid_df.columns.to_list()
    # analyte_list = headers[1:]
    # rt_array = []
    # mz_array = []
    # for analyte in analyte_list:
    #     rt_array.append(analyte_dict[int(analyte)]["rt"])
    #     mz_array.append(analyte_dict[int(analyte)]["mz"])

    flat_grid_df["mz_mine_sample_name"] = flat_grid_df["sample_name"] + " Peak area"
    flat_grid_df = flat_grid_df.drop(["sample_name"], axis=1)
    col = flat_grid_df.pop("mz_mine_sample_name")
    flat_grid_df.insert(0, col.name, col)
    print(flat_grid_df)

    transposed_grid = flat_grid_df.transpose()
    print(transposed_grid)

    headers = transposed_grid.index.to_list()
    experiment_analyte_ids = headers[1:]

    rt_array = [None]
    mz_array = [None]

    for entry in experiment_analyte_ids:
        id = int(entry)
        rt_array.append(analyte_dict[id]["rt"])
        mz_array.append(analyte_dict[id]["mz"])

    transposed_grid["row retention time"] = rt_array
    transposed_grid["row m/z"] = mz_array

    col = transposed_grid.pop("row retention time")
    transposed_grid.insert(0, col.name, col)
    col = transposed_grid.pop("row m/z")
    transposed_grid.insert(0, col.name, col)

    print(transposed_grid)

    output_path = Path("/Users/roger/Git/rgl_utilities/data/npanalyst/MS2Analyte_flat_files/Garcinia NP analyte_npanalyst_mzmine_input.csv")
    transposed_grid.to_csv(output_path, index=False, line_terminator='\n')


def analyte_to_dict(analyte_list):
    """Convert Experiment analyte list to dict of rt and ms values"""
    analyte_dict = {}

    for analyte in analyte_list:
        analyte_inserted = False
        for mass_peak in analyte.experiment_analyte_spectrum:
            if mass_peak.relative_intensity == 100:
                analyte_dict[analyte.experiment_analyte_id] = {"rt": round(analyte.rt, 4), "mz": round(mass_peak.average_mass, 4)}
                analyte_inserted = True
                break
        if not analyte_inserted:
            analyte_dict[analyte.experiment_analyte_id] = {"rt": 0, "mz": 0}

    return analyte_dict


def create_np_analyst_input():
    """wrapper for ms2analyte_to_npanalyst function"""

    grid_path = Path("/Users/roger/Git/rgl_utilities/data/npanalyst/MS2Analyte_flat_files/Garcinia NP analyte_experiment_analyte_peak_area_grid.csv")
    analyte_path = Path("/Users/roger/Git/rgl_utilities/data/npanalyst/MS2Analyte_flat_files/Garcinia NP analyte_Samples_experiment_analytes.pickle")

    with open(analyte_path, "rb") as f:
        analyte_list = pickle.load(f)

    analyte_dict = analyte_to_dict(analyte_list)

    with open(grid_path) as g:
        grid_df = pd.read_csv(g)

    ms2analyte_to_npanalyst(grid_df, analyte_dict)


if __name__ == "__main__":

    create_np_analyst_input()

