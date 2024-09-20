#!usr/bin/env python3

"""Tools to reformat output CSV files from BioTek Neo plate reader in HTCB"""


import csv
import string

import pandas as pd
from pathlib import Path
import numpy as np


def extract_raw_values(input_df: pd.DataFrame, timepoint):
    """Extract raw plate values from dataframe from initial BioTek plate reader output (i.e. pd import of .xlsx"""
    plate_reader_values = input_df.iloc[26:, 2:26].values.tolist()
    list_data = []
    for row_index, row_data in enumerate(plate_reader_values):
        for column_index, value in enumerate(row_data):
            list_data.append([f"{string.ascii_uppercase[row_index]}{column_index + 1:02d}", value])

    plate_data = pd.DataFrame(list_data, columns=["well", f"{timepoint}_value"])

    return plate_data


def plate_reader_to_df(file_path: Path) -> pd.DataFrame:
    """Convert standard BioTek Neo output to flat file format

    Args:
        file_path (Path): The path to the original .xlsx BioTek output file

    Returns:
        pandas DataFrame"""

    input_df = pd.read_excel(file_path)
    plate_info = file_path.stem.split("_")
    plate_number = plate_info[0]
    strain = plate_info[1]

    plate_data = extract_raw_values(input_df, "t0")
    plate_data["plate_number"] = plate_number
    plate_data["strain"] = strain

    plate_data = plate_data[["strain", "plate_number", "well", "t0_value"]]

    t20_path = Path(file_path.parents[1], "t20", file_path.name)
    t20_df = pd.read_excel(t20_path)
    t20_data = extract_raw_values(t20_df, "t20")

    plate_data = pd.merge(plate_data, t20_data, on="well")
    plate_data["difference"] = plate_data["t20_value"] - plate_data["t0_value"]

    return plate_data


def normalize_absorbance_data(input_df: pd.DataFrame, filepath: Path) -> pd.DataFrame:
    """Normalize plate data based on positive and negative control values"""

    raw_df = pd.read_excel(filepath)
    t0_positive_control_values = [x for xs in raw_df.iloc[26:34, 24:26].values.tolist() for x in xs]
    t0_negative_control_values = [x for xs in raw_df.iloc[34:, 24:26].values.tolist() for x in xs]

    t20_path = Path(filepath.parents[1], "t20", filepath.name)
    t20_df = pd.read_excel(t20_path)
    t20_positive_control_values = [x for xs in t20_df.iloc[26:34, 24:26].values.tolist() for x in xs]
    t20_negative_control_values = [x for xs in t20_df.iloc[34:, 24:26].values.tolist() for x in xs]

    normalized_positive_control = np.array(t20_positive_control_values) - np.array(t0_positive_control_values)
    normalized_negative_control = np.array(t20_negative_control_values) - np.array(t0_negative_control_values)
    mean_positive_value = normalized_positive_control.mean()
    mean_negative_value = normalized_negative_control.mean()

    input_df["percent_g"] = ((input_df["difference"] - mean_positive_value)/
                             (mean_negative_value - mean_positive_value)) * 100

    return input_df


def add_percent_inhibition(input_df: pd.DataFrame) -> pd.DataFrame:
    """Convert percent growth into percent inhibition value. Add column to normalize negative inhibition values to
    zero"""
    input_df["percent_inhibition"] = 100 - input_df["percent_g"]
    input_df["normalized_percent_inhibition"] = input_df["percent_inhibition"].clip(lower=0)
    input_df["normalized_percent_inhibition"] = input_df["percent_inhibition"].clip(upper=100)

    return input_df


def add_prefraction_numbers(source_dir):
    """Add the RL extract code to each row in the dataframe and overwrite csv file"""

    plate_name_path = Path(source_dir, "Plate_names.xlsx")
    plate_names = pd.read_excel(plate_name_path)
    plate_name_dict = plate_names.set_index('plate_label').to_dict()

    platemap_dict = {}

    platemap_dir = Path(source_dir, "plate_maps")
    for platemap_file in platemap_dir.glob("*.csv"):
        plate_map_df = pd.read_csv(platemap_file)
        plate_name = platemap_file.stem.split("-")[0]
        platemap_dict[plate_name_dict["plate_number"][plate_name]] = plate_map_df

    for output_file in Path(source_dir, "output", "plate_data").glob("*.csv"):
        results_df = pd.read_csv(output_file)
        plate_number = int(results_df["plate_number"].unique().tolist()[0][2:])
        output_df = pd.merge(results_df, platemap_dict[plate_number][["Well", "Prefraction"]],
                             left_on="well",
                             right_on="Well")
        output_df.rename(columns={'Prefraction':'prefraction'}, inplace=True)
        output_df = output_df[["strain",
                               "plate_number",
                               "well",
                               "prefraction",
                               "t0_value",
                               "t20_value",
                               "difference",
                               "percent_g",
                               "percent_inhibition",
                               "normalized_percent_inhibition"]]
        output_df = output_df.round(3)
        output_df.to_csv(output_file)


def create_npanalyst_input(source_dir: Path):
    """Create flat file with prefractions as rows and assay results as columns"""
    strain_list = ["EMPTY", "GES5", "IMP1", "KPC2", "NDM1", "OXA48", "VIM1"]
    df_list = []
    for output_file in Path(source_dir, "output", "plate_data").glob("*.csv"):
        results_df = pd.read_csv(output_file)
        df_list.append(results_df)

    full_data_df = pd.concat(df_list, ignore_index=True)
    filtered_df = full_data_df[full_data_df["prefraction"].notna()]
    npanalyst_df = pd.pivot_table(filtered_df,
                                  values="normalized_percent_inhibition",
                                  index=["prefraction"],
                                  columns=["strain"],
                                  aggfunc="sum")
    npanalyst_df.index.name = None
    npanalyst_df.rename(columns={"strain": "prefraction"}, inplace=True)
    output_path = Path(source_dir, "output", "npanalyst_bioassay_data.csv")
    npanalyst_df.to_csv(output_path)


if __name__ == "__main__":

    source_directory = Path("/Users/roger/Documents/Students/Emily_McMann/Carbapenemase_inhibitors/HPLC_subfractions_screen/Screening data")
    output_dir = Path(source_directory, "output")

    for file in Path(source_directory, "t0").glob("*.xlsx"):
        print(f"Starting file {file}")
        plate_data = plate_reader_to_df(file)
        plate_data = normalize_absorbance_data(plate_data, file)
        plate_data = add_percent_inhibition(plate_data)
        output_path = Path(output_dir, "plate_data", f"{file.stem}_processed.csv")
        plate_data.to_csv(output_path, encoding='UTF-8', lineterminator='\n', index=False)

    add_prefraction_numbers(source_directory)
    create_npanalyst_input(source_directory)