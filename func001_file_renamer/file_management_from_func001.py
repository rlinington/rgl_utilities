"""Tools to manage copying and file renaming for func001 files from the Linington Lab for use with MS2Analyte

Designed for use with legacy func001 input files in the general format RLxx_9999x-1_func001.csv

"""

import glob
import os
import sys
import csv
from shutil import copy
from collections import defaultdict
from os.path import expanduser


def replace_sample_name(filepath):

    filename = str(os.path.basename(filepath).rsplit(".")[0])
    new_filename = filename.replace('_func001', '')
    new_filename = new_filename.replace('-', '_R')

    return new_filename


def create_sw_dict(sw_csv_file):

    sw_dict = {}

    with open(sw_csv_file) as f:
        csv_f = csv.reader(f)
        next(f)

        for row in csv_f:
            sw_dict[row[0].split("-")[0]] = row[6].replace("-", "_")

    return sw_dict


if __name__ == "__main__":

    source_data_directory = os.path.expanduser(os.path.join("~", "Documents", "Collaborators", "Eustaquio", "BKD_MS2Analyte", "PKS_NRPS_BGCs", "Blanks"))
    destination_directory = os.path.expanduser(os.path.join("~", "Documents", "Collaborators", "Eustaquio", "BKD_MS2Analyte", "PKS_NRPS_BGCs", "Renamed_Blanks"))
    # sw_naming_file = os.path.join("rgl_utilities", "data", "func001_file_renamer", "sw_converter_list.csv")

    # sw_renaming_dict = create_sw_dict(sw_naming_file)

    for funcfile in glob.glob(os.path.join(source_data_directory, "*_func001.csv")):
        print(funcfile)

        revised_sample_name = replace_sample_name(funcfile)

        # if revised_sample_name[:2] == "SW":
        #     sample_id = revised_sample_name.split("_")[0]
        #     replicate_id = revised_sample_name.split("_")[1]
        #
        #     new_sample_name = sw_renaming_dict[sample_id]
        #     revised_sample_name = new_sample_name + "_" + replicate_id

        copy(funcfile, os.path.join(destination_directory, revised_sample_name + ".csv"))
