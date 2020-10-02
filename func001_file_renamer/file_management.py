"""Tools to manage copying and file renaming for func001 files from the Linington Lab for use with MS2Analyte"""

import glob
import os
import sys
import csv
from shutil import copy
from collections import defaultdict


def extract_original_sample_name(path_data):
    """Extract original sample name from directory name"""

    return os.path.basename(path_data).rsplit(".")[0]


def strip_datestamp(original_sample_name):
    """Remove leading datestamp from sample names"""

    return original_sample_name[9:]


def extract_sample_list(source_data_directory):
    """Create list of unique sample names from dataset"""

    unique_sample_list = []

    for directory in glob.glob(os.path.join(source_data_directory, "*.raw")):
        original_sample_name = extract_original_sample_name(directory)
        sample_name_minus_datestamp = strip_datestamp(original_sample_name)
        if sample_name_minus_datestamp not in unique_sample_list:
            unique_sample_list.append(sample_name_minus_datestamp)

    return unique_sample_list


def replicate_filename_counter(unique_sample_list):
    """Create dict of sample names to number replicates sequentially"""

    sample_replicate_counter = {}

    for sample_name in unique_sample_list:
        sample_replicate_counter[sample_name] = 1

    return sample_replicate_counter


def replicate_check(source_data_directory, unique_sample_list, replicate_count):
    """Checker to make sure all samples have the same number of replicates"""

    replicate_counter = defaultdict(int)
    replicate_count_error = False

    for directory in glob.glob(os.path.join(source_data_directory, "*.raw")):
        original_sample_name_minus_datestamp = strip_datestamp(extract_original_sample_name(directory))
        replicate_counter[original_sample_name_minus_datestamp] += 1

    for sample_name in unique_sample_list:
        if replicate_counter[sample_name] != replicate_count:
            print("ERROR: incorrect number of replicates for sample " + sample_name + ". Found " +
                  str(replicate_counter[sample_name]) + " expected " + str(replicate_count))
            replicate_count_error = True

    if not replicate_count_error:
        print("Replicate count check PASSED")
    else:
        print("ERROR: One or more samples have the wrong number of replicates. Please correct and rerun naming script")
        sys.exit()


def create_sw_rl_name_dict(sw_name_csv):
    """Create dict to convert legacy 'SW' coding back to our standard 'RL' naming for samples"""

    pass


def create_ms2analyte_func001_copy(source_data_directory, destination_directory, replicate_count, sw_naming_file):
    """Find every directory containing MS data in the source directory, and use the directory name to rename the func001
     file and copy it to a new location

     """

    if not os.path.isdir(destination_directory):
        os.mkdir(destination_directory)

    unique_sample_list = extract_sample_list(source_data_directory)
    sample_replicate_counter = replicate_filename_counter(unique_sample_list)
    replicate_check(source_data_directory, unique_sample_list, replicate_count)
    sw_name_dict = create_sw_rl_name_dict(sw_naming_file)

    renaming_log = []
    headers = [["original_sample_name", "ms2analyte_filename"]]

    for directory in glob.glob(os.path.join(source_data_directory, "*.raw")):
        original_sample_name = extract_original_sample_name(directory)
        sample_name_minus_datestamp = strip_datestamp(original_sample_name)

        if sample_name_minus_datestamp[:2] == "RL":
            new_filename = sample_name_minus_datestamp + "_R" + \
                           str(sample_replicate_counter[sample_name_minus_datestamp])
            sample_replicate_counter[sample_name_minus_datestamp] += 1
        elif sample_name_minus_datestamp[:2] == "SW":
            new_filename = sw_name_dict[sample_name_minus_datestamp] + "_R" + \
                           str(sample_replicate_counter[sample_name_minus_datestamp])
            sample_replicate_counter[sample_name_minus_datestamp] += 1
        else:
            print("ERROR: Unrecognized filename (not SW or RL) for " + sample_name_minus_datestamp)
            sys.exit()

        copy(os.path.join(directory, "func001.csv"), os.path.join(destination_directory, new_filename))
        renaming_log.append([original_sample_name, new_filename])

    with open(os.path.join(destination_directory, "file_renaming_log.csv", "w")) as f:
        csv_f = csv.writer(f)
        csv_f.writerows(headers + renaming_log)


if __name__ == "__main__":

    source_data_directory = os.path.join("text")
    destination_directory = os.path.join("text2")
    replicate_count = 3
    sw_naming_file = os.path.join("text3")

    create_ms2analyte_func001_copy(source_data_directory, destination_directory, replicate_count, sw_naming_file)
