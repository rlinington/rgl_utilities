#!/usr/bin/env python3

"""Tool to extract func001 files from MSeXpress output into a single directory, and rename them according to each
directory name"""

import os
import shutil
import sys


def define_paths():
    """Obtain paths for source and destination directories for func001 data"""

    source = os.path.normpath(input("Enter the path for the directory containing the MSeXpress results: "))
    destination = os.path.normpath(input("Enter the destination directory: "))
    if not os.path.exists(destination):
        create_output_dir = input("Destination directory does not exist. Create? [y/n]: ")
        if create_output_dir == "y":
            os.mkdir(destination)
        else:
            sys.exit()
    if source == destination:
        response = input("WARNING: Source and destination directories the same. "
                         "Are you sure you wish to continue? [y/n]: ")
        if response == "y":
            pass
        else:
            sys.exit()

    return source, destination


def extract_func001(source, destination):
    """Copy func001 files from individual directories into the destination directory and rename each file using the
    parent directory name"""

    for entry in os.listdir(source):
        if os.path.isdir(os.path.join(source, entry)):
            if "func001.csv" in os.listdir(os.path.join(source, entry)):
                shutil.copy(os.path.join(source, entry, "func001.csv"),
                            os.path.join(destination, entry + ".csv"))


if __name__ == "__main__":
    source, destination = define_paths()
    # source = os.path.normpath("/Users/roger/Git/rgl_utilities/data/func001_file_renamer/func001_extractor_test_data/source_dir")
    # destination = os.path.normpath("/Users/roger/Git/rgl_utilities/data/func001_file_renamer/func001_extractor_test_data/destination_dir")
    extract_func001(source, destination)
