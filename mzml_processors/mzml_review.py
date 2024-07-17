#!usr/bin/env python3

"""Tools for preliminary evaluation and assessment of data from mzML files"""

import pymzml
import sys
import pandas as pd
from pathlib import Path


def extract_scan_by_scan_number(file_path: Path, output_path: Path, scan_number: int):

    """Import mzML files derived from applying MSConvert to .raw files.

    Args:
        file_path (Path): path to mzML file
        output_path (Path): path to output_directory
        scan_number (int): The scan number for which MS data output is required

    Waters data includes the lockspray internal calibrant scans as 'MS1' data. These are differentiated from true
    MS1 data by the 'function' attribute in the spectrum element. Data MS1 scans are function 1. Lockspray scans are
    assigned the highest possible function number (floating, depends on how many DDA scans were permitted during
    acquisition setup). Commonly lockspray function=5. This is always 3 for MSe (DIA) data.
    In order to filter out the lockspray data it is therefore necessary to filter the ms_level 1 scans to only
    retain those where the function value is also 1.
    NOTE: For MS2 data this is not an issue, because all MS2 data have ms level = 2, and are therefore all
    legitimate for inclusion.
    """
    run = pymzml.run.Reader(str(file_path))
    spec = run[scan_number]
    scan_data = pd.DataFrame([spec.mz, spec.i])
    pivot_scan_data = scan_data.transpose()
    pivot_scan_data.columns = ["mz", "intensity"]
    pivot_scan_data = pivot_scan_data[pivot_scan_data["intensity"] > 0]

    output_filename = f"{file_path.stem}_scan_{scan_number}.csv"
    output_file_path = Path(output_path, output_filename)
    pivot_scan_data.to_csv(output_file_path, index=False, lineterminator='\n', encoding='utf-8')


def rt_to_scan_number(file_path: Path, rt: float, ms_type=1) -> int:
    """For a given retention time, return the scan number that is closest to the desired rt"""
    run = pymzml.run.Reader(str(file_path))

    rt_difference = 10000
    rt_match_scan = None
    for spec in run:
        # Handle MS1 and MS2 scans separately
        if spec.ms_level == ms_type:
            if ms_type == 1:
                # Skip lockspray or other functions if there are any
                # If not, this is probably not Waters data and should be fine...
                fn = spec.id_dict.get("function")
                if fn is not None:
                    if fn != 1:
                        continue
            scan_rt = spec.scan_time_in_minutes()
            if abs(rt - scan_rt) < rt_difference:
                rt_match_scan = spec.ID
            if scan_rt > rt:
                break

    return rt_match_scan


def extract_rt_range(file_path: Path, output_path: Path, min_rt: float, max_rt: float):

    """Tool to extract MS1 data for rt range and export as csv.

    Args:
        file_path (Path): path to mzML file
        output_path (Path): path to output_directory
        min_rt (float): minimum rt value to be included
        max_rt (float): maximum rt value to be included
    """

    run = pymzml.run.Reader(str(file_path))
    ms1_input_data = []
    for spec in run:
        # Process MS1 only
        if spec.ms_level == 1:
            # Skip lockspray or other functions if there are any
            # If not, this is probably not Waters data and should be fine...
            fn = spec.id_dict.get("function")
            if fn is not None:
                if fn != 1:
                    continue
            if min_rt <= spec.scan_time_in_minutes() <= max_rt:
                ms1_input_data.append(extract_scan_data(spec))

    output_df = pd.concat(ms1_input_data, ignore_index=True)
    output_df[output_df["intensity"] > 0].to_csv(Path(output_path, f"{file_path.stem}_{min_rt}_{max_rt}.csv"),
                     lineterminator='\n',
                     index=False)


def extract_scan_data(spec):
    """Extract data from individual scans from mzML files including drift time"""

    scan_data = pd.DataFrame(spec.mz, columns=["mz"])
    scan_data["intensity"] = spec.i
    scan_data["scan"] = spec.ID
    scan_data["rt"] = round(spec.scan_time_in_minutes(), 4)

    scan_data = scan_data.astype({'scan': 'int64',
                                  'rt': 'float64',
                                  'mz': 'float64',
                                  'intensity': 'int64'})

    return scan_data


if __name__ == "__main__":
    mzml_file = Path("/Users/roger/Git/rgl_utilities/data/mzml/mrt/Megapolipeptin_A_mz_958.mzML")
    output_dir = Path("/Users/roger/Git/rgl_utilities/data/mzml/mrt/output")

    target_rt = 3.489
    target_window = 0.5

    target_scan = rt_to_scan_number(mzml_file, target_rt)
    print(target_scan)
    extract_scan_by_scan_number(mzml_file, output_dir, target_scan)

    # extract_rt_range(mzml_file, output_dir, target_rt - target_window, target_rt + target_window)
