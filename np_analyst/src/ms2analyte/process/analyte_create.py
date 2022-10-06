#!/usr/bin/env python3

"""Tools to build analytes from peak output from peak_create.py"""

class MassPeak:
    """Class containing all relevant data for each individual mass peak from peak_create"""
    def __init__(self, dataframe, peak_id, max_intensity, peak_area, min_scan, max_scan, scan_array, average_mass, rt,
                 average_drift, peak_assigned):
        self.dataframe = dataframe
        self.peak_id = peak_id
        self.max_intensity = max_intensity
        self.peak_area = peak_area
        self.min_scan = min_scan
        self.max_scan = max_scan
        self.scan_array = scan_array
        self.average_mass = average_mass
        self.rt = rt
        self.average_drift = average_drift
        self.peak_assigned = peak_assigned
        self.dda_data = None                # [(mz_array), (intensity_array)]. NOTE: If the same mass is selected for
                                            # DDA multiple times, the scan with the highest DDA intensity is retained


class Analyte:
    """Class containing data on peak ids for each assigned analyte"""
    def __init__(self, analyte_id, peak_list):
        self.analyte_id = analyte_id
        self.ms1_peak_list = peak_list          # List of MS1 MassPeak objects
        self.ms2_peak_list = None               # List of MS2 MassPeak objects from DIA experiments
        self.max_peak_id = None                 # peak_id of MS1 peak with max intensity
        self.max_peak_intensity_mass = None     # mz value for MS1 peak with max intensity
        self.max_peak_intensity = None          # intensity of MS1 peak with max intensity
        self.max_peak_area = None
        self.sum_peak_area = None
        self.max_peak_scan = None               # Scan number for max intensity data point from max peak
        self.analyte_rt = None
        self.replicate_match = None
        self.replicate_analyte_id = None
        self.experiment_match = None
        self.experiment_analyte_id = None
        self.blank_match = None
