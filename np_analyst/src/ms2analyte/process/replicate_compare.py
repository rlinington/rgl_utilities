#!/usr/bin/env python3

"""Tools to replicate compare analytes from analyte_create.py"""


class ReplicateMaxMassData:
    """Class containing all mass data from maximum intensity scan for one analyte from one replicate"""
    def __init__(self, replicate_number, sample_name, peak_list, max_scan, max_intensity, normalized_mass_data):
        self.replicate_number = replicate_number
        self.sample_name = sample_name
        self.peak_list = peak_list
        self.max_scan = max_scan
        self.max_intensity = max_intensity
        self.normalized_mass_data = normalized_mass_data


class ReplicateMassPeak:
    """Class containing data about each mass peak in a replicate analyte"""
    def __init__(self, mass_peak_data,
                 mass_peak_members,
                 average_mass,
                 average_intensity,
                 relative_intensity,
                 replicate_peak_id):
        self.mass_peak_data = mass_peak_data    # [[mz, intensity, replicate number]]
        self.mass_peak_members = mass_peak_members  #[[peak_id, replicate]]
        self.average_mass = average_mass
        self.average_intensity = average_intensity
        self.relative_intensity = relative_intensity
        # peak id for consensus mass peaks within a replicate analyte. Starts at 1 for each replicate analyte
        self.replicate_peak_id = replicate_peak_id


class ReplicateAnalyte:
    """Class containing data on analyte ids for each assigned replicate analyte (i.e. analytes that pass replicate
    comparison)

    """
    def __init__(self, replicate_analyte_id, replicate_analyte_members, analyte_list, max_peak_intensity_mass,
                 max_intensity, max_peak_area_mean, max_peak_area_rsd, sum_peak_area_mean, sum_peak_area_rsd,
                 replicate_analyte_rt, replicate_analyte_ms1_spectrum, replicate_analyte_ms2_spectrum, sample_name,
                 is_basketed):
        self.replicate_analyte_id = replicate_analyte_id
        self.replicate_analyte_members = replicate_analyte_members
        self.analyte_list = analyte_list                                # [[replicate, analyte id]]
        self.max_peak_intensity_mass = max_peak_intensity_mass
        self.max_intensity = max_intensity
        self.max_peak_area_mean = max_peak_area_mean
        self.max_peak_area_rsd = max_peak_area_rsd
        self.sum_peak_areas_mean = sum_peak_area_mean
        self.sum_peak_areas_rsd = sum_peak_area_rsd
        self.replicate_analyte_rt = replicate_analyte_rt
        self.replicate_analyte_ms1_spectrum = replicate_analyte_ms1_spectrum    # [ReplicateMassPeak list]
        self.replicate_analyte_ms2_spectrum = replicate_analyte_ms2_spectrum    # [ReplicateMassPeak list]
        self.sample_name = sample_name
        self.is_basketed = is_basketed   # This is used in basketing step to remove matched analytes from available pool
