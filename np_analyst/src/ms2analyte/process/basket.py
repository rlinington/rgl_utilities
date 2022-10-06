#!/usr/bin/env python3

"""Tools to basket replicate_analytes from replicate_compare.py"""


class ExperimentAnalyte:
    """Class for experiment analytes (i.e. an analyte found across a set of samples)"""
    def __init__(self, experiment_analyte_id, experiment_analyte_members, rt, experiment_analyte_spectrum):
        self.experiment_analyte_id = experiment_analyte_id
        self.experiment_analyte_members = experiment_analyte_members
        self.rt = rt
        self.experiment_analyte_spectrum = experiment_analyte_spectrum  # [ExperimentAnalyteMassPeak list]
        self.experiment_analyte_is_blank = False


class ExperimentAnalyteMassPeak:
    """Class for mass peaks for each experiment analyte (i.e. average relative intensities and masses)"""
    def __init__(self, average_mass, relative_intensity, contributing_peak_list):
        self.average_mass = average_mass
        self.relative_intensity = relative_intensity
        self.contributing_peak_list = contributing_peak_list    # (sample_name, replicatemasspeak)


class SampleExperimentId:
    """Class that defines the experiment analyte ids for each sample by sample analyte id. Used in
    input_data_sample_annotate

    """
    def __init__(self, sample_name, replicate_analyte_id_to_experiment_id):
        self.sample_name = sample_name
        self.replicate_analyte_id_to_experiment_id = replicate_analyte_id_to_experiment_id


class ExperimentAnalyteSpectrum:
    """Class for relative intensity spectra for experiment analyte id (averaged from all contributing sample analytes
    in that experiment analyte

    """
    def __init__(self, experiment_analyte_id, relative_experiment_mass_spectrum, experiment_analyte_is_blank):
        self.experiment_analyte_id = experiment_analyte_id
        self.relative_experiment_mass_spectrum = relative_experiment_mass_spectrum  # [ExperimentAnalyteMassPeak list]
        self.experiment_analyte_is_blank = experiment_analyte_is_blank
