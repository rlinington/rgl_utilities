#!usr/bin/env python3

from pathlib import Path


COMPOUND_PATH = Path("/Users/roger/Git/nmr-ase/data/output/full_compounds/full_compound_list.pickle")
DATA_PATH = Path("/Users/roger/Git/rgl_utilities/data/np_mrd")


def shifts_to_fingerprint(chemical_shifts: list, shift_variance: float):
    """Convert list of 13C chemical shifts to standardized bitstring that includes defined variance"""

    shift_bitstring = [0] * 2300
    for shift in chemical_shifts:
        if 0 <= shift <= 229.9:
            for i in range(int(-shift_variance*10), int(shift_variance*10), 1):
                shift_bitstring[int(shift*10) + i] = 1

    return shift_bitstring

