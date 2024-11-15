#!usr/bin/env python3

import src.nmr_ase.tools.compounds
from rdkit import Chem
import pickle
from pathlib import Path
import csv
import numpy as np
from operator import itemgetter

from np_mrd import nmr_search

COMPOUND_PATH = Path("/Users/roger/Git/nmr-ase/data/output/full_compounds/full_compound_list.pickle")
DATA_PATH = Path("/Users/roger/Git/rgl_utilities/data/np_mrd")
SHIFT_VARIANCE = 0.3

if __name__ == "__main__":
    with open(COMPOUND_PATH, "rb") as f:
        compound_list = pickle.load(f)

    bitstring_list = []
    for compound in compound_list:
        results = [f"{compound.source_type}: {compound.source_id}"]
        if compound.c_shifts is not None:
            results.append(nmr_search.shifts_to_fingerprint(compound.c_shifts.values(), SHIFT_VARIANCE))
        else:
            results.append([0] * 2300)
        bitstring_list.append(results)

    # with open(Path(DATA_PATH, "spectrum_bitstring_list.csv"), "w") as g:
    #     csv_g = csv.writer(g)
    #     csv_g.writerows(bitstring_list)

    # Source: jeol, id: 8931, smiles: O=CN(O)CCCC(NC(=O)C(CO)NC(=O)C1COC(c2ccccc2O)=N1)C(=O)NCCC(=O)NC1CCCN(O)C1=O
    # shifts: dict_values([51.2, 164.9, 49.42, 26.63, 20.25, 169.98, 35.22, 35.38, 171.12, 52.6, 28.55, 22.52, 48.77,
    # 161.74, 169.35, 55.32, 169.09, 67.24, 69.59, 166.07, 109.33, 159.07, 116.61, 134.0, 119.07, 128.07, 61.72])
    # test_spectrum = [51.2, 164.9, 49.42, 26.63, 20.25, 169.98, 35.22, 35.38, 171.12, 52.6, 28.55, 22.52, 48.77, 161.74,
    #                  169.35, 55.32, 169.09, 67.24, 69.59, 166.07, 109.33, 159.07, 116.61, 134.0, 119.07, 128.07, 61.72]

    # Source: https://doi.org/10.1021/acs.jnatprod.0c00444
    # smiles: O=C(C(C(C1(C)C2C(CC(C)CC2)C=C(C)C1CCC)=O)=C3O)NC3=O
    # spectrum: [200.2, 49.0, 45.2, 126.7, 130.1, 38.7, 42.4, 33.6, 35.8, 28.5, 40.1, 14.1, 131.0, 127.3, 18.0, 22.6,
    # 179.6, 100.6, 191.2, 65.8, 67.8, 19.7]
    # predicted spectrum from NP-MRD: [17.90008545, 127.2003174, 130.7998657, 44.99969482, 126.600647, 130.0003052,
    # 38.49945068, 42.1005249, 33.39996338, 22.50061035, 35.60028076, 28.30047607, 39.99938965, 49.00054932,
    # 13.89923096, 99.90081787, 176.1001587, 67.79937744, 67.39959717, 19.74029541, 192.7993774]
    test_spectrum = [200.2, 49.0, 45.2, 126.7, 130.1, 38.7, 42.4, 33.6, 35.8, 28.5, 40.1, 14.1, 131.0, 127.3, 18.0,
                     22.6, 179.6, 100.6, 191.2, 65.8, 67.8, 19.7]

    test_bitstring = nmr_search.shifts_to_fingerprint(test_spectrum, SHIFT_VARIANCE)

    match_list = []
    for reference_compound in bitstring_list:
        match_matrix = np.corrcoef(np.array(test_bitstring), np.array(reference_compound[1]))
        match_score = match_matrix[0, 1]
        if match_score > 0.5:
            match_list.append([reference_compound[0], round(float(match_score), 2)])

    print(sorted(match_list, key=itemgetter(1), reverse=True))
