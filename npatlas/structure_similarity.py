#!/usr/bin/env python3

"""Tools to compare new molecules to the chemical space in the NP Atlas DB"""

import pandas as pd
from pathlib import Path

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs


def generate_fingerprint(smiles):
    """Calculate molecular fingerprint using Morgan fingerprint radius 2"""

    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprint(mol, 2)


def calculate_similarity(atlas_fp, input_fp):
    """Calculate similarity score between two rdkit fingerprints using Dice similarity scoring"""

    return round(DataStructs.DiceSimilarity(atlas_fp, input_fp), 2)


if __name__ == "__main__":

    atlas_tsv_path = Path("/Users/roger/Git_new/rgl_utilities/data/NP_Atlas/NPAtlas_download_2023_06.tsv")
    new_compound_name = "Megapolypeptin B"
    new_compound_smiles = "CC(C)C(=O)C(=O)NC(CC(=O)NC(C(=O)NC(C(=O)NC(CO)C(O)CC(N)=O)C(C)O)C(C)O)CC(O)CC/C=C/CC/C=C/CC/C=C/CCCC(C)OC(CC(=O)CCC(=O)O)C(=O)O"
    new_compound_fingerprint = generate_fingerprint(new_compound_smiles)

    input_df = pd.read_csv(atlas_tsv_path, delimiter='\t')
    input_df["molecule_fingerprint"] = input_df["compound_smiles"].map(generate_fingerprint)
    input_df["Dice_similarity_score"] = input_df["molecule_fingerprint"].apply(calculate_similarity,
                                                                               args=(new_compound_fingerprint, ))
    input_df = input_df.drop(["molecule_fingerprint"], axis=1)

    output_path = Path(atlas_tsv_path.parent, new_compound_name + "_atlas_structure_similarity.tsv")
    input_df.to_csv(output_path, index=False, sep='\t', lineterminator='\n')
