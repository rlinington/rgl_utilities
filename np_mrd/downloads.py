#!usr/bin/env python3

"""Tools for examining the quality of download files form NP-MRD"""

import pandas as pd
from pathlib import Path
from rdkit import Chem
import csv


def import_csv_files(csv_dir: Path, report_path: Path) -> pd.DataFrame:
    """Import csv download files and return as single Dataframe"""

    csv_dfs = []

    for file in Path.glob(csv_dir, "*.csv"):
        input_data = pd.read_csv(file)
        csv_dfs.append(input_data)

    full_input_data = pd.concat(csv_dfs, ignore_index=True)

    with open(report_path, "a") as f:
        f.write(f"Number of csv compounds: {len(full_input_data)}\n")

    return full_input_data


def smiles_to_inchikey(smiles: str) -> str:
    """Convert smiles to inchikey or None if smiles not valid"""

    try:
        mol = Chem.MolFromSmiles(smiles)
    except:
        print("Failed to create mol")
        return None

    try:
        return Chem.MolToInchiKey(mol)
    except:
        print("Failed to create InChIKey")
        return None


def csv_add_inchikey(csv_df: pd.DataFrame) -> pd.DataFrame:
    """Add a column containing inchikey to df"""

    csv_df["inchikey"] = csv_df["SMILES"].apply(smiles_to_inchikey)

    return csv_df


def import_sdf_files(sdf_dir: Path, report_path: Path) -> list[Chem.Mol]:
    """Import sdf files and return list of rdkit mol objects"""

    mol_objects = []
    inchikey_mismatch = 0

    for file in Path.glob(sdf_dir, "*.sdf"):
        input_data = Chem.SDMolSupplier(str(file))
        for entry in input_data:
            if entry and entry.GetProp("INCHI_KEY") != Chem.MolToInchiKey(entry):
                inchikey_mismatch += 1
            mol_objects.append(entry)

    with open(report_path, "a") as f:
        f.write(f"InChIKey mismatch between SDF InChIKey and calculated InChIKey from mol block = {inchikey_mismatch}\n")
        f.write(f"Number of sdf compounds: {len(mol_objects)}\n")

    return mol_objects


def compare_sdf_smiles_inchikey(sdf_data: list[Chem.Mol], csv_dict: dict, report_path):
    """Compare InChIKey from sdf to InChIKey derived from SMILES. Count matches and mismatches"""
    match_count = 0
    mismatch_count = 0
    missing_smiles = 0
    failed_inchikey = 0
    bad_npmrd_id = 0

    for entry in sdf_data:
        try:
            inchikey = Chem.MolToInchiKey(entry)
        except:
            failed_inchikey += 1
            continue
        try:
            npmrd_id = entry.GetProp("NP-MRD_ID")
            if len(npmrd_id) != 9:
                bad_npmrd_id += 1
                continue
        except:
            bad_npmrd_id += 1
            continue
        try:
            smiles_inchikey = csv_dict[npmrd_id]
        except KeyError:
            missing_smiles += 1
            continue
        if inchikey == smiles_inchikey:
            match_count += 1
        else:
            mismatch_count += 1

    with open(report_path, "a") as f:
        f.write(f"SDF to SMILES InChIKey match count: {match_count}\n")
        f.write(f"SDF to SMILES InChIKey mismatch count: {mismatch_count}\n")
        f.write(f"SMILES NP-MRD ID missing: {missing_smiles}\n")
        f.write(f"SDF InChIKey creation fail: {failed_inchikey}\n")
        f.write(f"SDF bad NP-MRD ID: {bad_npmrd_id}\n")


def count_mismatched_inchikeys_in_atlas(sdf_dir: Path, atlas_path: Path):
    """Examine each entry in sdf. If saved InChIKey and 3D mol calculated InChIKey do not match then search both
    InChIKeys against Atlas. Count matches in each category"""

    atlas_df = pd.read_csv(atlas_path, sep='\t')

    mismatch_count = 0
    stored_match = 0
    calculated_match = 0
    both_match = 0
    atlas_npmrd_matches = [["npmrd_id",
                            "npmrd_stored_structure",
                            "npmrd_3d_mol_structure",
                            "npaid_match_list",
                            "match_type"]]

    for file in Path.glob(sdf_dir, "*.sdf"):
        input_data = Chem.SDMolSupplier(str(file))
        for compound in input_data:
            match_type = None
            atlas_match = None
            try:
                sdf_stored_inchikey = compound.GetProp("INCHI_KEY")
                sdf_3d_inchikey = Chem.MolToInchiKey(compound)
            except:
                print("ERROR with compound import")
                continue
            if sdf_stored_inchikey != sdf_3d_inchikey:
                mismatch_count += 1
                stored_atlas_match = atlas_df[atlas_df["compound_inchikey"] == sdf_stored_inchikey]
                calculated_atlas_match = atlas_df[atlas_df["compound_inchikey"] == sdf_3d_inchikey]
                stored_atlas_match_count = len(stored_atlas_match)
                calculated_atlas_match_count = len(calculated_atlas_match)
                if stored_atlas_match_count > 0 and calculated_atlas_match_count == 0:
                    stored_match += 1
                    match_type = "stored_match"
                    atlas_match = stored_atlas_match["npaid"].values.tolist()
                elif stored_atlas_match_count == 0 and calculated_atlas_match_count > 0:
                    calculated_match += 1
                    match_type = "calculated_match"
                    atlas_match = calculated_atlas_match["npaid"].values.tolist()
                elif stored_atlas_match_count > 0 and calculated_atlas_match_count > 0:
                    both_match += 1
                    match_type = "both"
                    atlas_match = (stored_atlas_match["npaid"].values.tolist() +
                                   calculated_atlas_match["npaid"].values.tolist())
                if match_type is not None:
                    try:
                        atlas_npmrd_matches.append([compound.GetProp("NP-MRD_ID"),
                                                    compound.GetProp("SMILES"),
                                                    Chem.MolToSmiles(compound),
                                                    atlas_match,
                                                    match_type])
                    except:
                        print("Failed to store match data")
    print(f"NP-MRD InChIKey mismatch count: {mismatch_count}")
    print(f"Stored InChIKey match to Atlas: {stored_match}")
    print(f"Calculated InChIKey from 3D mol match to Atlas: {calculated_match}")
    print(f"Both stored and calculated InChIKey match to Atlas: {both_match}")

    output_path = Path(sdf_dir.parent, "npmrd_inchikey_mismatches_atlas_matching_compounds.tsv")

    with open(output_path, "w", encoding='utf-8') as g:
        csv_g = csv.writer(g, delimiter='\t', lineterminator='\n')
        csv_g.writerows(atlas_npmrd_matches)


if __name__ == "__main__":

    data_path = Path("/Users/roger/Git/rgl_utilities/data/np_mrd/data_dump_20240920")
    report_path = Path(data_path, "quality_report_for_" + data_path.stem + ".txt")
    atlas_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.tsv")

    # Perform basic quality control checks on csv and sdf download files frm NP-MRD
    # with open(report_path, "w") as f:
    #     f.write(f"Structure Quality Report for {data_path.stem}\n")
    #
    # csv_data = import_csv_files(Path(data_path, "smiles"), report_path)
    # sdf_data = import_sdf_files(Path(data_path, "3d_mol"), report_path)
    #
    # csv_data = csv_add_inchikey(csv_data)
    # filtered_csv_data = csv_data.dropna()
    # with open(report_path, "a") as f:
    #     f.write(f"Number of csv compounds after InChIKey generation: {len(filtered_csv_data)}\n")
    # csv_npmrd_inchikey_dict = dict(zip(filtered_csv_data["NP_MRD_ID"], filtered_csv_data["inchikey"]))
    # compare_sdf_smiles_inchikey(sdf_data, csv_npmrd_inchikey_dict, report_path)

    # Compare compounds with mismatching InChIKeys to Atlas and count cases where structures match
    count_mismatched_inchikeys_in_atlas(Path(data_path, "3d_mol"),
                                        atlas_path)
