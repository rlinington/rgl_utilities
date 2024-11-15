#!usr/bin/env python3

"""Tools to explore how SMARTS scaffolds map to structures in the NP Atlas"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from pathlib import Path
from collections import defaultdict
from operator import itemgetter
import csv

ATLAS_PATH = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.tsv")


def match_smarts_to_structure(structure, smarts):
    """Match a SMILES string to an rdkit SMARTS object. Return match state (True/ False)"""
    try:
        test_mol = Chem.MolFromSmiles(structure)
    except:
        print(f"Failed to generate mol for {structure}")
        return
    if len(test_mol.GetSubstructMatches(smarts)) > 0:
        return True
    else:
        return False


def generate_scaffold_match_image(atlas_df: pd.DataFrame, smarts: str, smarts_label: str, data_path: Path):
    """Create an image containing all the structures that have SMARTS match to provided smarts string"""

    smarts_mol = Chem.MolFromSmarts(smarts)
    atlas_df["smarts_match"] = atlas_df["compound_smiles"].apply(match_smarts_to_structure, args=(smarts_mol,))
    smarts_match_df = atlas_df[atlas_df["smarts_match"] == True]
    output_path = Path(data_path, f"Atlas_SMARTS_matches_{smarts_label}.tsv")
    smarts_match_df.to_csv(output_path, lineterminator="\n", sep="\t", encoding="utf-8", index=False)

    if len(smarts_match_df) > 0:
        smiles_list = smarts_match_df[["compound_smiles", "npaid"]].values.tolist()
        mol_list = [Chem.MolFromSmiles(x[0]) for x in smiles_list]
        for index, mol in enumerate(mol_list):
            mol.SetProp("_Name", smiles_list[index][1])

        AllChem.Compute2DCoords(smarts_mol)
        for mol in mol_list:
            _ = AllChem.GenerateDepictionMatching2DStructure(mol, smarts_mol)

        img = Draw.MolsToGridImage(mol_list, molsPerRow=4, subImgSize=(400, 400),
                                   legends=[x.GetProp("_Name") for x in mol_list])

        img.save(str(Path(data_path, f"Atlas_SMARTS_matches_{smarts_label}_images.png")))


def annotate_pks_substitution(compound_df: pd.DataFrame, smarts: str, target_smarts: str, output_dir: Path):
    """Map SMARTS to each structure and create new flat file containing the index position of each atom in the SMARTS,
    and the corresponding substitution at that position for each molecule"""

    position_annotations = defaultdict(defaultdict)

    for entry in compound_df.values.tolist():
        name = entry[2]
        smiles = entry[10]
        compound_mol = Chem.MolFromSmiles(smiles)
        smarts_mol = Chem.MolFromSmarts(smarts)
        if compound_mol.HasSubstructMatch(smarts_mol):
            labeled_mol = create_smarts_match_images(name,
                                                     compound_mol,
                                                     smarts_mol,
                                                     Path(output_dir, "smarts_structures", target_smarts))
            for atom in labeled_mol.GetAtoms():
                try:
                    label = atom.GetProp("atomLabel")
                    position_type = []
                    for neighbor in atom.GetNeighbors():
                        try:
                            neighbor.GetProp("atomLabel")
                        except:
                            neighbor_atom = neighbor.GetAtomicNum()
                            if neighbor_atom > 1:
                                position_type.append(str(neighbor_atom))
                    position_annotations[name][label] = ("").join(sorted(position_type))
                except:
                    pass
        else:
            print(f"No SMARTS match for {name}")

    output_data = [["compound_name"]]
    sorted_results_list = []

    for key1, value1 in position_annotations.items():
        results_list = []
        for key2, value2 in value1.items():
            results_list.append([int(key2), value2])
        annotations = [key1]
        sorted_results_list = sorted(results_list, key=itemgetter(0))
        for entry in sorted_results_list:
            annotations.append(entry[1])
        output_data.append(annotations)
    for entry in sorted_results_list:
        output_data[0].append(entry[0])

    with open(str(output_dir) + "/" + target_smarts + ".csv", "w") as f:
        csv_f = csv.writer(f)
        csv_f.writerows(output_data)


def create_smarts_match_images(compound_name: str,
                               compound_mol: Chem.Mol,
                               smarts_mol: Chem.Mol,
                               output_dir: Path):
    """Compare SMARTS against the compound_mol. Annotate each matching atom and bond in mol and generate image
    Args:
        compound_name (str): name of compound being analyzed
        compound_mol (Chem.Mol): compound being analyzed
        smarts_mol (Chem.Mol): smarts being analyzed
        output_dir (Path): The full directory path to which the image is to be saved.
    Returns:
        Image if SMARTS matched structure with SMRATS atom indices"""

    # Index information for atoms and bonds for SMARTS matches
    smarts_bonds = []
    smarts_atoms = []
    # Index information for SMARTS matches that also have data in shift dictionary
    shifts_bonds = []
    shifts_atoms = []

    if compound_mol.HasSubstructMatch(smarts_mol):
            matches = list(compound_mol.GetSubstructMatches(smarts_mol))
            for hit_ats in matches:
                smarts_atoms = smarts_atoms + list(hit_ats)
                for bond in smarts_mol.GetBonds():
                    aid1 = hit_ats[bond.GetBeginAtomIdx()]
                    aid2 = hit_ats[bond.GetEndAtomIdx()]
                    smarts_bonds.append(compound_mol.GetBondBetweenAtoms(aid1, aid2).GetIdx())

    atom_colors = {}
    for atom_index in smarts_atoms:
        if atom_index in shifts_atoms:
            atom_colors[atom_index] = (0, 0.7, 1)
        else:
            atom_colors[atom_index] = (1, 0, 0)

    bond_colors = {}
    for bond_index in smarts_bonds:
        if bond_index in shifts_bonds:
            bond_colors[bond_index] = (0, 0.7, 1)
        else:
            bond_colors[bond_index] = (1, 0, 0)

    # Find the atom indices and SMARTS indices for matching atoms in the structure
    ind_map = {}
    for index, atom in enumerate(smarts_mol.GetAtoms()):
        ind_map[atom.GetIdx()] = index
    # Add SMARTS index numbers to structure
    match_list = compound_mol.GetSubstructMatches(smarts_mol)
    for index, match in enumerate(match_list[0]):
        try:
            compound_mol.GetAtomWithIdx(match).SetProp("atomLabel", f"{index}")
        except:
            print("Error with SMARTS numbering")

    d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    rdMolDraw2D.PrepareAndDrawMolecule(d,
                                       compound_mol,
                                       highlightAtoms=smarts_atoms,
                                       highlightAtomColors=atom_colors,
                                       highlightBonds=smarts_bonds,
                                       highlightBondColors=bond_colors)
    d.FinishDrawing()
    d.GetDrawingText()
    # d = rdMolDraw2D.MolDraw2DCairo(500, 500)
    # rdMolDraw2D.PrepareAndDrawMolecule(d,
    #                                    compound_mol,
    #                                    highlightAtoms=smarts_atoms,
    #                                    highlightAtomColors=atom_colors,
    #                                    highlightBonds=smarts_bonds,
    #                                    highlightBondColors=bond_colors)
    # d.FinishDrawing()
    # d.GetDrawingText()
    output_path = Path(output_dir,
                       f"{compound_name}_{Chem.MolToSmiles(smarts_mol)}.png")
    d.WriteDrawingText(str(output_path))

    return compound_mol


if __name__ == "__main__":
    with open(ATLAS_PATH) as f:
        atlas_df = pd.read_csv(f, sep="\t")

    data_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas")
    # 14-member macrolide SMARTS string
    macrolide_smarts = "C1(=O)OCCCCCCCCCCCC1"
    generate_scaffold_match_image(atlas_df, macrolide_smarts, macrolide_smarts, data_path)
