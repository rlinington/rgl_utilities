#!usr/bin/env python3

import pandas as pd
from pathlib import Path
import shutil

from rdkit import Chem
from rdkit.Chem import Draw, AllChem


def create_example_node_images(min_compounds: int, max_compounds: int, atlas_path: Path, output_dir: Path):
    """Find all nodes with compound counts between min and max limits and create image of first compound in node

    Args:
        min_compounds (int): The minimum number of compounds in a node for it to be included
        max_compounds (int): The maximum number of compounds in a node for it to be included
        atlas_path (Path): The path to the atlas tsv download file
        output_dir (Path): The path to the output directory

    Returns:
        One image for each node within selected limits
        """

    atlas_df = pd.read_csv(atlas_path, sep='\t')

    for node_id, node_data in atlas_df.groupby("compound_node_id"):
        if min_compounds <= len(node_data.index) <= max_compounds:
            mol = Chem.MolFromSmiles(node_data['compound_smiles'].iloc[0])
            output_path = Path(output_dir,
                               f"node_{node_id:03d}_compound_{node_data['compound_id'].iloc[0]:06d}.png")
            Draw.MolToFile(mol, str(output_path))


def create_all_node_images(node_id: int, atlas_path: Path, output_dir: Path, substructure=None):
    """Create an image for each compound in a given node. If output dir exists then overwrite contents

    Args:
        node_id (int): the id number for the node of interest
        atlas_path (Path): The path to the atlas tsv download file
        output_dir (Path): The path to the output directory. Note that a node directory will be created at this
        destination and images stored in this node directory.
        substructure: An optional SMARTS string for aligning structures in set. If substructure is provided then an
        additional image is generated that includes all compounds in a standard orientation with respect to the SMARTS
        substructure

    Returns:
        One image for each node within selected limits
        """

    atlas_df = pd.read_csv(atlas_path, sep='\t')
    node_dir = Path(output_dir, f"node_{node_id:03d}")

    try:
        shutil.rmtree(str(node_dir))
    except OSError as e:
        print(f"ERROR: {e.strerror} for {e.filename}")

    Path.mkdir(node_dir, exist_ok=True)

    node_df = atlas_df[atlas_df['compound_node_id'] == node_id]
    node_compounds = list(zip(node_df['compound_id'], node_df['compound_smiles']))
    node_mols = []
    for node_compound in node_compounds:
        mol = Chem.MolFromSmiles(node_compound[1])
        mol.SetProp("_Name", f"NPAID{node_compound[0]:06d}")
        node_mols.append(mol)
        output_path = Path(node_dir, f"node_{node_id:03d}_compound_{node_compound[0]:06d}.png")
        Draw.MolToFile(mol, str(output_path))

    if substructure:
        p = Chem.MolFromSmarts(substructure)
        subms = [x for x in node_mols if x.HasSubstructMatch(p)]
        print(len(subms))
        AllChem.Compute2DCoords(p)
        for m in subms:
            _ = AllChem.GenerateDepictionMatching2DStructure(m, p)
        img = Draw.MolsToGridImage(subms, molsPerRow=4, subImgSize=(600, 600),
                                   legends=[x.GetProp("_Name") for x in subms])
        output_path = Path(node_dir, f"node_{node_id:03d}_compounds.png")
        img.save(str(output_path))
