#!usr/bin/env python3

"""Tools to extract and group biosynthetic classes based on structural and biosynthetic attributes"""

from pathlib import Path
from npatlas import biosynthetic_classes, scaffold_matching, atlas_to_df

ATLAS_JSON_PATH = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.json")
ATLAS_PATH = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_download_2024_09.tsv")
DATA_PATH = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/Macrolides")

if __name__ == "__main__":

    atlas_df = atlas_to_df.convert_tsv(ATLAS_PATH)
    # atlas_df = biosynthetic_classes.convert_json(ATLAS_JSON_PATH)
    # biosynthetic_class_df = biosynthetic_classes.extract_biosynthetic_class(atlas_df,
    #                                                                         "pathway_results",
    #                                                                         "Polyketides")

    # for i in range(14):
    #     smarts = "C1(=O)OCCCCCCC" + "CC" * i + "C1"
    #     macrolide_size = 10 + 2*i
    #     print(f"Macrolide size: {macrolide_size}. SMARTS: {smarts}")
    #     scaffold_matching.generate_scaffold_match_image(atlas_df,
    #                                                     smarts,
    #                                                     f"{macrolide_size}_membered_ester_pks",
    #                                                     DATA_PATH)

    target_macrolide = "14_membered_ester_pks"

    target_df = atlas_to_df.convert_tsv(Path(DATA_PATH, "Atlas_SMARTS_matches_" + target_macrolide + ".tsv"))
    scaffold_matching.annotate_pks_substitution(target_df,
                                                "C1(=O)OCCCCCCCCCCCC1",
                                                target_macrolide,
                                                DATA_PATH)