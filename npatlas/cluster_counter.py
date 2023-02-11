#!/usr/bin/env python3

"""Tools to extract statistics on cluster frequencies, distributions etc. from NP Atlas tsv download"""

import pandas as pd
from pathlib import Path
import sys


def count_clusters(atlas_df, db_type="all", compound_percentage=90):
    """Return number of clusters that account for the given percentage of compounds. Used to say (for example) that
    '90% of compounds are described by just 840 clusters'

    db_type can be 'all', 'bacterial' or 'fungal'
    compound_percentage can be from 0 - 100"""

    if db_type == "all":
        filtered_df = atlas_df
    elif db_type in ("Bacterium", "Fungus"):
        filtered_df = atlas_df[atlas_df["origin_type"] == db_type]
    else:
        print("ERROR. db_type must be 'all', 'Bacterium', or 'Fungus'")
        sys.exit()

    total_compound_count = len(filtered_df.index)
    print("Total compound count = " + str(total_compound_count))
    cluster_counts = filtered_df.groupby(['compound_cluster_id']).size().reset_index(name='counts')
    sorted_cluster_counts = cluster_counts.sort_values(by="counts", ascending=False)["counts"].to_list()

    compound_count = 0
    for i, count in enumerate(sorted_cluster_counts):
        compound_count += count
        if i == 999:
            print(str(round(compound_count * 100/total_compound_count, 1)) + "% of compounds contained in 1000 clusters")
        if compound_count > total_compound_count * (compound_percentage/100):
            print(str(compound_percentage) + "% of compounds contained in " + str(i + 1) + " clusters")
            break


if __name__ == "__main__":

    db_path = Path("/Users/roger/Git/rgl_utilities/data/NP_Atlas/NPAtlas_2022_09.tsv")

    input_df = pd.read_csv(db_path, delimiter='\t')
    count_clusters(input_df)
    count_clusters(input_df, db_type='Bacterium')
    count_clusters(input_df, db_type='Fungus')