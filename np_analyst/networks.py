#!usr/bin/env python3

"""Tools to reprocess output data from NP Analyst 1.0"""

import pandas as pd
import networkx as nx
from pathlib import Path


def table_to_network(table_path: Path, activity_score: float, cluster_score: float, output_dir: Path):
    """Reprocess output csv table from NP Analyst and create new network with specified activity and cluster scores.
    Generates graphML file with the same attributes and column headers as the download graph from the online system
    Args:
        table_path (Path): Path to NP Analyst csv download
        activity_score (float): Desired Activity Score cutoff for network
        cluster_score (float): Desired Cluster Score cutoff for network
        output_dir (Path): Path to output directory"""

    table_df = pd.read_csv(table_path)
    filtered_table_df = table_df[(table_df["ACTIVITY_SCORE"] >= activity_score) & (table_df["CLUSTER_SCORE"] >= cluster_score)]
    sample_names = filtered_table_df["SampleList"].to_list()
    sample_name_list = []
    for sample_name in sample_names:
        name_list = sample_name.split("|")
        sample_name_list = sample_name_list + name_list
    unique_names_list = list(set(sample_name_list))

    npanalyst_graph = nx.Graph()
    npanalyst_graph.add_nodes_from(unique_names_list)
    node_attributes = {}
    for node in npanalyst_graph.nodes(data=True):
        node_attributes[node[0]] = {"id": node[0],
                                    "type_": "sample",
                                    "x": "",
                                    "y": "",
                                    "community": ""}
    nx.set_node_attributes(npanalyst_graph, node_attributes)
    npanalyst_graph.add_nodes_from(filtered_table_df["BasketID"].to_list())
    filtered_list = filtered_table_df.values.tolist()
    node_attributes = {}
    for row in filtered_list:
        if row[9] > 0:
            community = int(row[9])
            for sample in row[6].split("|"):
                node_attributes[sample] = {"community": community}
        else:
            community = ""
        node_attributes[row[0]] = {"id": row[0],
                                   "type_": "basket",
                                   "x": "",
                                   "y": "",
                                   "community": community,
                                   "freq": len(row[6].split("|")),
                                   "samples": row[6],
                                   "PrecMz": row[1],
                                   "RetTime": row[5],
                                   "PrecIntensity": row[2],
                                   "MinPrecIntensity": row[3],
                                   "MaxPrecIntensity": row[4],
                                   "activity_score": row[7],
                                   "cluster_score": row[8],}
        for sample in row[6].split("|"):
            npanalyst_graph.add_edge(row[0], sample)
    nx.set_node_attributes(npanalyst_graph, node_attributes)

    output_path = Path(output_dir, f"{table_path.stem}_activity_{activity_score}_cluster_{cluster_score}.graphML")
    nx.write_graphml(npanalyst_graph, str(output_path))


if __name__ == "__main__":
    table_path = Path("/Users/roger/Git/rgl_utilities/data/npanalyst/csv_tables/npanalyst_1ad81e2e-e2cc-4805-a0a3-fea86ad0af2d.csv")
    output_dir_path = Path("/Users/roger/Git/rgl_utilities/data/npanalyst/output")
    activity_score = 2
    cluster_score = -0.01

    table_to_network(table_path, activity_score, cluster_score, output_dir_path)