#!/usr/bin/env python3

"""Tools to modify the structure of graphML files"""


import networkx as nx
import os


def filter_by_activity(graph, min_activity_score, min_cluster_score):
    """Create filtered graph from comprehensive graph from NP Analyst by removing all mz nodes with activity or cluster
    scores below defined cutoffs."""

    # Remove baskets with low activity or cluster scores
    drop_basket_list = []

    for node, data in graph.nodes(data=True):
        if data["type_"] == "basket":
            if data["cluster_score"] < min_cluster_score or data["activity_score"] < min_activity_score:
                drop_basket_list.append(node)

    for node in drop_basket_list:
        graph.remove_node(node)

    # Remove sample nodes that have no attached baskets
    drop_sample_list = []

    for node, data in graph.nodes(data=True):
        if data["type_"] == "sample" and len(graph.edges(node)) == 0:
            drop_sample_list.append(node)

    for node in drop_sample_list:
        graph.remove_node(node)

    # Save filtered graph
    with open(os.path.expanduser(os.path.join("~",
                                              "Git",
                                              "rgl_utilities",
                                              "data",
                                              "network_management",
                                              "np_analyst",
                                              "network_2_03_filtered.graphml")), "wb") as g:
        nx.write_graphml(graph, g)


if __name__ == "__main__":

    with open(os.path.expanduser(os.path.join("~",
                                              "Git",
                                              "rgl_utilities",
                                              "data",
                                              "network_management",
                                              "np_analyst",
                                              "network_2_03_communities.graphml"))) as f:
        input_graph = nx.read_graphml(f)

    filter_by_activity(input_graph, 2, 0.3)