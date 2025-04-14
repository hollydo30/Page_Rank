"""
COMP 614
Homework 6: DFS + PageRank
"""

import comp614_module6

def bfs_dfs(graph, start_node, rac_class):
    """
    Performs a breadth-first search on graph starting at the given node.
    Returns a two-element tuple containing a dictionary mapping each visited
    node to its distance from the start node and a dictionary mapping each
    visited node to its parent node.

    Inputs:
    graph - a graph class instance 
    start_node - a node
    rac_class - a reference to either the provided Queue class, or the provided Stack class 

    Outputs:
    parent - mapping of nodes to their parents

    """
    queue = rac_class()
    dist = {}
    parent = {}
    
    # Initialize distances and parents; no nodes have been visited yet
    for node in graph.nodes():
        dist[node] = float("inf")
        parent[node] = None

    # Initialize start node's distance to 0
    dist[start_node] = 0
    queue.push(start_node)

    # Continue as long as there are new reachable nodes
    while queue:
        node = queue.pop()
        nbrs = graph.get_neighbors(node)

        for nbr in nbrs:
            # Only update neighbors that have not been seen before
            if dist[nbr] == float('inf'):
                dist[nbr] = dist[node] + 1
                parent[nbr] = node
                queue.push(nbr)

    return parent


def recursive_dfs(graph, start_node, parent):
    """
    Given a graph, a start node from which to search, and a mapping of nodes to
    their parents, performs a recursive depth-first search on graph from the 
    given start node, populating the parents mapping as it goes.

    Inputs:
    graph - a graph class instance 
    start_node - a node from which to search 
    parent - mapping of nodes to parents 

    Outputs:
    None but under the hood, parent has been modified appropriately

    """
    nbrs = graph.get_neighbors(start_node)
    
    for nbr in nbrs:
        # base case
        if nbr in parent.keys():
            pass 
        else:
            parent[nbr] = start_node
            # recursive case
            recursive_dfs(graph, nbr, parent)
    #pass


def get_inbound_nbrs(graph):
    """
    Given a directed graph, returns a mapping of each node n in the graph to
    the set of nodes that have edges into n.

    Inputs:
    graph - a Graph class instance

    Outputs:
    dict_in - mapping of nodes to their set of inbound nodes

    """
    dict_in = {}
    nodes = graph.nodes()
    for node in nodes:
        dict_in[node] = set()

    # document inbounds of each node
    for node in nodes:
        for node_current in graph.nodes():
            for nbr in graph.get_neighbors(node_current):
                if (nbr == node):
                    dict_in[node].add(node_current)
    return dict_in


def get_outbound_nbrs(graph):
    """
    Given a directed graph, returns a mapping of each node n in the graph to
    the set of nodes that have edges out from n.

    Inputs:
    graph - a Graph class instance

    Outputs:
    dict_out - mapping of nodes to their set of outbound nodes

    """
    dict_out = {}

    # document outbounds of each node
    for node in graph.nodes():
        nbrs = graph.get_neighbors(node)
        dict_out[node] = nbrs

    return dict_out


def remove_sink_nodes(graph):
    """
    Given a directed graph, returns a new copy of the graph where every node that
    was a sink node in the original graph now has an outbound edge linking it to 
    every other node in the graph (excluding itself).

    Inputs:
    graph - a Graph class instance 

    Outputs:
    graph - a new copy of the graph, modified graph with added edges from sink nodes

    """

    mod_graph = graph.copy()

    parent_graph = get_outbound_nbrs(graph)

    # adding edges from sink nodes to every other nodes
    for node, set1 in parent_graph.items():
        if len(set1) == 0:
            for node2 in mod_graph.nodes():
                if node2 != node:
                    mod_graph.add_edge(node, node2)

    return mod_graph


def page_rank(graph, damping):
    """
    Given a directed graph and a damping factor, implements the PageRank algorithm
    -- continuing until delta is less than 10^-8 -- and returns a dictionary that 
    maps each node in the graph to its page rank.

    Inputs:
    graph - a Graph class instance
    damping - probability of following the last page rank image

    Outputs:
    page_rank - mapping from each node of the graph to it's rank number

    """
    page_rank1 = {}
    num = len(graph.nodes())
    delta = float('inf')

    # disable sink nodes by adding edges
    mod_graph = remove_sink_nodes(graph)

    # document inbound and outbound in dicts
    inbound = get_inbound_nbrs(mod_graph)
    outbound = get_outbound_nbrs(mod_graph)

    # init page rank to 1/num
    for node in mod_graph.nodes():
        page_rank1[node] = 1/num

    # iterate until delta is < 10^-8
    while delta >= 10**-8:
        page_rank_image = {}
        # page_rank k pi = (1-d)/n + d * sum inbounds pj:(page_rank k-1 (pj)/out degree pj)
        for node, _ in page_rank1.items():
            page_rank_image[node] = page_rank1[node] 
        for p_i in page_rank1:
            sum_p_j = 0
            for p_j in inbound[p_i]:
                sum_p_j += page_rank_image[p_j] / len (outbound[p_j])
            #print(sum_p_j)
            rank_p_i = (1-damping)/num + (damping * sum_p_j)

            page_rank1[p_i] = rank_p_i
        
        # delta is sum of differences
        delta = sum(abs(rank - page_rank_image[node]) for node, rank in page_rank1.items())

    return page_rank1
