import numpy as np 
import networkx as nx 

from networkx.drawing.nx_agraph import to_agraph
import os

# import PyBoolNet

import argparse

from neet.boolean import LogicNetwork
from neet import sensitivity

# At the moment, this can generate networks with cyclic attractors and networks with a single attractor

def build_network_structure(num_in=2,
    num_core=5,
    num_out=5,
    max_in_connections=3,
    num_output_nodes=2,
    seed=0):

    np.random.seed(seed)

    N = num_in + num_core + num_out

    in_nodes = ["n{}".format(i) for i in range(num_in)]
    core_nodes = ["n{}".format(i) for i in range(num_in, num_in+num_core)]
    out_nodes = ["n{}".format(i) for i in range(num_in+num_core, N)]

    graph = nx.DiGraph()

    graph.add_nodes_from(  in_nodes + core_nodes + out_nodes )


    ## add edges from in nodes to core

    for in_node in in_nodes:
        
        num_in_connections = np.random.choice(range(1, max_in_connections+1))

        for _ in range(num_in_connections):

            core_node = np.random.choice(core_nodes)

            graph.add_edge(in_node, core_node, weight=np.random.choice((-1, 1)))

    ## core edges

    while not nx.is_strongly_connected(graph.subgraph(core_nodes)):

        n1 = np.random.choice(core_nodes)

        n2 = np.random.choice(core_nodes + out_nodes)

        if n1 == n2:
            graph.add_edge(n1, n2, weight=1.)
        else:
            graph.add_edge(n1, n2, weight=np.random.choice((-1, 1)))

    num_cycles = 0
    for _ in nx.simple_cycles(graph.subgraph(core_nodes)):
        num_cycles += 1
        if num_cycles > 1:
            break

    assert num_cycles > 1

        # loop_size = np.random.choice(range(1, num_core))
        # loop_nodes = list(np.random.choice(core_nodes, loop_size, replace=False))

        # for n1, n2 in zip(loop_nodes, loop_nodes[1:] + loop_nodes[:1]):
        #     graph.add_edge(n1, n2, weight=np.random.choice((-1, 1)))


    # for n1, n2 in zip(core_nodes, core_nodes[1:] + core_nodes[:1]):
    #     graph.add_edge(n1, n2, weight=np.random.choice((-1, 1)))

    ## out nodes

    while not nx.is_weakly_connected(graph.subgraph(core_nodes + out_nodes[:-num_output_nodes])):

        n1 = np.random.choice(core_nodes)

        n2 = np.random.choice(out_nodes[:-num_output_nodes])

        graph.add_edge(n1, n2, weight=np.random.choice((-1, 1)))

    while not nx.is_weakly_connected(graph.subgraph(out_nodes)):

        n1 = np.random.choice(out_nodes[:-num_output_nodes])
        n2 = np.random.choice(out_nodes[-num_output_nodes:])

        graph.add_edge(n1, n2, weight=np.random.choice((-1, 1)))

    assert nx.is_weakly_connected(graph)

    return graph


def build_rules(
    graph,
    and_p=.5):

    # nodes = set(graph)

    # core_nodes = list(max(nx.strongly_connected_component_subgraphs(graph), key=len))

    # in_component = set([
    #     n 
    #     for n in nodes - set(core_nodes)
    #     if nx.has_path(graph, n, core_nodes[0])
    # ])

    # out_component = nodes - set(core_nodes) - in_component
    
    # print ("number of nodes in graph", len(graph))
    # print ("number of nodes in in_component", len(in_component))
    # print ("number of nodes in core", len(core_nodes))
    # print ("number of nodes in out_component", len(out_component))

    rules = {}

    for n in graph:

        in_nodes = list(graph.in_edges(n, data="weight"))

        if len(in_nodes) == 0:
            ## handle input nodes
            rule = n

        else:

            rule = ""

            for u, _, w in in_nodes:

                if rule != "":

                    connective = np.random.choice(["&", "|"], p=[and_p, 1-and_p])

                    rule += connective
                
                if w == -1:

                    u = "!" + u

                rule += u


        rules.update({n: rule})

    return rules


def write_bnet(rules, 
    filename):

    if not filename.endswith(".bnet"):
        filename += ".bnet"

    
    with open(filename, "w") as f:

        for k, v in rules.items():

            f.write("{},\t{}\n".format(k, v))

def write_neet_logic(rules,
    filename):

    if not filename.endswith(".txt"):
        filename += ".txt"

    with open(filename, "w") as f:

        for k, v in rules.items():

            v = v.replace("|", " OR ")
            v = v.replace("&", " AND ")
            v = v.replace("!", "NOT ")
            f.write("{} = {}\n".format(k, v))


def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(description="Generate bow tie networks with dynamics descrbed by .bnet file.")

    parser.add_argument("--num_in", 
        dest="num_in", type=int, default=2,
        help="Number of nodes in the in component.")
    parser.add_argument("--num_core", 
        dest="num_core", type=int, default=5,
        help="Number of nodes in the core.")
    parser.add_argument("--num_out", 
        dest="num_out", type=int, default=5,
        help="Number of nodes in the out component.")

    parser.add_argument("--max_in_connections", 
        dest="max_in_connections", type=int, default=3,
        help="Maximum number of edges from nodes in the in component.")
    parser.add_argument("--num_output_nodes", 
        dest="num_output_nodes", type=int, default=2,
        help="Number of nodes with no outgoing edges.")

    parser.add_argument("--root_directory", 
        dest="root_directory", type=str, 
        default="datasets/synthetic_bow_tie_networks",
        help="Directory to save networks.")

    parser.add_argument("--seed", 
        dest="seed", type=int, default=0,
        help="Random seed.")

    return parser.parse_args()

def main():

    args = parse_args()

    update = "synchronous"

    root_directory = args.root_directory
    if not os.path.exists(root_directory):
        os.mkdir(root_directory)

    # num_single_attractors = 0
    # num_single_steady_state_attractors = 0

    num_in = args.num_in
    num_core = args.num_core
    num_out = args.num_out

    max_in_connections = args.max_in_connections
    num_output_nodes = args.num_output_nodes

    num_seeds = 5

    for seed in range(num_seeds):

        output_directory = os.path.join(root_directory, 
            "{:03d}".format(seed))
        
        if not os.path.exists(output_directory):
            os.mkdir(output_directory)

        graph = build_network_structure(
            num_in,
            num_core,
            num_out,
            max_in_connections,
            num_output_nodes,
            seed=seed)

        nx.set_edge_attributes(graph, name="arrowhead",
		values={(u, v): ("normal" if w>0 else "tee") 
			for u, v, w in graph.edges(data="weight")})

        # draw whole network
        # plot_filename = os.path.join(output_directory, "whole_network.png")
        # graph.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
        # graph.graph['graph'] = {'scale': '3'}

        # a = to_agraph(graph)
        # a.layout('dot')   
        # a.draw(plot_filename)

        nx.write_edgelist(graph, os.path.join(output_directory, "edgelist.tsv"), data=["weight"], delimiter="\t")

        rules = build_rules(graph)

        bnet_filename = os.path.join(output_directory, "network.bnet")
        write_bnet(rules, bnet_filename)

        neet_logic_filename = os.path.join(output_directory, "network.txt")
        write_neet_logic(rules, neet_logic_filename)


        # net = LogicNetwork.read_logic(neet_logic_filename)

        # print ("read neet logic file:", neet_logic_filename)
        # avg_sensitivity = sensitivity.average_sensitivity(net)
        # print ("average sensitivity", avg_sensitivity)




        # print ("reading", bnet_filename)
        # primes = PyBoolNet.FileExchange.bnet2primes(bnet_filename)  

        # attrs_original = PyBoolNet.Attractors.compute_json(primes, 
        #     update, 
        #     Silent=True, )

        # if len(attrs_original["attractors"]) == 1:
        #     num_single_attractors += 1
        #     if attrs_original["attractors"][0]["is_steady"]:
        #         num_single_steady_state_attractors += 1

        # print ("competed seed {:03d}, number of attractors: {}".format(seed, len(attrs_original["attractors"])))

    # print ("proportion of num single attractors:", num_single_attractors/num_seeds)
    # print ("proportion of num single steady state attractors:", num_single_steady_state_attractors/num_seeds)


if __name__ == "__main__":
    main()