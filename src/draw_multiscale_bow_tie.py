import os
import argparse

import numpy as np 
import networkx as nx 
import pandas as pd 

from networkx.drawing.nx_agraph import to_agraph

import glob

def get_one_hop_neighbour_subgraph(g, nodes):

	if not isinstance(nodes, set):
		nodes = set(nodes)

	g_ = g.to_undirected()

	neighbours = set().union(*(g_.neighbors(n)
		for n in nodes)) - nodes

	return nx.DiGraph(g.subgraph(nodes.union(neighbours)))

def draw_bow_tie(g, 
	core_of_core,
	output_dir, 
	filename):

	assert isinstance(core_of_core, set)

	# draw all SCCs to file (to use as images on nodes)
	for i, g_ in enumerate(\
		nx.strongly_connected_component_subgraphs(g)):
				
		plot_filename = os.path.join(output_dir,
			"scc_{}.png".format(i))
		g_.graph['edge'] = {'arrowsize': '.8', 
			'splines': 'curved'}
		g_.graph['graph'] = {'scale': '3'}

		a = to_agraph(g_)

		# # add core of core subgraph
		# if all([c in g_ for c in core_of_core]):
		# 	a.add_subgraph(core_of_core, 
		# 		name="cluster_core",
		# 		color="black",
		# 		bgcolor="red")

		a.layout('dot')   
		a.draw(plot_filename)

		print ("plotted", plot_filename)

	# use nx.DiGraph to avoid multiple edges between nodes
	h = nx.MultiDiGraph(g.copy())

	# collapse all SCCs into single metanodes
	map_ = {}
	for i, scc in enumerate(nx.strongly_connected_components(g)):
		if len(scc) == 1:
			h.node[list(scc)[0]]["image"] = \
				os.path.join(output_dir, 
				"scc_{}.png".format(i))
		else:
			m = frozenset(scc)
			h.add_node(m, 
			image=os.path.join(output_dir,
				"scc_{}.png".format(i)))

			for n in scc:
				map_.update({n: m})
				for u, _, w in \
				filter(lambda x: x[0] not in scc, 
				g.in_edges(n, data="weight")):
					if u in map_:
						u = map_[u]
					assert u in h, u
					h.add_edge(u, m, 
					weight=w,
					arrowhead=("normal" if w>0 else "tee"),
					label="to_{}".format(n))
				for _, v, w in \
				filter(lambda x: x[1] not in scc, 
				g.out_edges(n, data="weight")):
					if v in map_:
						v = map_[v]
					assert v in h
					h.add_edge(m, v, 
					weight=w,
					arrowhead=("normal" if w>0 else "tee"),
					label="from_{}".format(n))
				h.remove_node(n)

	# remove self loops to clean up drawing
	print ("removing self loop edges")
	h.remove_edges_from(list(nx.selfloop_edges(h))) 

	# remove label from nodes so the images can be seen
	nx.set_node_attributes(h, name="label", values="")

	plot_filename = os.path.join(output_dir,
		filename)
	h.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
	h.graph['graph'] = {'scale': '3'}

	a = to_agraph(h)

	a.layout('dot')   
	a.draw(plot_filename)

	print ("plotted", plot_filename)


	## cleanup of directory
	print ("cleaning up directory")
	for f in glob.iglob(os.path.join(output_dir, "scc_*.png")):
		print ("removing", f)
		os.remove(f)


def parse_args():
	'''
	Parse from command line
	'''
	parser = argparse.ArgumentParser(description="Read in edgelist and write bow-tie structure to file")

	parser.add_argument("--edgelist", 
		dest="edgelist", type=str, default=None,
		help="edgelist to load.")
	parser.add_argument("--mapping", 
		dest="mapping", type=str, default=None,
		help="mapping file of node ids to names.")

	parser.add_argument("--merge_depths", 
		dest="merge_depths", type=str, 
		help="Path to merge depth file for this network.")

	parser.add_argument("--output", dest="output", 
		type=str, default=None,
		help="Directory to save images.")

	parser.add_argument("--max-hops", 
		dest="max_hops", type=int, default=2, 
		help="Maximum number of hops to consider.")

	return parser.parse_args()

def main():
	args = parse_args()

	edgelist_file = args.edgelist

	print ("decomposing", edgelist_file)

	g = nx.read_weighted_edgelist(edgelist_file, 
		create_using=nx.DiGraph(), 
		delimiter="\t")

	mapping_file = args.mapping
	if mapping_file:
		print ("relabeling nodes using", mapping_file)
		mapping = pd.read_csv(mapping_file, 
			index_col=0, header=None, dtype=str)[1].to_dict()
		mapping = {str(k): v for k, v in mapping.items()}
		g = nx.relabel_nodes(g, mapping=mapping)

	print ("found graph with", len(g), 
		"nodes and", len(g.edges()), "edges")

	nx.set_edge_attributes(g, name="arrowhead",
		values={(u, v): ("normal" if w>0 else "tee") 
			for u, v, w in g.edges(data="weight")})

	merge_depths_filename = args.merge_depths
	print ("reading merge depths from", 
		merge_depths_filename)
	merge_depth_df = pd.read_csv(merge_depths_filename, 
		index_col=0, header=None )
	max_merge_depth = merge_depth_df[1].max()
	core_of_core = set(merge_depth_df.index[
		merge_depth_df[1] == max_merge_depth].tolist())

	print ("core of core genes are", core_of_core)

	for gene in core_of_core:
		assert gene in g

	output_dir = args.output
	if not os.path.exists(output_dir):
		print ("making directory", output_dir)
		os.makedirs(output_dir, exist_ok=True)

	nodes = core_of_core

	for i in range(1, args.max_hops+1):

		print ("getting", i, "hop subgraph")
		subgraph = get_one_hop_neighbour_subgraph(g, nodes)

		# nx.set_node_attributes(subgraph, name="color",
		# 	values={n: "red" for n in core_of_core})

		nx.set_node_attributes(subgraph, name="penwidth",
			values={n: 5.0 for n in core_of_core})

		print ("found subgraph with", len(subgraph), "nodes")

		print ("plotting")
		draw_bow_tie(subgraph, 
			core_of_core,
			output_dir, 
			"bow_tie_{}_hop.png".format(i))

		# update nodelist to current nodes in subgraph
		nodes = subgraph
	
if __name__ == "__main__":
	main()