import os
import argparse

import numpy as np 
import networkx as nx 
import pandas as pd 

from networkx.drawing.nx_agraph import to_agraph

import glob

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
	parser.add_argument("--output", dest="output", 
		type=str, default=None,
		help="Directory to save images.")

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

	output_dir = args.output
	if not os.path.exists(output_dir):
		os.mkdir(output_dir)

	# # draw whole network
	# plot_filename = os.path.join(output_dir,
	# 	"whole_network.png")
	# g.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
	# g.graph['graph'] = {'scale': '3'}

	# a = to_agraph(g)
	# a.layout('dot')   
	# a.draw(plot_filename)

	# draw whole network in bow-tie form 
	h = g.copy()
	plot_filename = os.path.join(output_dir,
		"bow_tie_with_edges.png")
	h.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
	h.graph['graph'] = {'scale': '3'}

	a = to_agraph(h)

	for i, scc in enumerate(filter(lambda x: len(x) > 1, 
		nx.strongly_connected_components(h))):
		a.add_subgraph(scc, 
			name="cluster_{}".format(i)) 

	a.layout('dot')   
	a.draw(plot_filename)

	print ("plotted", plot_filename)

	# draw all SCCs to file (to use as images on nodes)
	for i, g_ in enumerate(\
		nx.strongly_connected_components(g)):
		g_ = g.subgraph(g_)

		plot_filename = os.path.join(output_dir,
			"scc_{}.png".format(i))
		g_.graph['edge'] = {'arrowsize': '.8', 
			'splines': 'curved'}
		g_.graph['graph'] = {'scale': '3'}

		a = to_agraph(g_)
		a.layout('dot')   
		a.draw(plot_filename)

		print ("plotted", plot_filename)

	# use nx.DiGraph to avoid multiple edges between nodes
	h = nx.MultiDiGraph(g.copy())

	# collapse all SCCs into single metanodes
	map_ = {}
	for i, scc in enumerate(nx.strongly_connected_components(g)):
		if len(scc) == 1:
			h.nodes[list(scc)[0]]["image"] = \
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
					arrowhead=("normal" if w>0 else "tee"))
				for _, v, w in \
				filter(lambda x: x[1] not in scc, 
				g.out_edges(n, data="weight")):
					if v in map_:
						v = map_[v]
					assert v in h
					h.add_edge(m, v, 
					weight=w,
					arrowhead=("normal" if w>0 else "tee"))
				h.remove_node(n)

	# remove self loops to clean up drawing
	print ("removing self loop edges")
	h.remove_edges_from(list(nx.selfloop_edges(h))) 

	# remove label from nodes so the images can be seen
	nx.set_node_attributes(h, name="label", values="")

	plot_filename = os.path.join(output_dir,
		"bow_tie.png")
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

if __name__ == "__main__":
	main()