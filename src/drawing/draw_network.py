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

	# draw whole network
	plot_filename = os.path.join(output_dir,
		"whole_network.png")
	g.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
	g.graph['graph'] = {'scale': '3'}

	a = to_agraph(g)
	a.layout('dot')   
	a.draw(plot_filename)

if __name__ == "__main__":
    main()