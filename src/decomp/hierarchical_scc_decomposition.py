import os
import argparse

import numpy as np 
import networkx as nx 
import pandas as pd 
import itertools

from networkx.drawing.nx_agraph import write_dot, graphviz_layout, to_agraph

import matplotlib.pyplot as plt

import glob

def find_cycles(u, n, g, start, l=set()):   

	if isinstance(g, nx.MultiDiGraph):
		g = nx.DiGraph(g)
	
	if n==0:
		assert u == start
		return [[u]]
	
	l_ = l .union( {u} )
	if n > 1:
		neighbors = set(g.neighbors(u)) - l_
	else:
		neighbors = set(g.neighbors(u)) . intersection({start})
		
	paths = ( [u] + cycle
		for neighbor in neighbors
		for cycle in find_cycles(neighbor, n-1, g, start, l_) )
	return paths

def score_subgraph_module(g, groups):
	subgraph = g.subgraph(groups)
	# n = len(subgraph)

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v])

	k_self = len([(u, v) for 
		u, v, w in subgraph.edges(data="weight") 
		if u == v])

	k_all = sum((len(list(g.neighbors(u))) for u in subgraph))

	return (k_in + k_self) / k_all

def score_subgraph_module_positive(g, groups):
	subgraph = g.subgraph(groups)
	# n = len(subgraph)
	# assert n < 4

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v and w > 0])

	k_self = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u == v and w>0])
	k_self = 0

	k_all = sum((len(list(g.neighbors(u))) for u in subgraph))

	return (k_in + k_self) / k_all

def score_subgraph_density(g, groups):
	
	'''
	score a subgraph of g given by groups
	'''

	subgraph = g.subgraph(groups)
	n = len(subgraph)

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v in subgraph.edges() 
		if u != v])

	k_self = len([(u, v) 
		for u, v in subgraph.edges() 
		if u == v])

	return  (k_in + k_self) / n ** 2

def score_subgraph_density_positive(g, groups):
	
	'''
	score a subgraph of g given by groups
	'''

	subgraph = g.subgraph(groups)
	n = len(subgraph)

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v and w > 0])

	k_self = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u == v and w > 0])

	return  (k_in + k_self) / n ** 2

def score_subgraph_num_loops(g, groups):
	subgraph = g.subgraph(groups)
	n = len(subgraph)
	# assert n < 4
	d = (1 if n==1
		else 3 if n == 2
		else 8)
	return len(list(nx.simple_cycles(subgraph))) / d

def bottom_up_partition(g, 
	score_function=score_subgraph_density,
	subgraph_sizes=(2, 3)):

	'''
	perform partition in bottom-up manner
	'''

	g = nx.MultiDiGraph(g.copy())
	# g = nx.MultiGraph(g.copy())
	
	h = nx.DiGraph()

	h.add_nodes_from( g.nodes() )

	# ## handle self loops
	for u, _ in nx.selfloop_edges(g):
		h.add_edge(u, frozenset([u]))
	g = nx.relabel_nodes(g, 
		mapping={n: frozenset([n]) 
		for n, _ in nx.selfloop_edges(g)})
	
	subgraph_scores = {}

	i = 0
	for s in subgraph_sizes:
		print (s)
		for n in g.nodes():
			for cycle in map(frozenset, 
				find_cycles(n, s, g, start=n)):
				if cycle in subgraph_scores:
					continue
				subgraph_scores[cycle] = \
					score_function(g, cycle)
			print ("processed", i, "/", 
				len(g) * len(subgraph_sizes))
			i += 1

	# repeat until g has collapsed into a single node    
	while len(g) > 1:

		print ("number of nodes in g:", len(g), )


		if len(subgraph_scores) > 0:

			# determine all highest scoring subgraphs
			sorted_subgraphs = sorted(subgraph_scores, 
				key=lambda x: subgraph_scores[x],#(subgraph_scores[x], len(x)),
				reverse=True)
			chosen_subgraph = sorted_subgraphs.pop(0)
			chosen_subgraph_score = subgraph_scores[chosen_subgraph]
			# assert chosen_subgraph_score > 0, chosen_subgraph_score
			if chosen_subgraph_score > 0:

				chosen_subgraphs = [chosen_subgraph]

				for subgraph in sorted_subgraphs:
					if subgraph_scores[subgraph] < chosen_subgraph_score:
						break
					chosen_subgraphs.append(subgraph)

				# combine any chosen subgraphs that contain the 
				# same nodes
				if len(chosen_subgraphs) > 1:

					print ("number of chosen subgraphs before collapse", 
						len(chosen_subgraphs))

					overlaps = np.array([[len(x.intersection(y)) 
						for y in chosen_subgraphs]
						for x in chosen_subgraphs])
					overlap_g = nx.Graph(overlaps)
					chosen_subgraphs = [frozenset().union([x 
						for c in cc for x in chosen_subgraphs[c]]) 
						for cc in nx.connected_components(overlap_g)]

					print ("number of chosen subgraphs after collapse", 
						len(chosen_subgraphs))
			else:
				# could not find a subgraph with positive score
				print ("no subgraphs with positive score", chosen_subgraph_score)
				chosen_subgraphs = [frozenset().union([x for x in g])]

		else:
			# could not find a subgraph of selected sizes
			print ("no subgraphs of sizes", subgraph_sizes)
			chosen_subgraphs = [frozenset().union([x for x in g])]

		for chosen_subgraph in chosen_subgraphs:

			# remove scores associated with merged nodes
			print ("removing scores associated with", 
				"chosen subgraph")
			subgraph_scores = {k: v 
				for k, v in subgraph_scores.items()
				if not any([x in k for x in chosen_subgraph])}
			print ("done")
			
			# merge subgraph into super-node
			g.add_node(chosen_subgraph)

			for n in chosen_subgraph:

				new_edges = []

				if isinstance(g, nx.DiGraph):
					for u, _, w in g.in_edges(n, data="weight"):
						# if u == chosen_subgraph:
						# 	continue
						if u in chosen_subgraph:
							u = chosen_subgraph
						new_edges.append(
							(u, chosen_subgraph, {"weight": w}))
						# g.add_edge(u, chosen_subgraph, 
							# weight=w)
					for _, v, w in g.out_edges(n, data="weight"):
						# if v == chosen_subgraph:
						# 	continue
						if v in chosen_subgraph:
							v = chosen_subgraph
						new_edges.append(
							(chosen_subgraph, v, {"weight": w}))
						# g.add_edge(chosen_subgraph, v,
							# weight=w)
				else:
					
					# undirected case
					# assert False 
					for _, v, w in g.edges(n, data="weight"):
						# if v == chosen_subgraph:
							# continue
						if v in chosen_subgraph:
							v = chosen_subgraph
						new_edges.append(
							(chosen_subgraph, v, {"weight": w}))
						# g.add_edge(chosen_subgraph, v, 
							# weight=w)

				g.add_edges_from(new_edges)
				g.remove_node(n)

			# add chosen subgraph to h
			h.add_node(chosen_subgraph)
			for n in chosen_subgraph:
				h.add_edge(n, chosen_subgraph)

			# add cycles containing new node
			print ("determining cycles containing new node")
			for s in subgraph_sizes:
				for cycle in map(frozenset, 
				find_cycles(chosen_subgraph, s, 
					g, start=chosen_subgraph)):
					# assert nx.is_strongly_connected(g.subgraph(cycle))
					if cycle in subgraph_scores:
						assert False
						continue
					subgraph_scores[cycle] = \
						score_function(g, cycle)
			print ("done")

	return h

def decompose_all_sccs(g, 
	score_function=score_subgraph_density,
	subgraph_sizes=(2, 3)):
	'''
	run decomposition on each SCC in g
	'''
	h = nx.DiGraph()
	roots = []
	for scc in nx.strongly_connected_components(g):
		scc = g.subgraph(scc)
		scc_tree = bottom_up_partition(scc, 
			score_function=score_function,
			subgraph_sizes=subgraph_sizes)
		degrees = dict(scc_tree.out_degree())
		root = min(degrees, key=degrees.get)
		roots.append(root)
		h = nx.union(h, scc_tree)

	# add final root to represent whole network
	all_nodes = frozenset(g.nodes())
	if len(roots) > 1:
		for root in roots:
			h.add_edge(root, all_nodes)

	return h

def unpack(x):
	if isinstance(x, str):
		return [x]
	if not any([isinstance(x_, frozenset) for x_ in x]):
		return list(x)
	else:
		return [_x for x_ in x for _x in unpack(x_)]

def parse_args():
	'''
	Parse from command line
	'''
	parser = argparse.ArgumentParser(description="Read in edgelist and draw SCC decomposition")

	parser.add_argument("--edgelist", 
		dest="edgelist", type=str, default=None,
		help="edgelist to load.")
	parser.add_argument("--mapping", 
		dest="mapping", type=str, default=None,
		help="mapping file of node ids to names.")
	parser.add_argument("--output", dest="output", 
		type=str, default=None,
		help="Directory to save images/merge depths.")
	parser.add_argument("--score-function",
		dest="score_function", 
		type=str, default="density_pos", 
		choices=["density", "density_pos", 
			"module", "module_pos",
			"num_loops"],
		help="Scoring function.")
	parser.add_argument("--draw", action="store_true",
		help="Flag to specify to plot or not.")

	return parser.parse_args()

def main():

	args = parse_args()
	
	edgelist_file = args.edgelist
	
	print ("decomposing", edgelist_file)

	g = nx.read_weighted_edgelist(edgelist_file, 
		create_using=nx.DiGraph(), 
		delimiter="\t")

	mapping_file = args.mapping
	if mapping_file is not None:
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

	score_function = args.score_function
	print ("using", score_function, "score function")

	output_dir = args.output
	output_dir = os.path.join(output_dir, score_function)
	if not os.path.exists(output_dir):
		os.makedirs(output_dir, exist_ok=True)

	if score_function ==  "density":
		score_function = score_subgraph_density
	elif score_function == "density_pos":
		score_function = score_subgraph_density_positive
	elif score_function == "module":
		score_function = score_subgraph_module
	elif score_function == "module_pos":
		score_function = score_subgraph_module_positive
	else:
		score_function = score_subgraph_num_loops

	h = decompose_all_sccs(g, 
		score_function=score_function,
		subgraph_sizes=[2,3,] )
	
	if args.draw:

		h = h.reverse()

		map_ = {}
		for i, n in enumerate(h.nodes()):
			if isinstance(n, frozenset):
				g_ = nx.MultiDiGraph(g.subgraph(unpack(n)))
			else:
				g_ = nx.MultiDiGraph(g.subgraph([n]))
				g_.remove_edges_from(list(nx.selfloop_edges(g_)))
			nx.set_edge_attributes(g_, name="arrowhead",
				values={(u, v, w): ("normal" if w>0 else "tee") 
					for u, v, w in g_.edges(data="weight")})

			map_.update({n: i})
			
			if len(g_) > 1:
				for no, child in enumerate(h.neighbors(n)):
					# make metanode
					node = "metanode_{}".format(no)
					g_.add_node(node, label="", 
					image=os.path.join(output_dir, 
						"subgraph_{}.png".format(map_[child])))
					for n_ in unpack(child):
						for u, _, w in g_.in_edges(n_, data="weight"):
							if u == node:
								continue
							g_.add_edge(u, node, 
							weight=w, 
							arrowhead=("normal" if w>0 else "tee"))
						for _, v, w in g_.out_edges(n_, data="weight"):
							if v == node:
								continue
							g_.add_edge(node, v, 
							weight=w,
							arrowhead=("normal" if w>0 else "tee"))
						g_.remove_node(n_)
					
			plot_filename = os.path.join(output_dir,
				"subgraph_{}.png".format(i))
			g_.graph['edge'] = {'arrowsize': '.8', 
				'splines': 'curved'}
			g_.graph['graph'] = {'scale': '3'}

			a = to_agraph(g_)
			a.layout('dot')   
			a.draw(plot_filename)

			print ("plotted", plot_filename)
		
		h = h.reverse()

		nx.set_node_attributes(h, name="image", 
			values={n: 
			os.path.join(output_dir, 
			"subgraph_{}.png".format(i) )
			for i, n in enumerate(h.nodes())})
		nx.set_node_attributes(h, name="label", values="")

		tree_plot_filename = os.path.join(output_dir, 
			"scc_tree.png")
		h.graph['edge'] = {'arrowsize': '.8', 'splines': 'curved'}
		h.graph['graph'] = {'scale': '3'}

		a = to_agraph(h)
		a.layout('dot')   
		a.draw(tree_plot_filename)

		print ("plotted", tree_plot_filename)

	print ("determining merge depths")

	out_degrees = dict(h.out_degree())
	root = min(out_degrees, key=out_degrees.get)

	merge_depths = {node: \
		nx.shortest_path_length(h, node, root) 
		for node in g }

	merge_depth_filename = os.path.join(output_dir,
		"merge_depths.csv")
	print ("saving merge depths to", merge_depth_filename)
	pd.Series(merge_depths).to_csv(merge_depth_filename)

	max_merge_depth = max(merge_depths.values())
	genes_with_max_merge_depth = [k 
		for k, v in merge_depths.items()
		if v == max_merge_depth]
	print ("genes with maximum merge depth")
	print (genes_with_max_merge_depth)
	print ()


	# target input genes
	core = list(max(nx.strongly_connected_components(g),
		key=len))
	
	in_component = [
		n for n in g
		if n not in core
		and nx.has_path(g, n, core[0])
	]

	scores = {u: 
		np.mean([1. / nx.shortest_path_length(g, u, c) 
			for c in genes_with_max_merge_depth])
		for u in in_component
	}

	print ()
	print ("In-component scores:")
	for u in sorted(scores, key=scores.get, reverse=True):
		print ("{:6s}\t{:.03f}".format(u, scores[u]))

	score_filename = os.path.join(output_dir, 
		"in_component_scores.csv")
	print ("saving scores to", score_filename)
	pd.DataFrame.from_dict(scores, 
		orient="index").to_csv(score_filename)


	## cleanup of directory
	print ("cleaning up directory")
	for f in glob.iglob(os.path.join(output_dir, "subgraph_*.png")):
		print ("removing", f)
		os.remove(f)

if __name__ == "__main__":
	main()