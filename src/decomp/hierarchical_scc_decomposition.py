import os
import argparse

import numpy as np 
import networkx as nx 
import pandas as pd 
import itertools

from networkx.drawing.nx_agraph import write_dot, graphviz_layout, to_agraph

import matplotlib.pyplot as plt

import glob

def find_chains(u, n, g, l=set()):   

	if n==0:
		return [[u]]
	
	l_ = l .union( {u} )
	# if n > 1:
	neighbors = set(g.neighbors(u)) - l_
	# else:
	# 	neighbors = set(g.neighbors(u)) . intersection({start})
		
	paths = ( [u] + chain
		for neighbor in neighbors
		for chain in find_chains(neighbor, n-1, g, l_) )
	return paths

def find_cycles(u, n, g, start, l=set()):   
	
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

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v])

	k_self = len([(u, v) for 
		u, v, w in subgraph.edges(data="weight") 
		if u == v])

	k_all = sum((len(list(g.neighbors(u))) for u in subgraph))

	return (k_in + 0) / k_all

def score_subgraph_module_positive(g, groups):
	subgraph = g.subgraph(groups)

	# k_neg =  len([(u, v) 
	# 	for u, v, w in subgraph.edges(data="weight") 
	# 	if u != v and w < 0])
	# if k_neg > 0:
	# 	return 0

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v and w > 0])

	k_self = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u == v and w > 0])

	k_all = sum((len(list(g.neighbors(u))) 
		for u in subgraph))

	return (k_in + 0) / k_all

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

	return  (k_in + 0) / n ** 2

def score_subgraph_density_positive(g, groups):
	
	'''
	score a subgraph of g given by groups
	'''

	subgraph = g.subgraph(groups)
	n = len(subgraph)

	# k_neg =  len([(u, v) 
	# 	for u, v, w in subgraph.edges(data="weight") 
	# 	if u != v and w < 0])
	# if k_neg > 0:
	# 	return 0

	# all internal edges of subgraph
	k_in = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u != v and w > 0])

	k_self = len([(u, v) 
		for u, v, w in subgraph.edges(data="weight") 
		if u == v and w > 0])

	return  (k_in + 0) / n ** 2

def score_subgraph_bc(g, groups):
	
	bc = nx.betweenness_centrality(nx.DiGraph(g))

	return np.mean([bc[n] for n in groups])

def collapse_subgraph(g, subgraph):

	g.add_node(subgraph)

	for n in subgraph:

		new_edges = []

		if isinstance(g, nx.DiGraph):
			for u, _, w in g.in_edges(n, data="weight"):
				if u in subgraph:
					u = subgraph
				new_edges.append(
					(u, subgraph, {"weight": w}))
			for _, v, w in g.out_edges(n, data="weight"):
				if v == n: #  self loops already included
					continue
				if v in subgraph:
					v = subgraph
				new_edges.append(
					(subgraph, v, {"weight": w}))
		else:
			
			# undirected case
			# assert False 
			for _, v, w in g.edges(n, data="weight"):
				if v in subgraph:
					v = subgraph
				new_edges.append(
					(subgraph, v, {"weight": w}))

		g.add_edges_from(new_edges)
		g.remove_node(n)

def bottom_up_partition(g, 
	score_function=score_subgraph_density,
	subgraph_sizes=(2, 3)):

	'''
	perform partition in bottom-up manner
	'''

	g = nx.MultiDiGraph(g.copy())
	# g = nx.DiGraph(g.copy())
	# g = nx.MultiGraph(g.copy())
	num_edges = len(g.edges())
	h = nx.DiGraph()

	h.add_nodes_from( g.nodes() )

	# ## handle self loops
	for u, _ in nx.selfloop_edges(g):
		h.add_edge(u, frozenset([u]))
	g = nx.relabel_nodes(g, 
		mapping={n: frozenset([n]) 
		for n, _ in nx.selfloop_edges(g)})

	max_size = 15

	# repeat until g has collapsed into a single node    
	while len(g) > 1:
		print ("number of nodes in g:", len(g),
			"number of edges:", len(g.edges()))

		s = 2

		while s < max_size :
			print ("scoring subgraphs of size", s)
			subgraph_scores = {}

			# for cycle in map(frozenset, nx.simple_cycles(g)):

			# for s in subgraph_sizes:
			for n in g.nodes():
				for cycle in map(frozenset, 
					find_cycles(n, s, g, start=n)):
					# find_chains(n, s, g, )):
					if cycle in subgraph_scores:
						assert np.allclose(subgraph_scores[cycle], 
							score_function(g, cycle))
						# assert False
						continue
					subgraph_scores[cycle] = \
						score_function(g, cycle)


			if len(subgraph_scores) > 0:
				chosen_subgraph = max(subgraph_scores, 
					key=subgraph_scores.get)
				# chosen_subgraph = sorted_subgraphs.pop(0)
				chosen_subgraph_score = subgraph_scores[chosen_subgraph]

				if chosen_subgraph_score > 0:
					print ("found positive scoring subgraph of size", 
						len(chosen_subgraph))
					break
			s += 1

		if s == max_size:
			print ("failed to find any postive scorign subgraphs of size", max_size)
			chosen_subgraph = frozenset().union(g.nodes())

		collapse_subgraph(g, chosen_subgraph)

		# add chosen subgraph to h
		h.add_node(chosen_subgraph)
		for n in chosen_subgraph:
			h.add_edge(n, chosen_subgraph)

		assert len(g.edges()) == num_edges
		


		# if len(subgraph_scores) > 0:

		# 	# determine all highest scoring subgraphs
		# 	# sorted_subgraphs = sorted(subgraph_scores, 
		# 		# key=lambda x: subgraph_scores[x],
		# 		# reverse=True)
		# 	chosen_subgraph = max(subgraph_scores, 
		# 		key=subgraph_scores.get)
		# 	# chosen_subgraph = sorted_subgraphs.pop(0)
		# 	chosen_subgraph_score = subgraph_scores[chosen_subgraph]

		# 	if chosen_subgraph_score > 0:

		# 		chosen_subgraphs = [chosen_subgraph]

		# 		# for subgraph in sorted_subgraphs:
		# 		# 	if subgraph_scores[subgraph] < chosen_subgraph_score:
		# 		# 		break
		# 		# 	chosen_subgraphs.append(subgraph)

		# 		# combine any chosen subgraphs that contain the 
		# 		# same nodes
		# 		# if len(chosen_subgraphs) > 1:

		# 		# 	print ("number of chosen subgraphs before collapse", 
		# 		# 		len(chosen_subgraphs))

		# 		# 	overlaps = np.array([[len(x.intersection(y)) 
		# 		# 		for y in chosen_subgraphs]
		# 		# 		for x in chosen_subgraphs])
		# 		# 	overlap_g = nx.Graph(overlaps)
		# 		# 	chosen_subgraphs = [frozenset().union([x 
		# 		# 		for c in cc for x in chosen_subgraphs[c]]) 
		# 		# 		for cc in nx.connected_components(overlap_g)]

		# 		# 	print ("number of chosen subgraphs after collapse", 
		# 		# 		len(chosen_subgraphs))
		# 	else:
		# 		# could not find a subgraph with positive score
		# 		print ("no subgraphs with positive score", )
		# 		chosen_subgraphs = [frozenset().union([x for x in g])]

		# else:
		# 	# could not find a subgraph of selected sizes
		# 	print ("no subgraphs of sizes", subgraph_sizes)
		# 	chosen_subgraphs = [frozenset().union([x for x in g])]

		# for chosen_subgraph in chosen_subgraphs:

		# 	# remove scores associated with merged nodes
		# 	# print ("removing scores associated with", 
		# 	# 	"chosen subgraph")
		# 	# subgraph_scores = {k: v 
		# 	# 	for k, v in subgraph_scores.items()
		# 	# 	if not any([x in k for x in chosen_subgraph])}
		# 	# print ("done")
			
		# 	# merge subgraph into super-node
			
		# 	collapse_subgraph(g, subgraph_scores)

		# 	assert len(g.edges()) == num_edges

		# 	# add chosen subgraph to h
		# 	h.add_node(chosen_subgraph)
		# 	for n in chosen_subgraph:
		# 		h.add_edge(n, chosen_subgraph)

		# 	# add cycles containing new node
		# 	# print ("determining new cycles containing new node")

		# 	# for s in subgraph_sizes:
		# 	# 	for cycle in map(frozenset, 
		# 	# 	find_cycles(chosen_subgraph, s, g, start=chosen_subgraph)):
		# 	# 	# find_chains(chosen_subgraph, s, g,)):
		# 	# 		# assert nx.is_strongly_connected(g.subgraph(cycle))
		# 	# 		if cycle in subgraph_scores:
		# 	# 			assert np.allclose (subgraph_scores[cycle], score_function(g, cycle))
		# 	# 			# assert False
		# 	# 			continue
		# 	# 		subgraph_scores[cycle] = \
		# 	# 			score_function(g, cycle)
		# 	# print ("done")

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
		print ("processing SCC", scc)
		scc = g.subgraph(scc)
		scc_tree = bottom_up_partition(scc, 
			score_function=score_function,
			subgraph_sizes=subgraph_sizes)
		degrees = dict(scc_tree.out_degree())
		root = min(degrees, key=degrees.get)
		roots.append(root)
		h = nx.union(h, scc_tree)
		print ()

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
			"num_loops", "bc"],
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
	print ("using score function:", score_function, )

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
	elif score_function == "num_loops":
		raise Exception
		score_function = score_subgraph_num_loops
	elif score_function == "bc":
		score_function = score_subgraph_bc
	else:
		raise Exception

	h = decompose_all_sccs(g, 
		score_function=score_function,
		subgraph_sizes=range(2, 4))
	
	print ("determining merge depths")

	out_degrees = dict(h.out_degree())
	root = min(out_degrees, key=out_degrees.get)

	core = list(
		max(nx.strongly_connected_components(g),
		key=len))

	merge_depths = {node: \
		nx.shortest_path_length(h, node, root) 
		for node in core }

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

	score_filename = os.path.join(output_dir, 
		"in_component_scores.csv")
	print ("saving scores to", score_filename)
	pd.DataFrame.from_dict(scores, 
		orient="index").to_csv(score_filename)
	print ()


	### DRAWING
	if args.draw:

		draw_height = 2

		print ("DRAWING DECOMPOSITION AT HEIGHT",
			draw_height)

		h = h.reverse()

		node_id_map = {}
		node_height_map = {}

		for i, n in enumerate(h.nodes()):
			if isinstance(n, frozenset):
				g_ = nx.MultiDiGraph(g.subgraph(unpack(n)))
			else:
				g_ = nx.MultiDiGraph(g.subgraph([n]))
				g_.remove_edges_from(list(nx.selfloop_edges(g_)))
			nx.set_edge_attributes(g_, name="arrowhead",
				values={(u, v, w): ("normal" if w>0 else "tee") 
					for u, v, w in g_.edges(data="weight")})

			node_id_map.update({n: i})

			children = list(h.neighbors(n))
			if len(children) == 0:
				height = 0
			else:
				height = max([nx.shortest_path_length(h, n, 
					el) for el in unpack(n)])

			if height > draw_height:
				continue
			
			node_height_map.update({n: height})

			for no, child in enumerate(children):
				# make metanode
				node = "metanode_{}".format(no)
				image_filename = os.path.join(output_dir, 
					"subgraph_{}.png".format(node_id_map[child]))
				assert os.path.exists(image_filename)
				g_.add_node(node, label="", 
					image=image_filename)
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

			# print ("plotted", plot_filename)
		
		h = h.reverse()

		# skip drawing tree for now

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

		## cleanup of directory
		print ("cleaning up directory")
		nodes_to_remove = [node_id_map[n]
			for n, h in node_height_map.items()
			if h < draw_height]
		# for f in glob.iglob(os.path.join(output_dir, 
		# 	"subgraph_*.png")):
		for f in (os.path.join(output_dir, 
			"subgraph_{}.png".format(n) )
			for n in nodes_to_remove):
			os.remove(f)

if __name__ == "__main__":
	main()