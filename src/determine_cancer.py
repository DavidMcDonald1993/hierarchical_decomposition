import os

import numpy as np
import pandas as pd 


import itertools

import PyBoolNet
from PyBoolNet import StateTransitionGraphs as STGs

from egfr import build_STG_and_determine_attractors, compute_average_activation

def main():

	update = "synchronous"

	primes = PyBoolNet.FileExchange.read_primes("datasets/EGFR_full/egfr_primes.json")


	proliferation_grow_tfs = set(["elk1", "creb", "ap1", "cmyc", 
		"p70s6_2", "hsp27"])
	apoptosis_tfs = set(["pro_apoptotic"])

	# input_genes = set(PyBoolNet.PrimeImplicants.find_inputs(primes))
	potential_targets = {"erbb1", "erbb2", "erbb3", "erbb4", "erbb11", "erbb12", "erbb13", "erbb14", "erbb23", "erbb24", "erbb34", "erbb44", "shp1"}

	output_genes = proliferation_grow_tfs.union(apoptosis_tfs)

	for gene in output_genes.union(potential_targets):
		assert gene in primes, gene

	num_state_samples = 10000

	print ("state space is too large -- sampling", 
			num_state_samples, 
			"states")
	states = set()
	while len(states) < num_state_samples:
		state = tuple(np.random.randint(2, size=len(primes)))
		states.add(state)
	states = list(map(lambda state: 
		STGs.state2str({p: s for p, s in zip(primes, state)}), states))

	print ("completed sampling states -- determining attractors")

	attractors, _ = build_STG_and_determine_attractors(primes, states)
	
	print ("ORIGINAL NETWORK")

	print ("found", len(set(map(frozenset, attractors))), 
		"unique attractors")

	gene_counts_original = compute_average_activation(primes, 
		genes=output_genes,
		attractors=attractors)

	# make a dataframe for each output gene
	output_dfs = {gene: pd.DataFrame() for gene in output_genes}

	## add original network
	for output_gene in output_genes:
		output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_original[output_gene], name="original"))
		output_dfs[output_gene].to_csv(os.path.join("results", "{}_mutations.csv".format(output_gene)))


	print ("completed original network")

	for n_genes in range(1, 4):

		for mutated_genes in itertools.combinations(potential_targets, n_genes):

			print ("turning on", mutated_genes)

			modified_network = PyBoolNet.PrimeImplicants.\
				create_constants(primes, 
				{mutated_gene: 1 for mutated_gene in mutated_genes}, 
				Copy=True)
			print ("constants are", PyBoolNet.\
					PrimeImplicants.find_constants(modified_network))

			attractors, _ = build_STG_and_determine_attractors(modified_network, states)

			gene_counts_modified = compute_average_activation(modified_network, 
				genes=output_genes,
				attractors=attractors)

			for output_gene in output_genes:
				output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(mutated_genes)))
				output_dfs[output_gene].to_csv(os.path.join("results", "{}_mutations.csv".format(output_gene)))


	for output_gene in output_genes:
		output_dfs[output_gene].to_csv(os.path.join("results", "{}_mutations.csv".format(output_gene)))

if __name__ == "__main__":
	main()