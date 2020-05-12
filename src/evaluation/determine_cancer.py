import os

import numpy as np
import pandas as pd 


import itertools

import PyBoolNet
from PyBoolNet import StateTransitionGraphs as STGs

from calculate_expressions import select_states, build_STG_and_determine_attractors, compute_average_activation

def main():

	output_dir = os.path.join("results", "cancer")
	if not os.path.exists(output_dir):
		print ("making", output_dir)
		os.makedirs(output_dir, exist_ok=True)

	update = "synchronous"

	primes = PyBoolNet.FileExchange.read_primes("datasets/EGFR_full/egfr_primes.json")

	proliferation_grow_tfs = ["elk1", "creb", "ap1", "cmyc", 
		"p70s6_2", "hsp27"]
	apoptosis_tfs = ["pro_apoptotic"]


	output_genes = proliferation_grow_tfs + apoptosis_tfs

	potential_targets = ["erbb1", "erbb2", "erbb3", "erbb4", 
		"erbb11", "erbb12", "erbb13", "erbb14", "erbb23", 
		"erbb24", "erbb34", "erbb44", "shp1"]

	for gene in output_genes + potential_targets:
		assert gene in primes, gene


	output_filenames = {output_gene: 
		os.path.join(output_dir,
		"{}_mutations.csv".format(output_gene))
		for output_gene in output_genes}

	states = select_states(primes)

	# make a dataframe for each output gene
	if os.path.exists(output_filenames[output_genes[0]]):
		output_dfs = {gene: pd.read_csv(output_filenames[gene], index_col=0) 
			for gene in output_genes}
		for output_gene in output_genes:
			output_dfs[output_gene].columns = \
				[int(c) for c in output_dfs[output_gene].columns]
	else:
		output_dfs = {gene: pd.DataFrame() for gene in output_genes}

	## add original
	if "original" not in output_dfs[output_genes[0]].index:

		print ("determining attractors for original network")

		attractors, _ = build_STG_and_determine_attractors(primes, states)
		
		gene_counts_original = compute_average_activation(primes, 
			genes=output_genes,
			attractors=attractors)

		for output_gene in output_genes:

			output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_original[output_gene], name="original"))
			output_dfs[output_gene].to_csv(output_filenames[output_gene])

	for n_genes in range(1, 4):

		for mutated_genes in itertools.combinations(potential_targets, n_genes):

	# for mutated_genes in [
	# 	["erbb11"],
	# 	["erbb1",  "erbb11", "erbb12", "erbb13", "erbb14"],
	# 	["erbb2", "erbb23", "erbb24"]
	# ]:

			if "_".join(mutated_genes) not in output_dfs[output_genes[0]].index:

				print ("turning on", mutated_genes)

				modified_network = PyBoolNet.PrimeImplicants.\
					create_constants(primes, 
					{mutated_gene: 1 for mutated_gene in mutated_genes}, 
					Copy=True)

				attractors, _ = build_STG_and_determine_attractors(modified_network, states)

				gene_counts_modified = compute_average_activation(modified_network, 
					genes=output_genes,
					attractors=attractors)

				for output_gene in output_genes:
					output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(mutated_genes)))
					output_dfs[output_gene].to_csv(output_filenames[output_gene])


if __name__ == "__main__":
	main()