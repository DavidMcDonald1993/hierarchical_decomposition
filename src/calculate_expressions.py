import os
from pathlib import Path

import numpy as np 
import networkx as nx
import pandas as pd

import itertools

import PyBoolNet
from PyBoolNet import Attractors
from PyBoolNet import StateTransitionGraphs as STGs

from utils import select_states, build_STG_and_determine_attractors
from utils import compute_average_activation, threadsafe_save_test_results

import argparse

def touch(filename):
    Path(filename).touch()

def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(description="Compute activations for boolean networks")

    parser.add_argument("-i", 
        dest="i", type=int, default=0,
        help="index of list of modification.")

    parser.add_argument("--edgelist", 
        dest="edgelist", type=str, default=None,
        help="edgelist to load.")

    parser.add_argument("--primes", 
        dest="primes", type=str, default=None,
        help="primes to load.")
    parser.add_argument("--cancer_mutation", 
        dest="cancer_mutation", type=str, default=None, nargs="*",
        help="genes to set to overexpressed.")
    

    parser.add_argument("--output", dest="output", 
        type=str, default=None,
        help="Directory to save results.")

    parser.add_argument("--output_genes", dest="output_genes", 
        type=str, default=None, nargs="+",
        help="Directory to save results.")

    return parser.parse_args()

def main():

    chunksize = 10

    args = parse_args()

    output_dir = args.output
    if not os.path.exists(output_dir):
        print ("making", output_dir)
        os.makedirs(output_dir, exist_ok=True)


    edgelist_filename = args.edgelist
    print ("loading interaction graph from", edgelist_filename)
    g = nx.read_weighted_edgelist(edgelist_filename, 
        delimiter="\t", create_using=nx.DiGraph())

    core = max(nx.strongly_connected_components(g), 
        key=len)

    # remove any canger genes from consideration
    cancer_mutuation = args.cancer_mutation
    if cancer_mutuation:
        core = core - set(cancer_mutuation)
    core = sorted(core)

    print ("core is", core)

    possible_candidates = [("cancer",)] + [genes for n_genes in [1, 2] for genes in itertools.combinations(core, n_genes)]

    print ("number of possible candidates is", len(possible_candidates))

    i = args.i
    chosen_candidates = possible_candidates[i*chunksize : 
        (i+1)*chunksize]
    # chosen_candidate_identifier = "_".join(chosen_candidate)
    # print ("chosen candidate is", chosen_candidate_identifier)

    # proliferation_grow_tfs = ["elk1", "creb", "ap1", "cmyc", 
    #     "p70s6_2", "hsp27"]
    # apoptosis_tfs = ["pro_apoptotic"]

    # output_genes = proliferation_grow_tfs + apoptosis_tfs
    output_genes = args.output_genes

    for gene in output_genes:
        assert gene in g, gene

    output_filenames = {output_gene: 
        os.path.join(output_dir,
        "{}_expressions_chunk_{}.csv".format(output_gene, i))
        for output_gene in output_genes}
    # lock_filenames = {output_gene: 
    #     os.path.join(output_dir,
    #     "{}_expressions_chunk_{}.lock".format(output_gene, i))
    #     for output_gene in output_genes}
    # for lock_filename in lock_filenames.values():
    #     touch(lock_filename)
        # assert os.path.exists(lock_filename), lock_filename
    # make a dataframe for each output gene

    # if False:#os.path.exists(output_filenames[output_genes[0]]):
    #     print ("loading dfs")
    #     output_dfs = {gene: pd.read_csv(output_filenames[gene], index_col=0) 
    #         for gene in output_genes}
    #     print ("loaded dfs")
    #     for output_gene in output_genes:
    #         output_dfs[output_gene].columns = \
    #             [int(c) for c in output_dfs[output_gene].columns]
    #         if chosen_candidate_identifier in output_dfs[output_gene].index:
    #             print (chosen_candidate_identifier, "is already done, terminating")
    #             return
    #     print ("mapped columns")
    # else:
    output_dfs = {gene: pd.DataFrame() for gene in output_genes}


    primes_filename = args.primes
    if primes_filename.endswith(".bnet"):
        print ("loading from bnet file", primes_filename)
        json_filename = primes_filename.replace(".bnet", ".json")
        print("saving primes json to", json_filename)
        primes = PyBoolNet.FileExchange.bnet2primes(primes_filename, 
            FnamePRIMES=json_filename)
    else:
        assert primes_filename.endswith(".json")
        print ("loading primes from json", primes_filename)
        primes = PyBoolNet.FileExchange.read_primes(primes_filename)

    if cancer_mutuation:
        print ("turning on", "_".join(cancer_mutuation))
        primes = PyBoolNet.PrimeImplicants.\
            create_constants(primes, 
            {mutation: 1 for mutation in cancer_mutuation}, 
            Copy=True)

    for gene in output_genes:
        assert gene in primes, gene

    states = select_states(primes)

    update = "synchronous"

    for chosen_candidate in chosen_candidates:
        chosen_candidate_identifier = "_".join(chosen_candidate)
        print ("chosen candidate is", chosen_candidate_identifier)


        # mad modification to network if necessary
        if chosen_candidate_identifier is not "cancer":
            
            print ("switching off", chosen_candidate)
            
            modified_network = PyBoolNet.PrimeImplicants.\
                create_constants(primes, 
                {gene: 0 for gene in chosen_candidate}, 
                Copy=True)

        else:
            # do nothing for cancer
            modified_network = primes

        print ("determining attractors")
        attractors = build_STG_and_determine_attractors(primes, states)
        
        print ("determing activations for output genes")
        gene_counts = compute_average_activation(primes, 
            genes=output_genes,
            attractors=attractors)

        # print ("done")

        # print ("writing results to file")

        for output_gene in output_genes:
            # threadsafe_save_test_results(lock_filenames[output_gene], 
            # output_filenames[output_gene],
            # chosen_candidate_identifier,
            # gene_counts[output_gene]
            # )
            output_dfs[output_gene] = output_dfs[output_gene]\
                .append(pd.Series(gene_counts[output_gene], name=chosen_candidate_identifier))



    print ("writing results to file")
    for output_gene in output_genes:
        output_dfs[output_gene].to_csv(output_filenames[output_gene])

    # raise SystemExit

    # ## add original
    # if "original" not in output_dfs[output_genes[0]].index:

    #     print ("determining attractors for original network")

    #     attractors = build_STG_and_determine_attractors(primes, states)
        
    #     gene_counts_original = compute_average_activation(primes, 
    #         genes=output_genes,
    #         attractors=attractors)

    #     for output_gene in output_genes:

    #         output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_original[output_gene], name="original"))
    #         output_dfs[output_gene].to_csv(output_filenames[output_gene])

    # # make cancer network
    # cancer_network =  PyBoolNet.PrimeImplicants.\
    #     create_constants(primes, 
    #     {mutation: 1 for mutation in mutations}, 
    #     Copy=True)
    
    # ## add cancer
    # if "cancer" not in output_dfs[output_genes[0]].index:

    #     print ("determining attractors for cancer network")

    #     attractors  = build_STG_and_determine_attractors(cancer_network, states)
        
    #     gene_counts_cancer = compute_average_activation(cancer_network, 
    #         genes=output_genes,
    #         attractors=attractors)

    #     for output_gene in output_genes:

    #         output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_cancer[output_gene], name="cancer"))
    #         output_dfs[output_gene].to_csv(output_filenames[output_gene])

    # for n_genes in [1, 2]:

    #     print ("n_genes", n_genes)

    #     for core_genes in itertools.combinations(core, n_genes):

    #         if "_".join(core_genes) not in output_dfs[output_genes[0]].index:

    #             print ("turning off", core_genes)

    #             modified_network = PyBoolNet.PrimeImplicants.\
    #                 create_constants(cancer_network, 
    #                 {core_gene: 0 for core_gene in core_genes}, 
    #                 Copy=True)

    #             attractors  = build_STG_and_determine_attractors(modified_network, states)

    #             gene_counts_modified = compute_average_activation(modified_network, 
    #                 genes=output_genes,
    #                 attractors=attractors)

    #             for output_gene in output_genes:
    #                 output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(core_genes)))
    #                 output_dfs[output_gene].to_csv(output_filenames[output_gene])


if __name__ == "__main__":
    main()