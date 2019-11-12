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

    output_genes = args.output_genes

    edgelist_filename = args.edgelist
    print ("loading interaction graph from", edgelist_filename)
    g = nx.read_weighted_edgelist(edgelist_filename, 
        delimiter="\t",
        create_using=nx.DiGraph())
    
    for gene in output_genes:
        assert gene in g, gene

    core = max(nx.strongly_connected_components(g), 
        key=len)

    # remove any canger genes from consideration
    cancer_mutuation = args.cancer_mutation
    if cancer_mutuation:
        core = core - set(cancer_mutuation)
    core -= set(output_genes)
    core = sorted(core)

    print ("core is", core)

    possible_candidates = [("cancer",)] + [genes for n_genes in [1, 2] for genes in itertools.combinations(core, n_genes)]

    print ("number of possible candidates is", len(possible_candidates))

    i = args.i
    chosen_candidates = possible_candidates[i*chunksize : 
        (i+1)*chunksize]



    output_filenames = {output_gene: 
        os.path.join(output_dir,
        "{}_expressions_chunk_{}.csv".format(output_gene, i))
        for output_gene in output_genes}

    output_dfs = {gene: pd.DataFrame() 
        for gene in output_genes}


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
        attractors = build_STG_and_determine_attractors(modified_network, 
            states)
        
        print ("determing activations for output genes")
        gene_counts = compute_average_activation(modified_network, 
            genes=output_genes,
            attractors=attractors)

        for output_gene in output_genes:
            output_dfs[output_gene] = output_dfs[output_gene]\
                .append(pd.Series(gene_counts[output_gene], name=chosen_candidate_identifier))

    print ("writing results to file")
    for output_gene in output_genes:
        output_dfs[output_gene].to_csv(output_filenames[output_gene])

if __name__ == "__main__":
    main()