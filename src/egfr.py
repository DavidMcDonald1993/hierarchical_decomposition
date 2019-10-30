import os

import numpy as np 
import networkx as nx
import pandas as pd

import itertools

import PyBoolNet
from PyBoolNet import Attractors
from PyBoolNet import StateTransitionGraphs as STGs

import matplotlib.pyplot as plt

def build_STG_and_determine_attractors(primes, states):

    ## assume synchronous

    assert isinstance(states[0], str)

    stg = nx.DiGraph()

    attractors = []

    for i, state in enumerate(states):

        next_state = STGs.state2str(STGs.successor_synchronous(primes, state))

        while next_state not in stg:

            stg.add_edge(state, next_state)
            state = next_state
            next_state = STGs.state2str(STGs.successor_synchronous(primes, state))

            # if len(stg) % 1000 == 0:
            #     print ("graph contains", len(stg), "states")

        assert next_state in stg
        stg.add_edge(state, next_state)


        if state == next_state:
            # print ("found steady state attractor", next_state)
            attractors.append([state])
                
        else:

            visited = [state]
            while next_state not in visited:
                visited.append(next_state)
                assert len(list(stg.neighbors(next_state))) == 1

                next_state = list(stg.neighbors(next_state))[0]

            idx = visited.index(next_state)
            attractor = visited[idx:]
            # print ("found cyclic attractor with period", len(attractor))
            attractors.append(attractor)

        # if i % 1000 == 0:
        #     print ("processed state {:04d}/{}".format(i, len(states)))


    return attractors, stg

def compute_average_activation(primes, genes, attractors):

    counts = {gene: [] for gene in genes}

    for attractor in attractors:

        attractor_counts = {gene: 0 for gene in genes}

        for state in attractor:

            state_dict = STGs.state2dict(primes, state)

            for gene in genes:
                attractor_counts[gene] += state_dict[gene]

        # attractor_counts = {gene: count/len(attractor) for gene, count in attractor_counts.items()}
        for gene in genes:
            counts[gene].append(attractor_counts[gene] / len(attractor))

    # return {gene: count/len(attractors) for gene, count in counts.items()}
    return counts


def plot_results_df(df, filename):

    num_genes = df.shape[1]

    ind = np.arange(num_genes)  * 10# the x locations for the groups
    width = 1  # the width of the bars

    fig, ax = plt.subplots(figsize=[10,10])

    for i, col in enumerate(df.columns):
        data = df[col]
        
        x = ind + (-num_genes//2 + i) * width

        print (x, )
        print (data)
        print (width)
        ax.bar(x, data, width, label=col)
        
    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Difference in activation from original network')
    ax.set_title('Test1')
    ax.set_xticks(ind)
    ax.set_xticklabels(df.index)
    ax.legend()

    plt.savefig(filename)

    plt.show()

def main():

    update = "synchronous"

    # primes = PyBoolNet.FileExchange.bnet2primes("datasets/EGFR_full/egfr.bnet", 
        # FnamePRIMES="datasets/EGFR_full/egfr_primes.json")
    primes = PyBoolNet.FileExchange.read_primes("datasets/EGFR_full/egfr_primes.json")


    core_of_the_core = set(["pi3k", "pip3", "gab1"])

    proliferation_grow_tfs = set(["elk1", "creb", "ap1", "cmyc", 
        "p70s6_2", "hsp27"])
    apoptosis_tfs = set(["pro_apoptotic"])
    # additional_genes_of_interest = set(["akt", "ras"])

    input_genes = set(PyBoolNet.PrimeImplicants.find_inputs(primes))
    output_genes = proliferation_grow_tfs.union(apoptosis_tfs)

    for gene in output_genes:
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

    ## add original
    for output_gene in output_genes:
        output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_original[output_gene], name="original"))



  


    cancer_network =  PyBoolNet.PrimeImplicants.\
        create_constants(primes, 
        {"erbb1": 1,}, 
        Copy=True)


    attractors, _ = build_STG_and_determine_attractors(cancer_network, states)
    
    print ("CANCER NETWORK")

    print ("found", len(set(map(frozenset, attractors))), 
        "unique attractors")

    gene_counts_cancer = compute_average_activation(cancer_network, 
        genes=output_genes,
        attractors=attractors)

    # for k, v in gene_counts_original.items():
    #     # print ("{}\t{}".format(k, v))
    #     plt.hist(v)
    #     plt.title(k)
    #     plt.show()

    # print ()
    # raise SystemExit

    # ## build dataframe
    # raw_value_df = pd.DataFrame()
    # # add original network
    # raw_value_df = raw_value_df.append(pd.Series(gene_counts_original, name="original"))

    # difference_df = pd.DataFrame()

    # make a dataframe for each output gene
    # output_dfs = {gene: pd.DataFrame() for gene in output_genes}

    ## add cancer
    for output_gene in output_genes:
        output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_cancer[output_gene], name="cancer"))


    for n_genes in range(1, len(core_of_the_core)+1):

        for core_genes in itertools.combinations(core_of_the_core, n_genes):

            print ("turning off", core_genes)

            modified_network = PyBoolNet.PrimeImplicants.\
                create_constants(cancer_network, 
                {core_gene: 0 for core_gene in core_genes}, 
                Copy=True)
            print ("constants are", PyBoolNet.\
                    PrimeImplicants.find_constants(modified_network))

            attractors, _ = build_STG_and_determine_attractors(modified_network, states)

            # print ("found", len(set(map(frozenset, attractors))), 
            #     "unique attractors")

            gene_counts_modified = compute_average_activation(modified_network, 
                genes=output_genes,
                attractors=attractors)

            # for k, v in gene_counts_modified.items():
            #     print ("{}\t{}".format(k, v))

            # print ()

            # raw_value_df = raw_value_df.append(pd.Series(gene_counts_modified, name="_".join(core_genes)))

            # gene_counts_difference = {gene: gene_counts_modified[gene] - gene_counts_original[gene] for gene in gene_counts_modified}

            # difference_df = difference_df.append(pd.Series(gene_counts_difference, name="_".join(core_genes)))
            for output_gene in output_genes:
                output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(core_genes)))

    # raw_value_df.to_csv("raw_values.csv")
    # difference_df.to_csv("differences.csv")

    # print ("plotting raw values")
    # plot_results_df(raw_value_df, "raw_values.png")

    # print ("plotting differences")
    # plot_results_df(difference_df, "differences.png")

    potential_targets = set(primes) - input_genes - output_genes - core_of_the_core - {"erbb1", }


    
    for n_genes in [1]:

        for non_core_genes in itertools.combinations(potential_targets, n_genes):

            print ("turning off", non_core_genes)

            modified_network = PyBoolNet.PrimeImplicants.\
                create_constants(cancer_network, 
                {non_core_gene: 0 for non_core_gene in non_core_genes}, 
                Copy=True)
            print ("constants are", PyBoolNet.\
                    PrimeImplicants.find_constants(modified_network))

            attractors, _ = build_STG_and_determine_attractors(modified_network, states)

            # print ("found", len(set(map(frozenset, attractors))), 
            #     "unique attractors")

            gene_counts_modified = compute_average_activation(modified_network, 
                genes=output_genes,
                attractors=attractors)

            # for k, v in gene_counts_modified.items():
            #     print ("{}\t{}".format(k, v))

            # print ()

            # raw_value_df = raw_value_df.append(pd.Series(gene_counts_modified, name="_".join(core_genes)))

            # gene_counts_difference = {gene: gene_counts_modified[gene] - gene_counts_original[gene] for gene in gene_counts_modified}

            # difference_df = difference_df.append(pd.Series(gene_counts_difference, name="_".join(core_genes)))
            for output_gene in output_genes:
                output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(non_core_genes)))

    for output_gene in output_genes:
        output_dfs[output_gene].to_csv(os.path.join("results", "{}_expressions.csv".format(output_gene)))
    


if __name__ == "__main__":
    main()