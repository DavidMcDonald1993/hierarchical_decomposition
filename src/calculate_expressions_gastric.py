import os

import numpy as np 
import networkx as nx
import pandas as pd

import itertools

import PyBoolNet
from PyBoolNet import Attractors
from PyBoolNet import StateTransitionGraphs as STGs

import matplotlib.pyplot as plt

def select_states(primes, num_state_samples=10000, seed=0):

    np.random.seed(seed)

    print ("sampling", 
            num_state_samples, 
            "states")
    states = set()
    while len(states) < num_state_samples:
        state = tuple(np.random.randint(2, size=len(primes)))
        states.add(state)
    states = list(map(lambda state: 
        STGs.state2str({p: s 
        for p, s in zip(primes, state)}), states))

    return states

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


# def plot_results_df(df, filename):

#     num_genes = df.shape[1]

#     ind = np.arange(num_genes)  * 10# the x locations for the groups
#     width = 1  # the width of the bars

#     fig, ax = plt.subplots(figsize=[10,10])

#     for i, col in enumerate(df.columns):
#         data = df[col]
        
#         x = ind + (-num_genes//2 + i) * width

#         print (x, )
#         print (data)
#         print (width)
#         ax.bar(x, data, width, label=col)
        
#     # Add some text for labels, title and custom x-axis tick labels, etc.
#     ax.set_ylabel('Difference in activation from original network')
#     ax.set_title('Test1')
#     ax.set_xticks(ind)
#     ax.set_xticklabels(df.index)
#     ax.legend()

#     plt.savefig(filename)

#     plt.show()

def main():

    output_dir = os.path.join("results", "gastric")
    if not os.path.exists(output_dir):
        print ("making", output_dir)
        os.makedirs(output_dir, exist_ok=True)

    update = "synchronous"

    primes = PyBoolNet.FileExchange.read_primes("datasets/gastric/gastric_primes.json")

  
    # core_of_the_core = set(["pi3k", "pip3", "gab1"])

    # proliferation_grow_tfs = ["elk1", "creb", "ap1", "cmyc", 
    #     "p70s6_2", "hsp27"]
    # apoptosis_tfs = ["pro_apoptotic"]
    anti_survival = ["Caspase8", "Caspase9", "FOXO"]
    pro_survival = ["RSK", "TCF", "cMYC"]

    output_genes = pro_survival + anti_survival


    g = nx.read_weighted_edgelist("datasets/gastric/edgelist.tsv", 
        delimiter="\t", create_using=nx.DiGraph())

    # remove output nodes from consideration
    core = sorted(max(nx.strongly_connected_components(g), 
        key=len) - set(output_genes))

    print ("core is", core)
    print ("core contains", len(core), "nodes")


    for gene in output_genes:
        assert not gene in core, gene
        assert gene in primes, gene

    output_filenames = {output_gene: 
        os.path.join(output_dir,
        "{}_expressions.csv".format(output_gene))
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

        attr = PyBoolNet.Attractors.find_attractor_state_by_randomwalk_and_ctl(primes, update)

        gene_counts_original = compute_average_activation(primes, 
            genes=output_genes,
            attractors=attractors)

        for output_gene in output_genes:

            output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_original[output_gene], name="original"))
            output_dfs[output_gene].to_csv(output_filenames[output_gene])    

    for n_genes in [1, 2]:
        for core_genes in itertools.combinations(core, n_genes):
    # for core_genes in [["MEK", "ERK"]]:

            if "_".join(core_genes) not in output_dfs[output_genes[0]].index:

                print ("turning off", core_genes)

                modified_network = PyBoolNet.PrimeImplicants.\
                    create_constants(primes, 
                    {core_gene: 0 for core_gene in core_genes}, 
                    Copy=True)

                attractors, _ = build_STG_and_determine_attractors(modified_network, states)

                gene_counts_modified = compute_average_activation(modified_network, 
                    genes=output_genes,
                    attractors=attractors)

                for output_gene in output_genes:
                    output_dfs[output_gene] = output_dfs[output_gene].append(pd.Series(gene_counts_modified[output_gene], name="_".join(core_genes)))
                    output_dfs[output_gene].to_csv(output_filenames[output_gene])


if __name__ == "__main__":
    main()