import numpy as np
import networkx as nx
from networkx.drawing.nx_agraph import write_dot, graphviz_layout, to_agraph

import pandas as pd 

import PyBoolNet
from PyBoolNet import StateTransitionGraphs as STGs
from PyBoolNet import AspSolver, Attractors

import itertools

import matplotlib.pyplot as plt

from scipy.stats import norm, normaltest

import random
random.seed(0)

import argparse

import os

from .egfr import build_STG_and_determine_attractors

def determine_attractors(primes, 
    update="synchronous", 
    num_state_samples=10000,
    states=None,
    verbose=False):

    if len(primes) < 16:

        if verbose:
            print ("state space is small enough -- computing all attractors")
        
            print ("computing stage transition graph")
        stg = STGs.primes2stg(primes, update)

        if verbose:
            print ("determining state trajectories")
        
        attrs = {}

        for basin in nx.weakly_connected_component_subgraphs(stg):
            assert len(list(nx.simple_cycles(basin))) == 1
            cycle = next(nx.simple_cycles(basin)) 
            for n in basin:
                attrs.update({n: cycle})

    else:

        if verbose:
            print ("state space is too large --  sampling", 
                num_state_samples, 
                "states")
        if states is None: 
            states = set()
            while len(states) < num_state_samples:
                state = tuple(np.random.randint(2, size=len(primes)))
                states.add(state)
            states = list(map(lambda state: 
                STGs.state2str({p: s for p, s in zip(primes, state)}), states))

        if verbose:
            print ("completed sampling states -- determining trajectories")
        # attrs = {
        #     state : Attractors.find_attractor_state_by_randomwalk_and_ctl(primes, 
        #     update, 
        #     InitialState=state)
        #     for state in states   
        # }
        attrs, _ = build_STG_and_determine_attractors(primes, states)
        attrs = {state: attr for state, attr in zip(states, attrs)}

    if verbose:
        num_attractors = len(set(map(frozenset, attrs.values())))

        print ("found", num_attractors, "attractors")

    return attrs

def determinine_trajectories(
    primes, 
    states,
    time_steps=25,
    update="synchronous", 
    verbose=False
    ):

    def traj2array(d):
        return np.array([
                [   
                    [int(i) for i in st]
                    for st in d[state]
                ]
                for state in sorted(d)
                ])

    if update == "synchronous":
        update = STGs.successor_synchronous
    else:
        update = STGs.successor_asynchronous

    d = {state: [state] for state in states}

    for _ in range(time_steps):
        for state in d:
            d[state].append(STGs.state2str(update(primes, d[state][-1])))

    return traj2array(d)

def determine_derrida_curve(
    primes,
    states,
    num_repeats=1000,
    update="synchronous",
    verbose=False
    ):

    def perturb_state(state, n):
        idx = np.random.choice(len(state), n)
        return "".join((str(1-int(x)) if i in idx else x for i, x in enumerate(state)))

    n = len(primes) + 1

    if verbose:
        print ("determining derrida curve")

    if update == "synchronous":
        update = STGs.successor_synchronous
    else:
        update = STGs.successor_asynchronous

    results = {n_perturbations: 
        [0] if n_perturbations==0 
        else [] 
        for n_perturbations in range(n)}

    states = list(states)

    for n_perturbations in range(1, n):

        for s1 in random.choices(states, k=num_repeats):

            s2 = perturb_state(s1, n_perturbations)
            
            s1_ = STGs.state2str(update(primes, s1))
            s2_ = STGs.state2str(update(primes, s2))

            distance = mean_hamming_distance(
                np.array([int(x) for x in s1_]),
                np.array([int(x) for x in s2_]),
            )

            results[n_perturbations].append(distance)

    results = {k / (n-1): np.mean(v) - k / (n-1) 
        for k, v in results.items()}

    return results

def mean_hamming_distance(x, y, axis=-1):
    return np.mean(np.abs(x - y), axis=axis)

def plot_change_in_attractor(control_df, test_df):

    print ("PLOTTING CHANGE IN ATTRACTOR")

    bins = np.linspace(0, 0.26, 26 )
    n, _, _ = plt.hist(x=control_df.T , 
        bins=bins, color="b",
        alpha=0.5, rwidth=0.85)
    m, _, _ = plt.hist(x=test_df.T , 
        bins=bins, color="r",
        alpha=0.5, rwidth=0.85)

    plt.grid(axis='y', alpha=0.75)
    plt.xlabel('Value')
    plt.ylabel('Frequency')
    plt.title('Histogram of relative change in basin')
    plt.legend(["Control", "Test"])
    maxfreq = np.append(n, m).max()
    plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 
        if maxfreq % 10 else maxfreq + 10)
    plt.show()

def plot_hamming_distance_over_time(control_df, test_df):

    print ("ATTRACTOR OVER TIME")

    control_df_mean = control_df.mean(0)
    test_df_mean = test_df.mean(0)

    mean_df = pd.DataFrame([test_df_mean, control_df_mean])

    control_df_std = control_df.std(0)
    test_df_std = test_df.std(0)

    std_df = pd.DataFrame([test_df_std, control_df_std])

    # num_control = control_df.shape[0]
    # num_test = test_df.shape[0]
    num_control, num_test = 1, 1

    # df = test_df.append(control_df)

    plt.plot(mean_df.T)

    plt.xlabel('t')
    plt.ylabel('Hamming distance')
    plt.title('Plot of mean hamming distance across initial conditions')
    plt.legend(["Test-{}".format(i) 
        for i in range(num_test)] + \
        ["Control-{}".format(i)
        for i in range(num_control)])
    plt.show()

def plot_derrida(control_df, test_df):
    
    print ("DERRIDA CURVE")

    control_df_mean = control_df.mean(0)
    test_df_mean = test_df.mean(0)

    mean_df = pd.DataFrame([test_df_mean, control_df_mean])

    control_df_std = control_df.std(0)
    test_df_std = test_df.std(0)

    std_df = pd.DataFrame([test_df_std, control_df_std])

    # num_control = control_df.shape[0]
    # num_test = test_df.shape[0]
    num_control, num_test = 1, 1

    # df = test_df.append(control_df)

    plt.plot(mean_df.T)

    # plt.errorbar(data[:,0], data[:,1], yerr=data[:,2])
    # plt.plot(data[:,0], data[:,0], "k")

    plt.xlabel('H(t)')
    plt.ylabel('H(t+1)')
    plt.title('Derrida plot')
    plt.xlim([0, 1])
    plt.ylim([-1, 1])
    plt.legend(["Test", "Control"])
    plt.show()

def check_cyclic_state_equality(list_of_states_1, list_of_states_2):
    assert isinstance(list_of_states_1, list)
    assert isinstance(list_of_states_2, list)
    if len(list_of_states_1) != len(list_of_states_2):
        return False
    start = list_of_states_1[0]
    if start not in list_of_states_2:
        return False
    idx = list_of_states_2.index(start)

    return all ([x==y for x, y in zip(list_of_states_1, list_of_states_2[idx:] + list_of_states_2[:idx])])
    

def parse_args():
    '''
    Parse from command line
    '''
    parser = argparse.ArgumentParser(description="Evaluate dynamics.")

    parser.add_argument("--bnet", 
        dest="bnet", type=str, 
        help="bnet file to read dynamics from.")
    parser.add_argument("--output", 
        dest="output", type=str, 
        help="Output directory to store results.")
   
    parser.add_argument("--merge_depths", 
        dest="merge_depths", type=str, 
        help="Path to merge depth file for this network.")


    parser.add_argument("--update", 
        dest="update", type=str, default="synchronous",
        help="update method.")

    return parser.parse_args()

def main():

    args = parse_args()

    update = args.update
    num_time_steps = 10

    assert update == "synchronous"

    bnet_filename = args.bnet
    assert os.path.exists(bnet_filename)

    output_dir = args.output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)

    merge_depths_filename = args.merge_depths

    primes = PyBoolNet.FileExchange.bnet2primes(bnet_filename)

    
    possible_targets = set(primes) - {"n0", "n1"} # TODO

    # primes = PyBoolNet.PrimeImplicants.\
    #         create_constants(primes, 
    #         {"n0": 1, "n1": 0}, 
    #         Copy=True)
    assert len(primes) < 16    

    original_attractors = determine_attractors(primes, 
        update,
        verbose=True)
    num_attractors = len(set(map(frozenset, original_attractors.values())))
    assert num_attractors > 1

    original_trajectories = determinine_trajectories(primes, 
        time_steps=num_time_steps,
        states=original_attractors)

    original_derrida = determine_derrida_curve(primes, 
        states=original_attractors,
        update=update,)

    # for k, v in original_derrida.items():
    #     print (k, np.mean(v))
    # raise SystemExit
    # plot_derrida(original_derrida)

    # attrs_original = PyBoolNet.Attractors.compute_json(primes, 
    #     update, 
    #     Silent=True, )

    from collections import Counter

    steady_state_basin_sizes = dict(Counter(filter(lambda x: len(x) == 1, map(frozenset, original_attractors.values()))))
    chosen_steady_state = list(min(steady_state_basin_sizes, key=steady_state_basin_sizes.get))
    chosen_steady_state_dict = STGs.state2dict(primes, chosen_steady_state[0])
    print ("chosen steady state attractor", steady_state_basin_sizes[frozenset(chosen_steady_state)])
    for k, v in chosen_steady_state_dict.items():
        print ("{}={}".format(k, v))
    print ()

    merge_depth_df = pd.read_csv(merge_depths_filename, index_col=0, header=None )
    max_merge_depth = merge_depth_df[1].max()
    test_genes = set(merge_depth_df.index[
        merge_depth_df[1] == max_merge_depth].tolist())

    print ("test genes are", test_genes)

    for test_gene in test_genes:
        assert test_gene in primes

    ## MAKE ALL DATA FRAMES FOR RESULTS

    ## change in attractor data frames
    # change_in_attractor_df_control = pd.DataFrame()
    # change_in_attractor_df_test = pd.DataFrame()


    ## basin size data frames
    basin_size_df_control = pd.DataFrame()
    basin_size_df_test = pd.DataFrame()

    ## hamming distance over time data frames
    hamming_distance_over_time_df_control = pd.DataFrame()
    hamming_distance_over_time_df_test = pd.DataFrame()

   
    ## derrida curve data frames
    derrida_curve_df_control = pd.DataFrame()
    derrida_curve_df_test = pd.DataFrame()

    n_test_genes = len(test_genes)

    # for n_genes in [n_test_genes]:
    # for n_genes in [1]:
    for n_genes in range(1, n_test_genes + 1):

        for _genes in itertools.combinations(possible_targets, n_genes):

            # is_test_set = sum([test_gene in _genes 
            #     for test_gene in test_genes]) >= n_test_genes - 1

            is_test_set = any([test_gene in _genes 
                for test_gene in test_genes])

            # is_test_set = set(_genes) <= test_genes

            # constants = {g: 1. for g in _genes}
            constants = {g: chosen_steady_state_dict[g] 
                for g in _genes}
            
            new_primes = PyBoolNet.PrimeImplicants.\
                create_constants(primes, 
                constants, 
                Copy=True)

            all_constants = PyBoolNet.\
                PrimeImplicants.find_constants(new_primes)
            print("constants", all_constants, )

            new_attractors = determine_attractors(
                new_primes, 
                update, 
                states=original_attractors.keys()
            )
            

            # EXPERIMENT 1: change of attractors for each state
                


            # change_in_attractor = np.mean([ 
            #     not check_cyclic_state_equality(original_attractors[state],  new_attractors[state])
            #     for state in original_attractors
            #     # if all([STGs.state2dict(primes, state)[gene] == constants[gene] for gene in _genes])
            # ])

            # ## add to correct dataframe 
            # if is_test_set:
            #     # add to test dataframe
            #     change_in_attractor_df_test = change_in_attractor_df_test.append(pd.Series([change_in_attractor], name="{}".format(constants)))
            # else:
            #     # add to control dataframe
            #     change_in_attractor_df_control = change_in_attractor_df_control.append(pd.Series([change_in_attractor], name="{}".format(constants)))

            # EXPERIMENT 2: change of attractor to chosen attractor

            # change_in_attractor_to_chosen_attractor = np.mean([ 
            #     check_cyclic_state_equality(new_attractors[state], chosen_steady_state)
            #     for state in original_attractors
            #     # if all([STGs.state2dict(primes, state)[gene] == constants[gene] for gene in _genes])
            #     # if not check_cyclic_state_equality(original_attractors[state], chosen_steady_state)
            # ])

            # ## add to correct dataframe 
            # if is_test_set:
            #     # add to test dataframe
            #     basin_size_df_test = basin_size_df_test.append(pd.Series([change_in_attractor_to_chosen_attractor], name="{}".format(constants)))
            # else:
            #     # add to control dataframe
            #     basin_size_df_control = basin_size_df_control.append(pd.Series([change_in_attractor_to_chosen_attractor], name="{}".format(constants)))


            ## EXPERIMENT 3 average (over all states) hamming distance over time

            new_trajectories = determinine_trajectories(new_primes, 
                original_attractors.keys(),
                time_steps=num_time_steps,
                update=update)

            distance = mean_hamming_distance(original_trajectories, 
                    new_trajectories)

            # average over all states 
            distance = distance.mean(axis=0)

            # ## add to correct dataframe 
            if is_test_set:
                # add to test dataframe
                hamming_distance_over_time_df_test = hamming_distance_over_time_df_test.append(\
                    pd.Series(distance, name="{}".format(constants)))
            else:
                # add to control dataframe
                hamming_distance_over_time_df_control = hamming_distance_over_time_df_control.append(\
                    pd.Series(distance, name="{}".format(constants)))
            
            
            # ## EXPERIMENT 4 derrida curve
            # new_derrida = determine_derrida_curve(new_primes, 
            #     states=original_attractors,
            #     update=update,)

            # # # ## add to correct dataframe 
            # if is_test_set:
            #     # add to test dataframe
            #     derrida_curve_df_test = derrida_curve_df_test.append(\
            #         pd.Series(new_derrida, name="{}".format(constants)))
            # else:
            #     # add to control dataframe
            #     derrida_curve_df_control = derrida_curve_df_control.append(\
            #         pd.Series(new_derrida, name="{}".format(constants)))

    # print (change_in_attractor_df_control)
    # print (change_in_attractor_df_test)

    # print (basin_size_df_control)
    # print (basin_size_df_test)

    # print (hamming_distance_over_time_df_control)
    # print (hamming_distance_over_time_df_test)

    # plt.subplot(1, 2, 1)
    # plt.plot(hamming_distance_over_time_df_control.T)
    # plt.title("CONTROL")
    # plt.ylim([0, .5])

    # plt.subplot(1, 2, 2)
    # plt.plot(hamming_distance_over_time_df_test.T)
    # plt.title("TEST")
    # plt.ylim([0, .5])


    # plt.show()

    # print (derrida_curve_df_control)
    # print (derrida_curve_df_test)

    # plt.subplot(1, 2, 1)
    # plt.plot(derrida_curve_df_control.T)
    # plt.title("CONTROL")
    # plt.ylim([-1, 1])

    # plt.subplot(1, 2, 2)
    # plt.plot(derrida_curve_df_test.T)
    # plt.title("TEST")
    # plt.ylim([-1, 1])

    # plt.show()

    # plot_change_in_attractor(
    #     basin_size_df_control,
    #     basin_size_df_test,
    # )

    plot_hamming_distance_over_time(
        hamming_distance_over_time_df_control,
        hamming_distance_over_time_df_test,
    )

    # plot_derrida(
    #     derrida_curve_df_control,
    #     derrida_curve_df_test,
    # )


if __name__ == "__main__":
    main()