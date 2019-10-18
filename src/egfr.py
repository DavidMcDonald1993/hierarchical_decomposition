import numpy as np 
import networkx as nx

import PyBoolNet
from PyBoolNet import Attractors
from PyBoolNet import StateTransitionGraphs as STGs

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

        if i % 1000 == 0:
            print ("processed state", i, "/", len(states))


    return attractors, stg


def main():

    update = "synchronous"

    # primes = PyBoolNet.FileExchange.bnet2primes("datasets/EGFR_full/egfr.bnet", FnamePRIMES="datasets/EGFR_full/egfr_primes.json")
    primes = PyBoolNet.FileExchange.read_primes("datasets/EGFR_full/egfr_primes.json")





    # primes = PyBoolNet.FileExchange.bnet2primes("synthetic_bow_tie_networks/000/network.bnet")

    # stg = STGs.primes2stg(primes, update)


    # states = list(stg.nodes())

    # print (len(states))


    # attractors, stg_computed = build_STG_and_determine_attractors(primes, states)


    # for attractor in set(map(frozenset, attractors)):
    #     print (attractor)

    # raise SystemExit


    core_of_the_core = set(["pi3k", "pip3", "gab1"])

    for gene in core_of_the_core:
        assert gene in primes

    cancer_network =  PyBoolNet.PrimeImplicants.\
        create_constants(primes, 
        {"erbb1": 1}, 
        Copy=True)

    num_state_samples = 10000

    print ("state space is too large --  sampling", 
            num_state_samples, 
            "states")
    states = set()
    while len(states) < num_state_samples:
        state = tuple(np.random.randint(2, size=len(cancer_network)))
        states.add(state)
    states = list(map(lambda state: 
        STGs.state2str({p: s for p, s in zip(cancer_network, state)}), states))

    print ("completed sampling states -- determining attractors")

    # attrs = []

    # for i, state in enumerate(states):

    #     attrs.append(STGs.state2str(Attractors.find_attractor_state_by_randomwalk_and_ctl(cancer_network, 
    #     update, 
    #     InitialState=state) ))
    #     if i % 1 == 0:
    #         print ("done", i ,"/", num_state_samples)

    attractors, stg = build_STG_and_determine_attractors(cancer_network, states)
    
    print ("done", len(stg))

    print ("found", len(set(map(frozenset, attractors))), "unique attractors")




if __name__ == "__main__":
    main()