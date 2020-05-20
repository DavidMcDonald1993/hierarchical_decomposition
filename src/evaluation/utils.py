
import os
import fcntl
import functools
import itertools

import numpy as np
import networkx as nx 
import pandas as pd

from PyBoolNet import StateTransitionGraphs as STGs

def select_states(primes, num_state_samples=10000, seed=0):

	n = len(primes)

	if n <= 16:
		print ("using entire state space")

		states = set(itertools.product([0, 1], repeat=n))

	else:

		print ("sampling", 
				num_state_samples, 
				"states")
		states = set()
		np.random.seed(seed)

		while len(states) < num_state_samples:
			state = tuple(np.random.randint(2, size=n))
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

		assert next_state in stg
		stg.add_edge(state, next_state)


		if state == next_state: # fixed point attractor
			attractors.append([state])
				
		else: # cyclic attractor

			visited = [state]
			while next_state not in visited:
				visited.append(next_state)
				assert len(list(stg.neighbors(next_state))) == 1
				next_state = list(stg.neighbors(next_state))[0]

			idx = visited.index(next_state)
			attractor = visited[idx:]
			attractors.append(attractor)

	return attractors

def compute_average_activation(primes, genes, attractors):

	counts = {gene: [] for gene in genes}

	for attractor in attractors:

		attractor_counts = {gene: 0 for gene in genes}

		for state in attractor:

			state_dict = STGs.state2dict(primes, state)

			for gene in genes:
				attractor_counts[gene] += state_dict[gene]

		for gene in genes:
			counts[gene].append(attractor_counts[gene] / \
				len(attractor))

	return counts



