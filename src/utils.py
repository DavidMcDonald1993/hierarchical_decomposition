
import os
import fcntl
import functools

import numpy as np
import networkx as nx 
import pandas as pd

from PyBoolNet import StateTransitionGraphs as STGs

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


	return attractors

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




def lock_method(lock_filename):
	''' Use an OS lock such that a method can only be called once at a time. '''

	def decorator(func):

		@functools.wraps(func)
		def lock_and_run_method(*args, **kwargs):

			# Hold program if it is already running 
			# Snippet based on
			# http://linux.byexamples.com/archives/494/how-can-i-avoid-running-a-python-script-multiple-times-implement-file-locking/
			fp = open(lock_filename, 'r+')
			done = False
			while not done:
				try:
					fcntl.lockf(fp, fcntl.LOCK_EX | fcntl.LOCK_NB)
					done = True
				except IOError:
					pass
			return func(*args, **kwargs)

		return lock_and_run_method

	return decorator 

def threadsafe_fn(lock_filename, fn, *args, **kwargs ):
	lock_method(lock_filename)(fn)(*args, **kwargs)

def save_test_results(filename, index, data, ):
	d = pd.DataFrame(index=[index], data=[data])
	if os.path.exists(filename):
		test_df = pd.read_csv(filename, sep=",", index_col=0)
		test_df.columns = [int(col) for col in test_df.columns]
		test_df = d.combine_first(test_df)
	else:
		test_df = d
	test_df.to_csv(filename, sep=",")

def threadsafe_save_test_results(lock_filename, filename, index, data):
	print ("obtained lock", lock_filename, "for file", filename)
	threadsafe_fn(lock_filename, save_test_results, filename=filename, index=index, data=data)
