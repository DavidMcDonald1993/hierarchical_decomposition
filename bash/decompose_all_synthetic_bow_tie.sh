#!/bin/bash

for seed in {0..4}
do
    d=$(printf synthetic_bow_tie_networks/%03d $seed)
    edgelist=${d}/edgelist.tsv

    python draw_bow_tie.py --edgelist $edgelist --output $d
    python hierarchical_scc_decomposition.py --edgelist $edgelist --output $d --draw
done