#!/bin/bash

dataset=egfr
score_function=module_pos

args=$(echo "--edgelist datasets/${dataset}/edgelist_with_genes.tsv 
--output decomposition_output/${dataset}
--score-function ${score_function}
--draw" )
python src/decomp/hierarchical_scc_decomposition.py $args