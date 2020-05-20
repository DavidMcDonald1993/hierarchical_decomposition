#!/bin/bash

dataset=narang
score_function=module_pos


args=$(echo "--edgelist datasets/${dataset}/edgelist.tsv 
--mapping datasets/${dataset}/gene_ids.csv
--output decomposition_output/${dataset}
--score-function ${score_function}
--draw" )
python src/decomp/hierarchical_scc_decomposition.py $args