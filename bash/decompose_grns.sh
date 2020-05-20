#!/bin/bash

grns=(yeast_cell_cycle mouse_myeloid_development mammalian_cortical_development arabidopsis_thaliana_development)

score_function=module_pos

for grn in ${grns[@]}
do 
    args=$(echo "--edgelist datasets/${grn}/edgelist.tsv 
    --mapping datasets/${grn}/gene_ids.csv 
    --output decomposition_output/${grn}
    --score-function ${score_function}
    --draw" )
    python src/decomp/hierarchical_scc_decomposition.py $args
done 