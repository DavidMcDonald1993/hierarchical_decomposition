#!/bin/bash

#SBATCH --job-name=expressionsEGFR
#SBATCH --output=expressionsEGFR_%A_%a.out
#SBATCH --error=expressionsEGFR_%A_%a.err
#SBATCH --array=0-61
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

# pip install --user numpy networkx pandas PyBoolNet

output_dir=results/EGFR/erbb11
output_genes="elk1 creb ap1 cmyc p70s6_2 hsp27 pro_apoptotic"
chunk_no=${SLURM_ARRAY_TASK_ID}

if [ ! -f ${output_dir}/chunks/ap1_expressions_chunk_${chunk_no}.csv ]
then
    args=$(echo "--edgelist datasets/EGFR_full/edgelist_with_genes.tsv 
    --output ${output_dir}
    --primes datasets/EGFR_full/egfr_primes.json 
    --cancer_mutation erbb11 
    --output_genes $output_genes
    -i ${chunk_no}")

    echo $args

    python src/evaluation/calculate_expressions.py ${args} 
fi