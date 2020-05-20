#!/bin/bash

#SBATCH --job-name=expressionsSynthetic
#SBATCH --output=expressionsSynthetic_%A_%a.out
#SBATCH --error=expressionsSynthetic_%A_%a.err
#SBATCH --array=0-59
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

# pip install --user numpy networkx pandas PyBoolNet
seeds=({000..29})

n_seeds=${#seeds[@]}
n_chunks=2

seed_id=$((SLURM_ARRAY_TASK_ID / n_chunks % n_seeds))
chunk_no=$((SLURM_ARRAY_TASK_ID % n_chunks))

seed=${seeds[$seed_id]}

output_dir=results/synthetic_bow_tie/${seed}
output_genes="n13 n14 n15"
# chunk_no=${SLURM_ARRAY_TASK_ID}

if [ ! -f ${output_dir}/chunks/n13_expressions_chunk_${chunk_no}.csv ]
then
    args=$(echo "--edgelist datasets/synthetic_bow_tie_networks/${seed}/edgelist.tsv 
    --output ${output_dir}
    --primes datasets/synthetic_bow_tie_networks/${seed}/network.bnet 
    --output_genes ${output_genes}
    -i ${chunk_no}")

    echo $args

    python src/evaluation/calculate_expressions.py ${args} 
fi