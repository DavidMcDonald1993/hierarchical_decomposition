#!/bin/bash

#SBATCH --job-name=expressionsGastric
#SBATCH --output=expressionsGastric_%A_%a.out
#SBATCH --error=expressionsGastric_%A_%a.err
#SBATCH --array=0-5
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

# pip install --user numpy networkx pandas PyBoolNet


output_dir=results/gastric/
output_genes="Caspase8 Caspase9 FOXO RSK TCF cMYC"
chunk_no=${SLURM_ARRAY_TASK_ID}

if [ ! -f ${output_dir}/chunks/Caspase8_expressions_chunk_${chunk_no}.csv ]
then

    args=$(echo "--edgelist datasets/gastric/edgelist.tsv 
    --output ${output_dir}
    --primes datasets/gastric/gastric.json  
    --output_genes ${output_genes}
    -i ${chunk_no}")

    echo ${args}

    python src/evaluation/calculate_expressions.py ${args} 
fi