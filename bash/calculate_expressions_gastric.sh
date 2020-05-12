#!/bin/bash

#SBATCH --job-name=expressionsGastric
#SBATCH --output=expressionsGastric_%A_%a.out
#SBATCH --error=expressionsGastric_%A_%a.err
#SBATCH --array=2-190
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

# pip install --user numpy networkx pandas PyBoolNet

args=$(echo "--edgelist datasets/gastric/edgelist.tsv 
--output results/gastric/ 
--primes datasets/gastric/gastric.json  
--output_genes Caspase8 Caspase9 FOXO RSK TCF cMYC
-i ${SLURM_ARRAY_TASK_ID}")

echo ${args}

python src/evaluation/calculate_expressions.py ${args} 
