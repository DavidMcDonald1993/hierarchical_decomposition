#!/bin/bash

#SBATCH --job-name=expressionsGastric
#SBATCH --output=expressionsGastric_%A_%a.out
#SBATCH --error=expressionsGastric_%A_%a.err
#SBATCH --array=0-2145
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5G

# module purge
# module load bluebear
# module load bear-apps/2019a
# module load Python/3.7.2-GCCcore-8.2.0

pip install --user numpy networkx pandas PyBoolNet

args=$(echo "--edgelist datasets/gastric/edgelist.tsv 
--output test-results/gastric/ 
--primes datasets/gastric/gastric.json  
--output_genes Caspase8 Caspase9 FOXO RSK TCF cMYC")

python src/calculate_expressions.py $args -i ${SLURM_ARRAY_TASK_ID}
