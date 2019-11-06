#!/bin/bash

#SBATCH --job-name=expressionsEGFR
#SBATCH --output=expressionsEGFR_%A_%a.out
#SBATCH --error=expressionsEGFR_%A_%a.err
#SBATCH --array=0-595
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --mem=5G

args=$(echo "--edgelist datasets/EGFR_full/edgelist_with_genes.tsv 
--output results/EGFR/erbb11/ 
--primes datasets/EGFR_full/egfr_primes.json 
--cancer_mutation erbb11 
--output_genes elk1 creb ap1 cmyc p70s6_2 hsp27 pro_apoptotic")

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

pip install --user numpy networkx pandas
pip install --user PyBoolNet-2.2.6_linux64.tar.gz 

python src/calculate_expressions.py $args -i ${SLURM_ARRAY_TASK_ID}
