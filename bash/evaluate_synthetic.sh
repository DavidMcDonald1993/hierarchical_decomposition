#!/bin/bash 


#SBATCH --job-name=evaluateSynthetic
#SBATCH --output=evaluateSynthetic_%A_%a.out
#SBATCH --error=evaluateSynthetic_%A_%a.err
#SBATCH --array=0-29
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

# pip install --user numpy networkx pandas PyBoolNet
seeds=({000..29})

seed=${seeds[$SLURM_ARRAY_TASK_ID]}

python src/evaluation/evaluate_expressions.py -d results/synthetic_bow_tie/${seed}