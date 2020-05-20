#!/bin/bash

#SBATCH --job-name=decomposeSynthetic
#SBATCH --output=decomposeSynthetic_%A_%a.out
#SBATCH --error=decomposeSynthetic_%A_%a.err
#SBATCH --array=0-29
#SBATCH --time=05:00:00
#SBATCH --ntasks=1
#SBATCH --mem=1G

module purge
module load bluebear
module load bear-apps/2019a
module load Python/3.7.2-GCCcore-8.2.0

score_function=density

seeds=({0..29})
seed=${seeds[$SLURM_ARRAY_TASK_ID]}

# for seed in {0..29}
# do
d=$(printf datasets/synthetic_bow_tie_networks/%03d $seed)
edgelist=${d}/edgelist.tsv

draw_args=$(echo "--edgelist ${edgelist}
--output ${d}")
python src/drawing/draw_bow_tie.py $draw_args

decomp_args=$(echo "--edgelist ${edgelist}
--output decomposition_output/synthetic_bow_tie/${seed}
--score-function ${score_function}
--draw")
python src/decomp/hierarchical_scc_decomposition.py ${decomp_args}
# done